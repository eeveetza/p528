function result = tl_p528(d__km, h_1__meter, h_2__meter, f__mhz,  T_pol, p)
% tl_p528 - computes basic transmission loss according to Recommendation
% ITU-R P.528-5 for aeronautical mobile and radionavigation services.
%      result = tl_p528(d__km, h_1__meter, h_2__meter, f__mhz, t_pol, p)
%
%  Description:  Annex 2, Section 3 of Recommendation ITU-R
%                P.528-5, "Propagation curves for aeronautical mobile and
%                radionavigation services using the VHF, UHF and SHF bands"
%
%        Input:  d__km             - Path distance, in km
%                h_1__meter        - Height of the low terminal, in meters
%                h_2__meter        - Height of the high terminal, in meters
%                f__mhz            - Frequency, in MHz
%                p                 - Time percentage
%
%      Outputs:  result            - Result structure containing various
%                                    computed parameters
%
% Translated and adapted to MATLAB/Octave by Ivica Stevanovic (OFCOM CH)
% starting from the original C++ code by William Kozma Jr (NTIA, USA)
%
%
%   C++ Author:  William Kozma Jr
%                wkozma@ntia.gov
%                US Dept of Commerce, NTIA/ITS
%                June 2021 : Geneva Study Group 3 Meetings
%

%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    18AUG21     Ivica Stevanovic, OFCOM         Initial version

% MATLAB Version '9.8.0.1417392 (R2020a) Update 4' used in the development/translation of this code


% Various constants definitions used by the code

Const.epsilon_0 = 8.854187817e-12;         % Vacuum permittivity (F/m)
Const.a_0__km = 6371.0;                    % Earth radius, in km
Const.a_e__km = 9257.0;                    % Effective Earth radius, in km
Const.N_s =                                341;
Const.epsilon_r =                          15.0;
Const.sigma =                              0.005;
Const.LOS_EPSILON =                        0.00001;
Const.THIRD =                              1.0 / 3.0;

Const.CONST_MODE__SEARCH =                 0;
Const.CONST_MODE__DIFFRACTION =            1;
Const.CONST_MODE__SCATTERING =             2;

Const.CASE_1 =                             1;
Const.CASE_2 =                             2;

Const.PROP_MODE__NOT_SET =                 0;
Const.PROP_MODE__LOS =                     1;
Const.PROP_MODE__DIFFRACTION =             2;
Const.PROP_MODE__SCATTERING =              3;

% List of valid polarizations
Const.POLARIZATION__HORIZONTAL =           0;
Const.POLARIZATION__VERTICAL =             1;

Const.Y_pi_99_INDEX =                      16;

% RETURN CODES
Const.SUCCESS =                            0;
Const.ERROR_VALIDATION__D_KM =             1;
Const.ERROR_VALIDATION__H_1 =              2;
Const.ERROR_VALIDATION__H_2 =              3;
Const.ERROR_VALIDATION__TERM_GEO =         4;
Const.ERROR_VALIDATION__F_MHZ_LOW =        5;
Const.ERROR_VALIDATION__F_MHZ_HIGH =       6;
Const.ERROR_VALIDATION__PERCENT_LOW =      7;
Const.ERROR_VALIDATION__PERCENT_HIGH =     8;
Const.ERROR_VALIDATION__POLARIZATION =     9;
Const.ERROR_HEIGHT_AND_DISTANCE =          10;
Const.WARNING__DFRAC_TROPO_REGION =        20;

% Related to P835

Const.RHO_0__M_KG =                        7.5;

Const.ERROR_HEIGHT_TOO_SMALL =             -1;
Const.ERROR_HEIGHT_TOO_LARGE =             -2;



% reset Results struct
result.A_fs__db = 0;
result.A_a__db = 0;
result.A__db = 0;
result.d__km = 0;
result.theta_h1__rad = 0;
result.propagation_mode = Const.PROP_MODE__NOT_SET;

err = ValidateInputs(d__km, h_1__meter, h_2__meter, f__mhz, T_pol, p, Const);

if (err  ~= Const.SUCCESS)
    
    if (err == Const.ERROR_HEIGHT_AND_DISTANCE)
        
        result.A_fs__db = 0;
        result.A_a__db = 0;
        result.A__db = 0;
        result.d__km = 0;
        rtn = Const.SUCCESS;
        result.rtn = rtn;
        return
        
    else
        rtn = err;
        result.rtn = rtn;
        return
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute terminal geometries
%

% Step 1 for low terminal
terminal_1.h_r__km = h_1__meter / 1000;
terminal_1 = TerminalGeometry(f__mhz, terminal_1, Const);

% Step 1 for high terminal
terminal_2.h_r__km = h_2__meter / 1000;
terminal_2 = TerminalGeometry(f__mhz, terminal_2, Const);

%
% Compute terminal geometries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2
path.d_ML__km = terminal_1.d_r__km + terminal_2.d_r__km;                     % [Eqn 3-1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth earth diffraction line calculations
%

% Step 3.1
d_3__km = path.d_ML__km + 0.5 * ( Const.a_e__km^2 / f__mhz ).^Const.THIRD;     % [Eqn 3-2]
d_4__km = path.d_ML__km + 1.5 * ( Const.a_e__km^2 / f__mhz ).^Const.THIRD;     % [Eqn 3-3]

% Step 3.2
A_3__db = SmoothEarthDiffraction(terminal_1.d_r__km, terminal_2.d_r__km, f__mhz, d_3__km, T_pol, Const);
A_4__db = SmoothEarthDiffraction(terminal_1.d_r__km, terminal_2.d_r__km, f__mhz, d_4__km, T_pol, Const);

% Step 3.3
M_d = (A_4__db - A_3__db) / (d_4__km - d_3__km);     % [Eqn 3-4]
A_d0 = A_4__db - M_d * d_4__km;                      % [Eqn 3-5]

% Step 3.4
A_dML__db = (M_d * path.d_ML__km) + A_d0;           % [Eqn 3-6]
path.d_d__km = -(A_d0 / M_d);                       % [Eqn 3-7]

%
% End smooth earth diffraction line calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_LOS = 0;

% Step 4.  If the path is in the Line-of-Sight range, call LOS and then exit
if (path.d_ML__km - d__km > 0.001)
    
    result.propagation_mode = Const.PROP_MODE__LOS;
    [result, los_params, K_LOS] = LineOfSight(path, terminal_1, terminal_2, f__mhz, -A_dML__db, p, d__km, T_pol, Const);
    
    rtn = Const.SUCCESS;
    result.rtn = rtn;
    return
    
else
    
    % get K_LOS
    [result, los_params, K_LOS] = LineOfSight(path, terminal_1, terminal_2, f__mhz, -A_dML__db, p, path.d_ML__km - 1, T_pol, Const);
    
    % Step 6.  Search past horizon to find crossover point between Diffraction and Troposcatter models
    
    [rtn, M_d, A_d0, d_crx__km, CASE] = TranshorizonSearch(path, terminal_1, terminal_2, f__mhz, A_dML__db, M_d, A_d0,Const);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute terrain attenuation, A_T__db
    %
    
    % Step 7.1
    A_d__db = M_d * d__km + A_d0;                    % [Eqn 3-14]
    
    % Step 7.2
    %tropo = Troposcatter(path, terminal_1, terminal_2, d__km, f__mhz, Const);
    tropo = Troposcatter(terminal_1, terminal_2, d__km, f__mhz, Const);
    
    % Step 7.3
    
    if (d__km < d_crx__km)
        
        % always in diffraction if less than d_crx
        A_T__db = A_d__db;
        result.propagation_mode = Const.PROP_MODE__DIFFRACTION;
        
    else
        
        if (CASE == Const.CASE_1)
            
            % select the lower loss mode of propagation
            if (tropo.A_s__db <= A_d__db)
                
                A_T__db = tropo.A_s__db;
                result.propagation_mode = Const.PROP_MODE__SCATTERING;
                
            else
                
                A_T__db = A_d__db;
                result.propagation_mode = Const.PROP_MODE__DIFFRACTION;
                
            end
        else % CASE_2
            
            A_T__db = tropo.A_s__db;
            result.propagation_mode = Const.PROP_MODE__SCATTERING;
            
        end
    end
    
    %
    % Compute terrain attenuation, A_T__db
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%
    % Compute variability
    %
    
    % f_theta_h is unity for transhorizon paths
    f_theta_h = 1;
    
    % compute the 50% and p% of the long-term variability distribution
    
    [Y_e__db, ~] = LongTermVariability(terminal_1.d_r__km, terminal_2.d_r__km, d__km, f__mhz, p, f_theta_h, -A_T__db, Const);
    [Y_e_50__db, ~] = LongTermVariability(terminal_1.d_r__km, terminal_2.d_r__km, d__km, f__mhz, 50, f_theta_h, -A_T__db, Const);
    
    % compute the 50% and p% of the Nakagami-Rice distribution
    ANGLE = 0.02617993878;   % 1.5 deg
    
    if (tropo.theta_s >= ANGLE)        % theta_s > 1.5 deg
        K_t__db = 20;
    elseif (tropo.theta_s <= 0.0)
        K_t__db = K_LOS;
    else
        K_t__db = (tropo.theta_s * (20.0 - K_LOS) / ANGLE) + K_LOS;
    end
    
    Y_pi_50__db = 0.0;       %  zero mean
    Y_pi__db = NakagamiRice(K_t__db, p, Const);
    
    % combine the long-term and Nakagami-Rice distributions
    Y_total__db = CombineDistributions(Y_e_50__db, Y_e__db, Y_pi_50__db, Y_pi__db, p);
    
    %
    % Compute variability
    %%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Atmospheric absorption for transhorizon path
    %
    result_v = SlantPathAttenuation(f__mhz / 1000, 0, tropo.h_v__km, pi / 2, Const);
    
    result.A_a__db = terminal_1.A_a__db + terminal_2.A_a__db + 2 * result_v.A_gas__db;   % [Eqn 3-17]
    
    %
    % Atmospheric absorption for transhorizon path
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute free-space loss
    %
    
    r_fs__km = terminal_1.a__km + terminal_2.a__km + 2 * result_v.a__km;         % [Eqn 3-18]
    result.A_fs__db = 20.0 * log10(f__mhz) + 20.0 * log10(r_fs__km) + 32.45;     % [Eqn 3-19]
    
    %
    % Compute free-space loss
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    result.d__km = d__km;
    result.A__db =  result.A_fs__db  +  result.A_a__db +   A_T__db -  Y_total__db;     % [Eqn 3-20]
    
    result.theta_h1__rad = -terminal_1.theta__rad;
    result.rtn = rtn;
    return;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = ValidateInputs(d__km, h_1__meter, h_2__meter, f__mhz, T_pol, p, Const)
%
%  Description:  Validate the model input values
%
%        Input:  d__km             - Path distance, in km
%                h_1__meter        - Height of the low terminal, in meters
%                h_2__meter        - Height of the high terminal, in meters
%                f__mhz            - Frequency, in MHz
%                T_pol             - Code indicating either polarization
%                                      + 0 : POLARIZATION__HORIZONTAL
%                                      + 1 : POLARIZATION__VERTICAL
%                p	                - Time percentage
%
%      Returns:  SUCCESS, or validation error code
%

if (d__km < 0)
    err = Const.ERROR_VALIDATION__D_KM;
    return
end

if (h_1__meter < 1.5 || h_1__meter > 20000)
    err = Const.ERROR_VALIDATION__H_1;
    return
end

if (h_2__meter < 1.5 || h_2__meter > 20000)
    err = Const.ERROR_VALIDATION__H_2;
    return
end

if (h_1__meter > h_2__meter)
    err = Const.ERROR_VALIDATION__TERM_GEO;
    return
end

if (f__mhz < 100)
    err = Const.ERROR_VALIDATION__F_MHZ_LOW;
    return
end

if (f__mhz > 30000)
    err = Const.ERROR_VALIDATION__F_MHZ_HIGH;
    return
end

if (T_pol ~= Const.POLARIZATION__HORIZONTAL && T_pol ~= Const.POLARIZATION__VERTICAL)
    err = Const.ERROR_VALIDATION__POLARIZATION;
    return
end

if (p < 1)
    err = Const.ERROR_VALIDATION__PERCENT_LOW;
    return
end

if (p > 99)
    err = Const.ERROR_VALIDATION__PERCENT_HIGH;
    return
end

if (h_1__meter == h_2__meter && d__km == 0)
    err = Const.ERROR_HEIGHT_AND_DISTANCE;
    return
end

err = Const.SUCCESS;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terminal = TerminalGeometry(f__mhz, terminal, Const)
%  Description:  This file computes the terminal geometry as described
%                in Annex 2, Section 4 of Recommendation ITU-R P.528-5,
%                "Propagation curves for aeronautical mobile and
%                radionavigation services using the VHF, UHF and SHF bands"
%
%        Input:  f__mhz    - Frequency, in MHz
%
%      Outputs:  terminal  - Structure containing parameters dealing
%                            with the geometry of the terminal
%
%
theta_tx__rad = 0;

result = SlantPathAttenuation(f__mhz / 1000, 0, terminal.h_r__km, pi / 2 - theta_tx__rad, Const);
terminal.theta__rad = pi / 2 - result.angle__rad;
terminal.A_a__db = result.A_gas__db;
terminal.a__km = result.a__km;

% compute arc distance
central_angle = ((pi / 2 - result.angle__rad) - theta_tx__rad + result.bending__rad);            % [Thayer, Equ 2], rearranged
terminal.d_r__km = Const.a_0__km * central_angle;

terminal.phi__rad = terminal.d_r__km / Const.a_e__km;                          % [Eqn 4-1]
terminal.h_e__km = (Const.a_e__km / cos(terminal.phi__rad)) - Const.a_e__km;   % [Eqn 4-2]

terminal.delta_h__km = terminal.h_r__km - terminal.h_e__km;        % [Eqn 4-3]

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A_d__db = SmoothEarthDiffraction(d_1__km, d_2__km, f__mhz, d_0__km, T_pol, Const)
%
%
%  Description:  This file computes the smooth earth diffraction loss
%                as described in Annex 2, Section 10 of
%                Recommendation ITU-R P.528-5, "Propagation curves for
%                aeronautical mobile and radionavigation services using
%                the VHF, UHF and SHF bands"
%
%        Input:  d_1__km   - Horizon distance of terminal 1, in km
%                d_2__km   - Horizon distance of terminal 2, in km
%                a_e__km   - Effective earth radius, in km
%                f__mhz    - Frequency, in MHz
%                d_0__km   - Path length of interest, in km
%                T_pol     - Code indicating either polarization
%                              + 0 : POLARIZATION__HORIZONTAL
%                              + 1 : POLARIZATION__VERTICAL
%
%      Returns:  A_d__db   - Diffraction loss, in dB
%$
s = 18000 * Const.sigma / f__mhz;


if (T_pol == Const.POLARIZATION__HORIZONTAL)
    K = 0.01778 * f__mhz^(-Const.THIRD) * ( (Const.epsilon_r - 1)^2 + s^2 )^(-0.25);
else
    K = 0.01778 * f__mhz^(-Const.THIRD) * ( (Const.epsilon_r^2 + s^2 ) / ( (Const.epsilon_r - 1)^2  + s^2 )^0.5  )^0.5;
end

B_0 = 1.607;

% [Vogler 1964, Equ 2] with C_0 = 1 due to "4/3" Earth assumption
x_0__km = (B_0 - K) * (f__mhz^Const.THIRD) * d_0__km;
x_1__km = (B_0 - K) * (f__mhz^Const.THIRD) * d_1__km;
x_2__km = (B_0 - K) * (f__mhz^Const.THIRD) * d_2__km;

% Compute the distance function for the path
G_x__db = DistanceFunction(x_0__km);

% Compute the height functions for the two terminals
F_x1__db = HeightFunction(x_1__km, K);
F_x2__db = HeightFunction(x_2__km, K);

% [Vogler 1964, Equ 1] with C_1(K, b^0) = 20, which is the approximate value for all K (see Figure 5)
A_d__db = G_x__db - F_x1__db - F_x2__db - 20.0;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G_x__db = DistanceFunction(x__km)

% [Vogler 1964, Equ 13]

G_x__db = 0.05751 * x__km - 10.0 * log10(x__km);

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_x__db = HeightFunction(x__km, K)

% [FAA-ES-83-3, Equ 73]
y__db = 40.0 * log10(x__km) - 117.0;

% [Vogler 1964, Equ 13]
G_x__db = DistanceFunction(x__km);

if (x__km <= 200.0)
    x_t__km = 450 / (-(log10(K))^3);       % [Eqn 109]
    
    % [Eqn 110]
    if (x__km >= x_t__km)
        
        if (abs(y__db) < 117)
            F_x__db = y__db;
        else
            F_x__db = -117;
        end
        
    else
        F_x__db = 20 * log10(K) - 15 + (0.000025 * (x__km^2.0) / K);
    end
elseif (x__km > 2000.0)
    
    % [Vogler 1964] F_x ~= G_x for large x (see Figure 7)
    F_x__db = G_x__db;
    
else % Blend y__db with G_x__db for 200 < x__km < 2000
    
    % [FAA-ES-83-3, Equ 72] weighting variable
    W = 0.0134 * x__km * exp(-0.005 * x__km);
    
    % [FAA-ES-83-3, Equ 75]
    F_x__db = W * y__db + (1.0 - W) * G_x__db;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, los_params, K_LOS] = LineOfSight(path, terminal_1, terminal_2, f__mhz, A_dML__db, p, d__km, T_pol, Const)
%
%  Description:  This function computes the total loss in the line-of-sight
%                region as described in Annex 2, Section 6 of
%                Recommendation ITU-R P.528-5, "Propagation curves for
%                aeronautical mobile and radionavigation services using
%                the VHF, UHF and SHF bands"
%
%        Input:  path          - Struct containing path parameters
%                terminal_1    - Struct containing low terminal parameters
%                terminal_2    - Struct containing high terminal parameters
%                f__mhz        - Frequency, in MHz
%                A_dML__db     - Diffraction loss at d_ML, in dB
%                p             - Time percentage
%                d__km         - Path length, in km
%                T_pol         - Code indicating either polarization
%                                  + 0 : POLARIZATION__HORIZONTAL
%                                  + 1 : POLARIZATION__VERTICAL
%
%      Outputs:  los_params    - Struct containing LOS parameters
%                result        - Struct containing P.528 results
%                K_LOS         - K-value
%


% 0.2997925 = speed of light, gigameters per sec
lambda__km = 0.2997925 / f__mhz;                             % [Eqn 6-1]
terminate = lambda__km / 1e6;

% determine psi_limit, where you switch from free space to 2-ray model
% lambda / 2 is the start of the lobe closest to d_ML
psi_limit = FindPsiAtDeltaR(lambda__km / 2, path, terminal_1, terminal_2, terminate, Const);

% "[d_y6__km] is the largest distance at which a free-space value is obtained in a two-ray model
%   of reflection from a smooth earth with a reflection coefficient of -1" [ES-83-3, page 44]
d_y6__km = FindDistanceAtDeltaR(lambda__km / 6, path, terminal_1, terminal_2, terminate, Const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine d_0__km distance
%

% In IF-73, the values for d_0 (d_d in IF-77) were found to be too small when both antennas are low,
% so this "heuristic" was developed to fix that
% [Eqns 8-2 and 8-3]
if (terminal_1.d_r__km >= path.d_d__km || path.d_d__km >= path.d_ML__km)
    
    if (terminal_1.d_r__km > d_y6__km || d_y6__km > path.d_ML__km)
        path.d_0__km = terminal_1.d_r__km;
    else
        path.d_0__km = d_y6__km;
    end
    
elseif (path.d_d__km < d_y6__km && d_y6__km < path.d_ML__km)
    path.d_0__km = d_y6__km;
else
    path.d_0__km = path.d_d__km;
    
end
% Determine d_0__km distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Tune d_0__km distance
%

% Now that we have d_0, lets carefully walk it forward, 1 meter at a time, to tune it to as precise as possible without
%      going beyond the LOS region (ie, beyond d_ML)
d_temp__km = path.d_0__km;

los_result = LineOfSightParams();

while (1)
    
    psi = FindPsiAtDistance(d_temp__km, path, terminal_1, terminal_2, Const);
    
    los_result = RayOptics(terminal_1, terminal_2, psi, los_result, Const);
    
    % if the resulting distance is beyond d_0 OR if we incremented again we'd be outside of LOS...
    if (los_result.d__km >= path.d_0__km || (d_temp__km + 0.001) >= path.d_ML__km)
        
        % use the resulting distance as d_0
        path.d_0__km = los_result.d__km;
        break;
    end
    
    d_temp__km = d_temp__km + 0.001;
end

%
% Tune d_0__km distance
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute loss at d_0__km
%

psi_d0 = FindPsiAtDistance(path.d_0__km, path, terminal_1, terminal_2, Const);

los_params = LineOfSightParams();

los_params = RayOptics(terminal_1, terminal_2, psi_d0, los_params, Const);

[los_params, R_Tg] = GetPathLoss(psi_d0, path, f__mhz, psi_limit, A_dML__db, 0, T_pol, los_params, Const);

%
% Compute loss at d_0__km
%%%%%%%%%%%%%%%%%%%%%%%%%%

% tune psi for the desired distance
psi = FindPsiAtDistance(d__km, path, terminal_1, terminal_2, Const);

los_params = RayOptics(terminal_1, terminal_2, psi, los_params, Const);

[los_params, R_Tg] = GetPathLoss(psi, path, f__mhz, psi_limit, A_dML__db, los_params.A_LOS__db, T_pol, los_params, Const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute atmospheric absorption
%

result_slant = SlantPathAttenuation(f__mhz / 1000, terminal_1.h_r__km, terminal_2.h_r__km, pi / 2 - los_params.theta_h1__rad, Const);

result.A_a__db = result_slant.A_gas__db;

%
% Compute atmospheric absorption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute free-space loss
%

result.A_fs__db = 20.0 * log10(los_params.r_0__km) + 20.0 * log10(f__mhz) + 32.45; % [Eqn 6-4]

%
% Compute free-space loss
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% Compute variability
%

% [Eqn 13-1]

if (los_params.theta_h1__rad <= 0.0)
    f_theta_h = 1.0;
elseif (los_params.theta_h1__rad >= 1.0)
    f_theta_h = 0.0;
else
    f_theta_h = max(0.5 - (1 / pi) * (atan(20.0 * log10(32.0 * los_params.theta_h1__rad))), 0);
end

[Y_e__db, A_Y] = LongTermVariability(terminal_1.d_r__km, terminal_2.d_r__km, d__km, f__mhz, p, f_theta_h, los_params.A_LOS__db, Const);
[Y_e_50__db, A_Y] = LongTermVariability(terminal_1.d_r__km, terminal_2.d_r__km, d__km, f__mhz, 50, f_theta_h, los_params.A_LOS__db, Const);

% [Eqn 13-2]

if (A_Y <= 0.0)
    F_AY = 1.0;
elseif (A_Y >= 9.0)
    F_AY = 0.1;
else
    F_AY = (1.1 + (0.9 * cos((A_Y / 9.0) * pi))) / 2.0;
end
% [Eqn 175]

if (los_params.delta_r__km >= (lambda__km / 2.0))
    F_delta_r = 1.0;
elseif (los_params.delta_r__km <= lambda__km / 6.0)
    F_delta_r = 0.1;
else
    F_delta_r = 0.5 * (1.1 - (0.9 * cos(((3.0 * pi) / lambda__km) * (los_params.delta_r__km - (lambda__km / 6.0)))));
end
R_s = R_Tg * F_delta_r * F_AY;       % [Eqn 13-4]

Y_pi_99__db = 10.0 * log10(f__mhz * (result_slant.a__km^3)) - 84.26;	% [Eqn 13-5]
K_t = FindKForYpiAt99Percent(Y_pi_99__db, Const);

W_a = 10.0^(K_t / 10.0);         % [Eqn 13-6]
W_R = R_s^2 + 0.01^2;            % [Eqn 13-7]
W = W_R + W_a;                   % [Eqn 13-8]

% [Eqn 13-9]
if (W <= 0.0)
    K_LOS = -40.0;
else
    
    K_LOS = 10.0 * log10(W);
end

if (K_LOS < -40.0)
    K_LOS = -40.0;
end

Y_pi_50__db = 0.0;   %  zero mean
Y_pi__db = NakagamiRice(K_LOS, p, Const);

Y_total__db = -CombineDistributions(Y_e_50__db, Y_e__db, Y_pi_50__db, Y_pi__db, p);

%
% Compute variability
%%%%%%%%%%%%%%%%%%%%%

result.d__km = los_params.d__km;
result.A__db = result.A_fs__db + result.A_a__db - los_params.A_LOS__db + Y_total__db;
result.theta_h1__rad = los_params.theta_h1__rad;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = FindPsiAtDistance(d__km, path, terminal_1, terminal_2, Const)

if (d__km == 0)
    psi = pi/2;
    return
end

% initialize to start at mid-point
psi = pi / 2;
delta_psi = -pi / 4;


while (1)
    
    psi = psi + delta_psi; % new psi
    params_temp = LineOfSightParams();
    
    params_temp = RayOptics(terminal_1, terminal_2, psi, params_temp, Const);
    
    d_psi__km = params_temp.d__km;
    
    % compute delta
    if (d_psi__km > d__km)
        delta_psi = abs(delta_psi) / 2;
    else
        delta_psi = -abs(delta_psi) / 2;
    end
    if (abs(d__km - d_psi__km) <= 1e-3 || (abs(delta_psi) <= 1e-12))  % get within 1 meter of desired delta_r value
        break
    end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = FindPsiAtDeltaR(delta_r__km, path, terminal_1, terminal_2, terminate, Const)

psi = pi / 2;
delta_psi = -pi / 4;

while(1)
    
    psi = psi + delta_psi;
    
    params_temp = LineOfSightParams();
    
    params_temp = RayOptics(terminal_1, terminal_2, psi, params_temp, Const);
    
    if (params_temp.delta_r__km > delta_r__km)
        delta_psi = -abs(delta_psi) / 2;
    else
        delta_psi = abs(delta_psi) / 2;
    end
    
    if (abs(params_temp.delta_r__km - delta_r__km) <= terminate)
        break
    end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d__km = FindDistanceAtDeltaR(delta_r__km, path, terminal_1, terminal_2, terminate, Const)

psi = pi / 2;
delta_psi = -pi / 4;


while(1)
    
    psi = psi + delta_psi;
    
    params_temp = LineOfSightParams();
    
    params_temp = RayOptics(terminal_1, terminal_2, psi, params_temp, Const);
    
    if (params_temp.delta_r__km > delta_r__km)
        delta_psi = -abs(delta_psi) / 2;
    else
        delta_psi = abs(delta_psi) / 2;
    end
    
    if (abs(params_temp.delta_r__km - delta_r__km) <= terminate)
        break
    end
end

d__km = params_temp.d__km;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = RayOptics(terminal_1, terminal_2, psi, params, Const)

%
%  Description:  This function computes the line-of-sight ray optics
%                as described in Annex 2, Section 7 of
%                Recommendation ITU-R P.528-5, "Propagation curves for
%                aeronautical mobile and radionavigation services using
%                the VHF, UHF and SHF bands"
%
%        Input:  terminal_1    - Structure holding low terminal parameters
%                terminal_2    - Structure holding high terminal parameters
%                psi           - Reflection angle, in radians
%
%      Outputs:  params        - Structure holding resulting parameters
%


z = (Const.a_0__km / Const.a_e__km) - 1;       % [Eqn 7-1]
k_a = 1 / (1 + z * cos(psi));                  % [Eqn 7-2]
params.a_a__km = Const.a_0__km * k_a;          % [Eqn 7-3]

delta_h_a1__km = terminal_1.delta_h__km * (params.a_a__km - Const.a_0__km) / (Const.a_e__km - Const.a_0__km);  % [Eqn 7-4]
delta_h_a2__km = terminal_2.delta_h__km * (params.a_a__km - Const.a_0__km) / (Const.a_e__km - Const.a_0__km);  % [Eqn 7-4]

H__km = [0 0];
H__km(1) = terminal_1.h_r__km - delta_h_a1__km;    % [Eqn 7-5]
H__km(2) = terminal_2.h_r__km - delta_h_a2__km;    % [Eqn 7-5]

Hprime__km = [0 0];
for i = 1:2
    
    params.z__km(i) = params.a_a__km + H__km(i);                                  % [Eqn 7-6]
    params.theta(i) = acos(params.a_a__km * cos(psi) / params.z__km(i)) - psi;    % [Eqn 7-7]
    params.D__km(i) = params.z__km(i) * sin(params.theta(i));                     % [Eqn 7-8]
    
    % [Eqn 7-9]
    if (psi > 1.56)
        Hprime__km(i) = H__km(i);
    else
        Hprime__km(i) = params.D__km(i) * tan(psi);
    end
end

delta_z = abs(params.z__km(1) - params.z__km(2));   % [Eqn 7-10]

params.d__km = max(params.a_a__km * (params.theta(1) + params.theta(2)), 0);  % [Eqn 7-11]

alpha = atan((Hprime__km(2) - Hprime__km(1)) / (params.D__km(1) + params.D__km(2)));  % [Eqn 7-12]
params.r_0__km = max(delta_z, (params.D__km(1) + params.D__km(2)) / cos(alpha));            % [Eqn 7-13]
params.r_12__km = (params.D__km(1) + params.D__km(2)) / cos(psi);                           % [Eqn 7-14]

params.delta_r__km = 4.0 * Hprime__km(1) * Hprime__km(2) / (params.r_0__km + params.r_12__km);  % [Eqn 7-15]

params.theta_h1__rad = alpha - params.theta(1);                % [Eqn 7-16]
params.theta_h2__rad = -(alpha + params.theta(2));             % [Eqn 7-17]
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params, R_Tg] = GetPathLoss(psi__rad, path, f__mhz, psi_limit, A_dML__db, A_d_0__db, T_pol, params, Const)

%
%  Description:  This function computes the line of sight loss
%                as described in Annex 2, Section 8 of
%                Recommendation ITU-R P.528-5, "Propagation curves for
%                aeronautical mobile and radionavigation services using
%                the VHF, UHF and SHF bands"
%
%        Input:  psi__rad      - Reflection angle, in rad
%                path          - Struct containing path parameters
%                f__mhz        - Frequency, in MHz
%                psi_limit     - Angular limit separating FS and 2-Ray, in rad
%                A_dML__db     - Diffraction loss at d_ML, in dB
%                A_d_0__db     - Loss at d_0, in dB
%                T_pol         - Code indicating either polarization
%                                  + 0 : POLARIZATION__HORIZONTAL
%                                  + 1 : POLARIZATION__VERTICAL
%
%      Outputs:  params        - Line of sight loss params
%                R_Tg          - Reflection parameter


[R_g, phi_g] = ReflectionCoefficients(psi__rad, f__mhz, T_pol, Const);

if (tan(psi__rad) >= 0.1)
    D_v = 1.0;
else
    
    r_1 = params.D__km(1) / cos(psi__rad);       % [Eqn 8-3]
    r_2 = params.D__km(2) / cos(psi__rad);       % [Eqn 8-3]
    R_r = (r_1 * r_2) / params.r_12__km;         % [Eqn 8-4]
    
    term_1 = (2 * R_r * (1 + (sin(psi__rad))^2)) / (params.a_a__km * sin(psi__rad));
    term_2 = (2 * R_r / params.a_a__km)^2;
    D_v = (1.0 + term_1 + term_2)^(-0.5);         % [Eqn 8-5]
end

% Ray-length factor, [Eqn 8-6]
F_r = min(params.r_0__km / params.r_12__km, 1);

R_Tg = R_g * D_v * F_r;                            % [Eqn 8-7]

params.A_LOS__db = 0;

if (params.d__km > path.d_0__km)
    
    % [Eqn 8-1]
    params.A_LOS__db = ((params.d__km - path.d_0__km) * (A_dML__db - A_d_0__db) / (path.d_ML__km - path.d_0__km)) + A_d_0__db;
    
else
    
    lambda__km = 0.2997925 / f__mhz;	% [Eqn 8-2]
    
    if (psi__rad > psi_limit)
        
        % ignore the phase lag; Step 8-2
        params.A_LOS__db = 0;
        
    else
        
        % Total phase lag of the ground reflected ray relative to the direct ray
        
        % [Eqn 8-8]
        phi_Tg = (2 * pi * params.delta_r__km / lambda__km) + phi_g;
        
        % [Eqn 8-9]
        cplx = R_Tg * cos(phi_Tg) -1j*R_Tg * sin(phi_Tg);
        
        % [Eqn 8-10]
        W_RL = min(abs(1.0 + cplx), 1.0);
        
        % [Eqn 8-11]
        W_R0 = W_RL^2;
        
        % [Eqn 8-12]
        params.A_LOS__db = 10.0 * log10(W_R0);
    end
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_g, phi_g] = ReflectionCoefficients(psi__rad, f__mhz, T_pol, Const)
%
%  Description:  This function computes the reflection coefficients
%                as described in Annex 2, Section 9 of
%                Recommendation ITU-R P.528-5, "Propagation curves for
%                aeronautical mobile and radionavigation services using
%                the VHF, UHF and SHF bands"
%
%        Input:  psi__rad  - Reflection angle, in rad
%                f__mhz    - Frequency, in MHz
%                T_pol     - Code indicating either polarization
%                              + 0 : POLARIZATION__HORIZONTAL
%                              + 1 : POLARIZATION__VERTICAL
%
%      Outputs:  R_g       - Real part
%                phi_g     - Imaginary part


if (psi__rad <= 0.0)
    
    psi__rad = 0.0;
    sin_psi = 0.0;
    cos_psi = 1.0;
    
elseif (psi__rad >= pi / 2)
    
    psi__rad = pi / 2;
    sin_psi = 1.0;
    cos_psi = 0.0;
    
else
    
    sin_psi = sin(psi__rad);
    cos_psi = cos(psi__rad);
end

X = (18000.0 * Const.sigma) / f__mhz;              % [Eqn 9-1]
Y = Const.epsilon_r - (cos_psi)^2;                 % [Eqn 9-2]
T = sqrt(Y^2 + X^2) + Y;                           % [Eqn 9-3]
P = sqrt(T * 0.5);                                 % [Eqn 9-4]
Q = X / (2.0 * P);                                 % [Eqn 9-5]

% [Eqn 9-6]

if (T_pol == Const.POLARIZATION__HORIZONTAL)
    B = 1.0 / (P^2 + Q^2);
else
    B = ( (Const.epsilon_r)^2 + X^2 ) / (P^2 + Q^2);
end
% [Eqn 9-7]

if (T_pol == Const.POLARIZATION__HORIZONTAL)
    A = (2.0 * P) / (P^2 + Q^2);
else
    A = (2.0 * (P * Const.epsilon_r + Q * X)) / (P^2 + Q^2);
end

% [Eqn 9-8]
R_g = sqrt( (1.0 + ( B * sin_psi^2 ) - (A * sin_psi)) / (1.0 + (B * sin_psi^2) + (A * sin_psi)));

% [Eqn 9-9]

if (T_pol == Const.POLARIZATION__HORIZONTAL)
    alpha = atan2(-Q, sin_psi - P);
else
    alpha = atan2((Const.epsilon_r * sin_psi) - Q, Const.epsilon_r * sin_psi - P);
end

% [Eqn 9-10]

if (T_pol == Const.POLARIZATION__HORIZONTAL)
    beta = atan2(Q, sin_psi + P);
else
    beta = atan2((X * sin_psi) + Q, Const.epsilon_r * sin_psi + P);
end

% [Eqn 9-11]
phi_g = alpha - beta;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_e__db, A_Y] = LongTermVariability(d_r1__km, d_r2__km, d__km, f__mhz, p, f_theta_h, A_T, Const)

d_qs__km = 65.0 * (100.0 / f__mhz)^Const.THIRD;              % [Eqn 14-1]
d_Lq__km = d_r1__km + d_r2__km;                              % [Eqn 14-2]
d_q__km = d_Lq__km + d_qs__km;                               % [Eqn 14-3]

% [Eqn 14-4]

if (d__km <= d_q__km)
    d_e__km = (130.0 * d__km) / d_q__km;
else
    d_e__km = 130.0 + d__km - d_q__km;
end

% [Eqns 14-5 and 14-6]

if (f__mhz > 1600.0)
    
    g_10 = 1.05;
    g_90 = 1.05;
    
else
    
    g_10 = (0.21 * sin(5.22 * log10(f__mhz / 200.0))) + 1.28;
    g_90 = (0.18 * sin(5.22 * log10(f__mhz / 200.0))) + 1.23;
end

% Data Source for Below Consts: Tech Note 101, Vol 2
% Column 1: Table III.4, Row A* (Page III-50)
% Column 2: Table III.3, Row A* (Page III-49)
% Column 3: Table III.5, Row Continental Temperate (Page III-51)

c_1 = [ 2.93e-4, 5.25e-4, 1.59e-5 ];
c_2 = [ 3.78e-8, 1.57e-6, 1.56e-11 ];
c_3 = [ 1.02e-7, 4.70e-7, 2.77e-8 ];

n_1 = [ 2.00, 1.97, 2.32 ];
n_2 = [ 2.88, 2.31, 4.08 ];
n_3 = [ 3.15, 2.90, 3.25 ];

f_inf = [ 3.2, 5.4, 0.0 ];
f_m = [ 8.2, 10.0, 3.9 ];

% [Y_0(90) Y_0(10) V(50)]

Z__db = [0  0  0 ];

for i = 1:3
    
    f_2 = f_inf(i) + ( (f_m(i) - f_inf(i)) * exp( -c_2(i) * (d_e__km^n_2(i)) ) );
    
    Z__db(i) = ( c_1(i) * (d_e__km^n_1(i)) - f_2 ) * exp(-c_3(i) * (d_e__km^n_3(i))) + f_2;
end


if (p == 50)
    Y_p__db = Z__db(3);
elseif (p > 50)
    
    z_90 = InverseComplementaryCumulativeDistributionFunction(90.0 / 100.0);
    z_p = InverseComplementaryCumulativeDistributionFunction(p / 100.0);
    c_p = z_p / z_90;
    
    Y = c_p * (-Z__db(1) * g_90);
    Y_p__db = Y + Z__db(3);
    
else
    
    
    if (p >= 10)
        
        z_10 = InverseComplementaryCumulativeDistributionFunction(10.0 / 100.0);
        z_p = InverseComplementaryCumulativeDistributionFunction(p / 100.0);
        c_p = z_p / z_10;
        
    else
        
        % Source for values p < 10: [15], Table 10, Page 34, Climate 6
        ps = [ 1, 2, 5, 10 ];
        c_ps = [ 1.9507, 1.7166, 1.3265, 1.0000 ];
        
        %auto upper = upper_bound(data::P.begin(), data::P.end(), p);
        %auto dist = distance(data::P.begin(), upper);
        dist = distance_upper(dataP(), p);
        c_p = LinearInterpolation(ps(dist - 1), c_ps(dist - 1), ps(dist), c_ps(dist), p);
    end
    
    Y = c_p * (Z__db(2) * g_10);
    Y_p__db = Y + Z__db(3);
end

Y_10__db = (Z__db(2) * g_10) + Z__db(3);       % [Eqn 14-20]
Y_eI__db = f_theta_h * Y_p__db;                % [Eqn 14-21]
Y_eI_10__db = f_theta_h * Y_10__db;            % [Eqn 14-22]

% A_Y "is used to prevent available signal powers from exceeding levels expected for free-space propagation by an unrealistic
%      amount when the variability about L_b(50) is large and L_b(50) is near its free-space level" [ES-83-3, p3-4]

A_YI = (A_T + Y_eI_10__db) - 3.0;              % [Eqn 14-23]
A_Y = max(A_YI, 0);                            % [Eqn 14-24]
Y_e__db = Y_eI__db - A_Y;                      % [Eqn 14-25]

% For percentages less than 10%, do a correction check to,
%    "prevent available signal powers from exceeding levels expected from free-space levels
%     by unrealistic amounts" [Gierhart 1970]
if (p < 10)
    
    c_Y = [ -5.0, -4.5, -3.7, 0.0 ];
    
    %auto upper = upper_bound(data::P.begin(), data::P.end(), p);
    %auto dist = distance(data::P.begin(), upper);
    getP = dataP();
    dist = distance_upper(getP, p);
    
    c_Yi = LinearInterpolation(getP(dist - 1), c_Y(dist - 1), getP(dist), c_Y(dist), p);
    
    Y_e__db = Y_e__db + A_T;
    
    if (Y_e__db > -c_Yi)
        Y_e__db = -c_Yi;
    end
    
    Y_e__db = Y_e__db - A_T;
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = LinearInterpolation(x1, y1, x2, y2, x)
%  Description:  This function performs linear interpolation between the
%                points (x1, y1) and (x2, y2).
%
%       Inputs:  (x1, y1)  - Point 1
%                (x2, y2)  - Point 2
%                x         - Value of the dependent variable
%
%      Returns:  y         - Linearly interpolated value
%


y = (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q_q = InverseComplementaryCumulativeDistributionFunction(q)
%
%  Description:  This function computes the inverse complementary
%                cumulative distribution function approximation as
%               described in Recommendation ITU-R P.1057.  This
%                approximation is sourced from Formula 26.2.23 in
%                Abramowitz & Stegun.  This approximation has an error
%                of abs(epsilon(p)) < 4.5e-4
%
%        Input:  q     - Probability, 0.0 < q < 1.0
%
%      Returns:  Q_q   - Q(q)^-1
%

C_0 = 2.515516;
C_1 = 0.802853;
C_2 = 0.010328;
D_1 = 1.432788;
D_2 = 0.189269;
D_3 = 0.001308;

x = q;
if (q > 0.5)
    x = 1.0 - x;
end
T_x = sqrt(-2.0 * log(x));

zeta_x = ((C_2 * T_x + C_1) * T_x + C_0) / (((D_3 * T_x + D_2) * T_x + D_1) * T_x + 1.0);

Q_q = T_x - zeta_x;

if (q > 0.5)
    Q_q = -Q_q;
end
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nkr = dataNakagamiRiceCurves()

nkr = [ ...
    % K = -40 distribution
    [ ...
    -0.1417,   -0.1252,   -0.1004,   -0.0784,   -0.0634, ...
    -0.0515,   -0.0321,   -0.0155,    0.0000,    0.0156,    0.0323, ...
    0.0518,    0.0639,    0.0791,    0.1016,    0.1271,    0.1441 ...
    ]; ...
    [ ...
    -0.7676,   -0.6811,   -0.5497,   -0.4312,   -0.3504, ...
    -0.2856,   -0.1790,   -0.0870,    0.0000,    0.0878,    0.1828, ...
    0.2953,    0.3651,    0.4537,    0.5868,    0.7390,    0.8420 ...
    ]; ...
    [ ...
    -1.3183,   -1.1738,   -0.9524,   -0.7508,   -0.6121, ...
    -0.5003,   -0.3151,   -0.1537,    0.0000,    0.1564,    0.3269, ...
    0.5308,    0.6585,    0.8218,    1.0696,    1.3572,    1.5544 ...
    ]; ...
    [ ...
    -1.6263,   -1.4507,   -1.1805,   -0.9332,   -0.7623, ...
    -0.6240,   -0.3940,   -0.1926,    0.0000,    0.1969,    0.4127, ...
    0.6722,    0.8355,    1.0453,    1.3660,    1.7417,    2.0014 ...
    ]; ...
    [ ...
    -1.9963,   -1.7847,   -1.4573,   -1.1557,   -0.9462, ...
    -0.7760,   -0.4916,   -0.2410,    0.0000,    0.2478,    0.5209, ...
    0.8519,    1.0615,    1.3326,    1.7506,    2.2463,    2.5931 ...
    ]; ...
    [ ...
    -2.4355,   -2.1829,   -1.7896,   -1.4247,   -1.1695, ...
    -0.9613,   -0.6113,   -0.3007,    0.0000,    0.3114,    0.6573, ...
    1.0802,    1.3505,    1.7028,    2.2526,    2.9156,    3.3872 ...
    ]; ...
    [ ...
    -2.9491,   -2.6507,   -2.1831,   -1.7455,   -1.4375, ...
    -1.1846,   -0.7567,   -0.3737,    0.0000,    0.3903,    0.8281, ...
    1.3698,    1.7198,    2.1808,    2.9119,    3.8143,    4.4714 ...
    ];  ...
    [ ...
    -3.5384,   -3.1902,   -2.6407,   -2.1218,   -1.7535, ...
    -1.4495,   -0.9307,   -0.4619,    0.0000,    0.4874,    1.0404, ...
    1.7348,    2.1898,    2.7975,    3.7820,    5.0373,    5.9833 ...
    ]; ...
    [ ...
    -4.1980,   -3.7974,   -3.1602,   -2.5528,   -2.1180, ...
    -1.7565,   -1.1345,   -0.5662,    0.0000,    0.6045,    1.2999, ...
    2.1887,    2.7814,    3.5868,    4.9288,    6.7171,    8.1319 ...
    ]; ...
    [ ...
    -4.9132,   -4.4591,   -3.7313,   -3.0306,   -2.5247, ...
    -2.1011,   -1.3655,   -0.6855,    0.0000,    0.7415,    1.6078, ...
    2.7374,    3.5059,    4.5714,    6.4060,    8.9732,   11.0973 ...
    ]; ...
    [ ...
    -5.6559,   -5.1494,   -4.3315,   -3.5366,   -2.9578, ...
    -2.4699,   -1.6150,   -0.8154,    0.0000,    0.8935,    1.9530, ...
    3.3611,    4.3363,    5.7101,    8.1216,   11.5185,   14.2546 ...
    ]; ...
    [ ...
    -6.3810,   -5.8252,   -4.9219,   -4.0366,   -3.3871, ...
    -2.8364,   -1.8638,   -0.9455,    0.0000,    1.0458,    2.2979, ...
    3.9771,    5.1450,    6.7874,    9.6276,   13.4690,   16.4251 ...
    ]; ...
    [ ...
    -7.0247,   -6.4249,   -5.4449,   -4.4782,   -3.7652, ...
    -3.1580,   -2.0804,   -1.0574,    0.0000,    1.1723,    2.5755, ...
    4.4471,    5.7363,    7.5266,   10.5553,   14.5401,   17.5511 ...
    ]; ...
    [ ...
    -7.5229,   -6.8862,   -5.8424,   -4.8090,   -4.0446, ...
    -3.3927,   -2.2344,   -1.1347,    0.0000,    1.2535,    2.7446, ...
    4.7144,    6.0581,    7.9073,   11.0003,   15.0270,   18.0526 ...
    ]; ...
    [ ...
    -7.8532,   -7.1880,   -6.0963,   -5.0145,   -4.2145, ...
    -3.5325,   -2.3227,   -1.1774,    0.0000,    1.2948,    2.8268, ...
    4.8377,    6.2021,    8.0724,   11.1869,   15.2265,   18.2566 ...
    ]; ...
    [ ...
    -8.0435,   -7.3588,   -6.2354,   -5.1234,   -4.3022, ...
    -3.6032,   -2.3656,   -1.1975,    0.0000,    1.3130,    2.8619, ...
    4.8888,    6.2610,    8.1388,   11.2607,   15.3047,   18.3361 ...
    ]; ...
    [ ...
    -8.2238,   -7.5154,   -6.3565,   -5.2137,   -4.3726, ...
    -3.6584,   -2.3979,   -1.2121,    0.0000,    1.3255,    2.8855, ...
    4.9224,    6.2992,    8.1814,   11.3076,   15.3541,   18.3864 ...
    ]; ...
    ];
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = dataK()

K = [-40, -25, -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 20];
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = dataP()

P = [1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 98, 99 ];
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = distance_upper(vector, value)
% Returns an index of the first element in the vector that is greater than value,
% or of the last element if no such element is found.
k = find(vector > value);
if (isempty(k))
    dist = length(vector);
else
    dist = k(1);
end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = distance_lower(vector, value)
% Returns an index of the first element in the vector that is greater or equal than value,
% or of the last element if no such element is found.

k = find(vector >= value);
if (isempty(k))
    dist = length(vector);
else
    dist = k(1);
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = FindKForYpiAt99Percent(Y_pi_99__db, Const)
%  Description:  This function returns the K-value of the Nakagami-Rice
%                distribution for the given value of Y_pi(99)
%
%        Input:  Y_pi_99__db   - Y_pi(99), in dB
%
%       Returns: K             - K-value
%


getNakagamiRiceCurves = dataNakagamiRiceCurves();
getK = dataK();

% is Y_pi_99__db smaller than the smallest value in the distribution data
%if (Y_pi_99__db < data::NakagamiRiceCurves.front()[Y_pi_99_INDEX])
if ( Y_pi_99__db < getNakagamiRiceCurves(1,Const.Y_pi_99_INDEX+1) )
    %return data::K.front();
    getK = dataK();
    K = getK(1);
    return
end

% search the distribution data and interpolate to find K (dependent variable)
for i = 1:length(getK)
    if (Y_pi_99__db - getNakagamiRiceCurves(i,Const.Y_pi_99_INDEX+1) < 0)
        %return (data::K[i] * (Y_pi_99__db - data::NakagamiRiceCurves[i - 1][Y_pi_99_INDEX]) - data::K[i - 1] * (Y_pi_99__db - data::NakagamiRiceCurves[i][Y_pi_99_INDEX])) / (data::NakagamiRiceCurves[i][Y_pi_99_INDEX] - data::NakagamiRiceCurves[i - 1][Y_pi_99_INDEX]);
        K = (getK(i) * (Y_pi_99__db - getNakagamiRiceCurves(i - 1,Const.Y_pi_99_INDEX+1)) - getK(i - 1) * (Y_pi_99__db - getNakagamiRiceCurves(i, Const.Y_pi_99_INDEX+1) ) ) / ( getNakagamiRiceCurves(i, Const.Y_pi_99_INDEX+1) - getNakagamiRiceCurves(i - 1, Const.Y_pi_99_INDEX+1));
        return
    end
end

% no match.  Y_pi_99__db is greater than the data contains.  Return largest K
K = getK(end);
%return data::K.back();
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y_pi__db = NakagamiRice(K, p, Const)
%
%  Description:  This function computes the value of the Nakagami-Rice
%                distribution for K and p%
%
%        Input:  K         - K-value
%                p         - Time percentage
%
%      Outputs:  Y_pi__db  - Variability, in dB

getNakagamiRiceCurves = dataNakagamiRiceCurves();
getP = dataP();
getK = dataK();

%auto lower_K = lower_bound(data::K.begin(), data::K.end(), K);
%auto d_K = distance(data::K.begin(), lower_K);
d_K = distance_lower(getK, K);

%auto lower_p = lower_bound(data::P.begin(), data::P.end(), p);
%auto d_p = distance(data::P.begin(), lower_p);
d_p = distance_lower(getP, p);

%fprintf(1,'mat: K = %f, dk = %d\n', K, d_K);
%fprintf(1,'mat: p = %f, dp = %d\n', p, d_p);

if (d_K == 1) % K <= -40
    
    if (d_p == 1)
        Y_pi__db = getNakagamiRiceCurves(1,1);
        return
        %return data::NakagamiRiceCurves[0][0];
        
    else
        Y_pi__db = LinearInterpolation(getP(d_p), getNakagamiRiceCurves(1,d_p), getP(d_p - 1), getNakagamiRiceCurves(1,d_p - 1), p);
        return
        %return LinearInterpolation(data::P[d_p], data::NakagamiRiceCurves[0][d_p], data::P[d_p - 1], data::NakagamiRiceCurves[0][d_p - 1], p);
    end
    
elseif (d_K == length(getK) ) % K > 20
    
    if (d_p == 1)
        Y_pi__db = getNakagamiRiceCurves(d_K,1);
        
        return
        %return data::NakagamiRiceCurves[d_K - 1][0];
    else
        Y_pi__db = LinearInterpolation(getP(d_p), getNakagamiRiceCurves(d_K,d_p), getP(d_p - 1), getNakagamiRiceCurves(d_K, d_p - 1), p);
        return
        %return LinearInterpolation(data::P[d_p], data::NakagamiRiceCurves[d_K - 1][d_p], data::P[d_p - 1], data::NakagamiRiceCurves[d_K - 1][d_p - 1], p);
    end
else
    
    if (d_p == 1)
        Y_pi__db = LinearInterpolation(getK(d_K), getNakagamiRiceCurves(d_K, 1), getK(d_K - 1), getNakagamiRiceCurves(d_K - 1, 1), K);
        
        return
        %return LinearInterpolation(data::K[d_K], data::NakagamiRiceCurves[d_K][0], data::K[d_K - 1], data::NakagamiRiceCurves[d_K - 1][0], K);
    else
        
        % interpolate between K's at constant p first
        %double v1 = LinearInterpolation(data::K[d_K], data::NakagamiRiceCurves[d_K][d_p], data::K[d_K - 1], data::NakagamiRiceCurves[d_K - 1][d_p], K);
        v1 = LinearInterpolation(getK(d_K), getNakagamiRiceCurves(d_K,d_p), getK(d_K - 1), getNakagamiRiceCurves(d_K - 1, d_p), K);
        %double v2 = LinearInterpolation(data::K[d_K], data::NakagamiRiceCurves[d_K][d_p - 1], data::K[d_K - 1], data::NakagamiRiceCurves[d_K - 1][d_p] - 1, K);
        
        v2 = LinearInterpolation(getK(d_K), getNakagamiRiceCurves(d_K, d_p - 1), getK(d_K - 1), getNakagamiRiceCurves(d_K - 1, d_p-1), K);
        
        Y_pi__db = LinearInterpolation(getP(d_p), v1, getP(d_p - 1), v2, p);
        return
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn, M_d, A_d0, d_crx__km, CASE] = TranshorizonSearch(path, terminal_1, terminal_2, f__mhz, A_dML__db, M_d, A_d0, Const)

%
%    Description:  This file computes Step 6 in Annex 2, Section 3 of
%                  Recommendation ITU-R P.528-5, "Propagation curves for
%                  aeronautical mobile and radionavigation services using
%                  the VHF, UHF and SHF bands"
%
%          Input:  path              - Structure containing parameters dealing
%                                      with the propagation path
%                  terminal_1        - Structure containing parameters dealing
%                                      with the geometry of the low terminal
%                  terminal_2        - Structure containing parameters dealing
%                                      with the geometry of the high terminal
%                  f__mhz            - Frequency, in MHz
%                  A_dML__db         - Diffraction loss at d_ML, in dB
%
%        Outputs:  M_d               - Slope of the diffraction line
%                  A_d0              - Intercept of the diffraction line
%                  d_crx__km         - Final search distance, in km
%                  CASE              - Case as defined in Step 6.5
%                  rtn               - SUCCESS or warning code


CASE = Const.CONST_MODE__SEARCH;
k = 0;

tropo.A_s__db = 0;

% Step 6.1.  Initialize search parameters
d_search__km(1) = path.d_ML__km + 3;       % d', [Eqn 3-8]
d_search__km(2) = path.d_ML__km + 2;       % d", [Eqn 3-9]

A_s__db = [0 0];
M_s = 0;

SEARCH_LIMIT = 100; % 100 km beyond starting point

for i_search = 1: SEARCH_LIMIT
    
    A_s__db(2) = A_s__db(1);
    
    % Step 6.2
    %tropo = Troposcatter(path, terminal_1, terminal_2, d_search__km(1), f__mhz, Const);
    tropo = Troposcatter(terminal_1, terminal_2, d_search__km(1), f__mhz, Const);
    A_s__db(1) = tropo.A_s__db;
    
    % if loss is less than 20 dB, the result is not within valid part of model
    if (tropo.A_s__db < 20.0)
        
        d_search__km(2) = d_search__km(1);
        d_search__km(1) = d_search__km(1)+1;
        continue;
    end
    
    k = k+1;
    if (k <= 1) % need two points to draw a line and we don't have them both yet
        
        d_search__km(2) = d_search__km(1);
        d_search__km(1) = d_search__km(1)+1;
        continue;
    end
    
    % Step 6.3
    M_s = (A_s__db(1) - A_s__db(2)) / (d_search__km(1) - d_search__km(2));      % [Eqn 3-10]
    
    if (M_s <= M_d)
        
        d_crx__km = d_search__km(1);
        
        % Step 6.6
        A_d__db = M_d * d_search__km(2) + A_d0;                            % [Eqn 3-11]
        
        if (A_s__db(2) >= A_d__db)
            CASE = Const.CASE_1;
        else
            
            % Adjust the diffraction line to the troposcatter model
            M_d = (A_s__db(2) - A_dML__db) / (d_search__km(2) - path.d_ML__km);     % [Eqn 3-12]
            A_d0 = A_s__db(2) - (M_d * d_search__km(2));                            % [Eqn 3-13]
            
            CASE = Const.CASE_2;
        end
        
        rtn = Const.SUCCESS;
        return
    end
    
    d_search__km(2) = d_search__km(1);
    d_search__km(1) = d_search__km(1)+1;
end

% M_s was always greater than M_d.  Default to diffraction-only transhorizon model
CASE = Const.CONST_MODE__DIFFRACTION;
d_crx__km = d_search__km(2);

rtn = Const.WARNING__DFRAC_TROPO_REGION;
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function tropo = Troposcatter(path, terminal_1, terminal_2, d__km, f__mhz, Const)
function tropo = Troposcatter(terminal_1, terminal_2, d__km, f__mhz, Const)

%
%   Description:  This file computes the Troposcatter loss
%                 as described in Annex 2, Section 11 of
%                 Recommendation ITU-R P.528-5, "Propagation curves for
%                 aeronautical mobile and radionavigation services using
%                 the VHF, UHF and SHF bands"
%
%          Input:  path          - Struct containing path parameters
%                  terminal_1    - Struct containing low terminal parameters
%                  terminal_2    - Struct containing high terminal parameters
%                  d__km         - Path distance, in km
%                  f__mhz        - Frequency, in MHz
%
%        Outputs:  tropo         - Struct containing resulting parameters


tropo.d_s__km = d__km - terminal_1.d_r__km - terminal_2.d_r__km;       % [Eqn 11-2]

if (tropo.d_s__km <= 0.0)
    
    tropo.d_z__km = 0.0;
    tropo.A_s__db = 0.0;
    tropo.d_s__km = 0.0;
    tropo.h_v__km = 0.0;
    tropo.theta_s = 0.0;
    tropo.theta_A = 0.0;
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the geometric parameters
    %
    
    tropo.d_z__km = 0.5 * tropo.d_s__km;                               % [Eqn 11-6]
    
    A_m = 1 / Const.a_0__km;                                           % [Eqn 11-7]
    dN = A_m - (1.0 / Const.a_e__km);                                  % [Eqn 11-8]
    gamma_e__km = (Const.N_s * 1e-6) / dN;                       % [Eqn 11-9]
    
    z_a__km = 1.0 / (2 * Const.a_e__km) * ((tropo.d_z__km / 2)^2);     % [Eqn 11-10]
    z_b__km = 1.0 / (2 * Const.a_e__km) * ((tropo.d_z__km ^ 2));       % [Eqn 11-11]
    
    Q_o = A_m - dN;                                              % [Eqn 11-12]
    
    Q_a = A_m - dN / exp(min(35.0, z_a__km / gamma_e__km));      % [Eqn 11-13]
    Q_b = A_m - dN / exp(min(35.0, z_b__km / gamma_e__km));      % [Eqn 11-13]
    
    Z_a__km = (7.0 * Q_o + 6.0 * Q_a - Q_b) * ((tropo.d_z__km^2) / 96.0);  % [Eqn 11-14]
    Z_b__km = (Q_o + 2.0 * Q_a) * ((tropo.d_z__km^2) / 6.0);               % [Eqn 11-15]
    
    Q_A = A_m - dN / exp(min(35.0, Z_a__km / gamma_e__km));                     % [Eqn 11-16]
    Q_B = A_m - dN / exp(min(35.0, Z_b__km / gamma_e__km));                     % [Eqn 11-16]
    
    tropo.h_v__km = (Q_o + 2.0 * Q_A) * ((tropo.d_z__km ^ 2) / 6.0);        % [Eqn 11-17]
    
    tropo.theta_A = (Q_o + 4.0 * Q_A + Q_B) * tropo.d_z__km / 6.0;          % [Eqn 11-18]
    
    tropo.theta_s = 2 * tropo.theta_A;                                      % [Eqn 11-19]
    
    %
    % Compute the geometric parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the scattering efficiency term
    %
    epsilon_1 = 5.67e-6 * (Const.N_s ^ 2) - 0.00232 * Const.N_s + 0.031;           % [Eqn 11-20]
    epsilon_2 = 0.0002 * (Const.N_s ^ 2) - 0.06 * Const.N_s + 6.6;                % [Eqn 11-21]
    
    gamma = 0.1424 * (1.0 + epsilon_1 / exp(min(35.0, (tropo.h_v__km / 4.0)^6)));   % [Eqn 11-22]
    
    S_e__db = 83.1 - epsilon_2 / (1.0 + 0.07716 * (tropo.h_v__km ^ 2)) + 20 * log10( (0.1424 / gamma)^ 2 * exp(gamma * tropo.h_v__km));    % [Eqn 11-23]
    
    %
    % Compute the scattering efficiency term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the scattering volume term
    %
    
    X_A1__km2 = (terminal_1.h_e__km)^ 2 + 4.0 * (Const.a_e__km + terminal_1.h_e__km) * Const.a_e__km * ( sin(terminal_1.d_r__km / (Const.a_e__km * 2))^2);     % [Eqn 11-24]
    X_A2__km2 = (terminal_2.h_e__km)^ 2 + 4.0 * (Const.a_e__km + terminal_2.h_e__km) * Const.a_e__km * ( sin(terminal_2.d_r__km / (Const.a_e__km * 2))^2);     % [Eqn 11-24]
    
    ell_1__km = sqrt(X_A1__km2) + tropo.d_z__km;                        % [Eqn 11-25]
    ell_2__km = sqrt(X_A2__km2) + tropo.d_z__km;                        % [Eqn 11-25]
    ell__km = ell_1__km + ell_2__km;                                    % [Eqn 11-26]
    
    s = (ell_1__km - ell_2__km) / ell__km;                              % [Eqn 11-27]
    eta = gamma * tropo.theta_s * ell__km / 2;                          % [Eqn 11-28]
    
    kappa = f__mhz / 0.0477;                                            % [Eqn 11-29]
    
    rho_1__km = 2.0 * kappa * tropo.theta_s * terminal_1.h_e__km;       % [Eqn 11-30]
    rho_2__km = 2.0 * kappa * tropo.theta_s * terminal_2.h_e__km;       % [Eqn 11-30]
    
    SQRT2 = sqrt(2);
    
    A = (1 - s^2)^2;                                                    % [Eqn 11-36]
    
    X_v1 = (1 + s)^2 * eta;                                             % [Eqn 11-32]
    X_v2 = (1 - s)^2 * eta;                                             % [Eqn 11-33]
    
    q_1 = X_v1^2  + rho_1__km^2;                                        % [Eqn 11-34]
    q_2 = X_v2^2  + rho_2__km^2;                                        % [Eqn 11-35]
    
    % [Eqn 11-37]
    B_s = 6 + 8 * s^2 ...
        + 8 * (1.0 - s) * X_v1^2 * rho_1__km^2 / (q_1^2) ...
        + 8 * (1.0 + s) * X_v2^2 * rho_2__km^2 / (q_2^2) ...
        + 2 * (1.0 - s^2) * (1 + 2 * X_v1^2 / q_1) * (1 + 2 * X_v2^2 / q_2);
    
    % [Eqn 11-38]
    C_s = 12 ...
        * ((rho_1__km + SQRT2) / rho_1__km)^2  ...
        * ((rho_2__km + SQRT2) / rho_2__km)^2 ...
        * (rho_1__km + rho_2__km) / (rho_1__km + rho_2__km + 2 * SQRT2);
    
    temp = (A * eta^2 + B_s * eta) * q_1 * q_2 / ( (rho_1__km^2) * (rho_2__km^2));
    
    S_v__db = 10 * log10(temp + C_s);
    
    %
    % Compute the scattering volume term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tropo.A_s__db = S_e__db + S_v__db + 10.0 * log10(kappa * (tropo.theta_s^3) / ell__km);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_p = CombineDistributions(A_M, A_p, B_M, B_p, p)
%
%
%    Description:  This function combines two distributions A and B, returning
%                  the resulting percentile.
%
%          Input:  A_M   - Mean of distribution A
%                  A_p   - p% of distribution A
%                  B_M   - Mean of distribution B
%                  B_p   - p% of distribution B
%                  p     - Percentage
%
%         Returns: C_p   - p% of resulting distribution C



C_M = A_M + B_M;

Y_1 = A_p - A_M;
Y_2 = B_p - B_M;

Y_3 = sqrt((Y_1^2) + (Y_2^2));

if (p < 50)
    C_p = C_M + Y_3;
    return
else
    C_p = C_M - Y_3;
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    P.835  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T__kelvin = GlobalTemperature(h__km, Const)
%
%    Description:  The mean annual global reference atmospheric temperature,
%                  in Kelvin.
%
%          Input:  h__km         - Geometric height, in km
%
%        Returns:  T__kelvin     - Temperature, in Kelvin.
%                                  Or error code (negative number).
%

if (h__km < 0)
    T__kelvin = Const.ERROR_HEIGHT_TOO_SMALL;
    return
end

if (h__km > 100)
    T__kelvin = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end

if (h__km < 86)
    
    h_prime__km = ConvertToGeopotentialHeight(h__km);
    T__kelvin = GlobalTemperature_Regime1(h_prime__km, Const);
    return
    
else
    T__kelvin = GlobalTemperature_Regime2(h__km, Const);
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T__kelvin = GlobalTemperature_Regime1(h_prime__km, Const)
%
%    Description:  The mean annual global reference atmospheric temperature,
%                  in Kelvin, for the first height regime.
%                  See Equations (2a-g).
%
%          Input:  h_prime__km   - Geopotential height, in km'
%
%        Returns:  T__kelvin     - Temperature, in Kelvin.
%                                  Or error code (negative number).
if (h_prime__km < 0)
    T__kelvin = Const.ERROR_HEIGHT_TOO_SMALL;
    return
elseif (h_prime__km <= 11)
    T__kelvin = 288.15 - 6.5 * h_prime__km;
    return
elseif (h_prime__km <= 20)
    T__kelvin = 216.65;
    return
elseif (h_prime__km <= 32)
    T__kelvin = 216.65 + (h_prime__km - 20);
    return
elseif (h_prime__km <= 47)
    T__kelvin = 228.65 + 2.8 * (h_prime__km - 32);
    return
elseif (h_prime__km <= 51)
    T__kelvin = 270.65;
    return
elseif (h_prime__km <= 71)
    T__kelvin = 270.65 - 2.8 * (h_prime__km - 51);
    return
elseif (h_prime__km <= 84.852)
    T__kelvin = 214.65 - 2.0 * (h_prime__km - 71);
    return
else
    T__kelvin = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T__kelvin = GlobalTemperature_Regime2(h__km, Const)
% Description:     The mean annual global reference atmospheric temperature,
%                  in Kelvin, for the second height regime.
%                  See Equations (4a-b).
%
%          Input:  h__km         - Geometric height, in km
%
%        Returns:  T__kelvin     - Temperature, in Kelvin.
%                                  Or error code (negative number).
%
if (h__km < 86)
    T__kelvin = ERROR_HEIGHT_TOO_SMALL;
    return
elseif (h__km <= 91)
    T__kelvin = 186.8673;
    return
elseif (h__km <= 100)
    T__kelvin = 263.1905 - 76.3232 * sqrt(1 - ((h__km - 91) / 19.9429).^2);
    return
else
    T__kelvin = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p__hPa = GlobalPressure(h__km, Const)
%    Description:  The mean annual global reference atmospheric pressure,
%                  in hPa.
%
%          Input:  h__km         - Geometric height, in km
%
%        Returns:  p__hPa        - Dry air pressure, in hPa.
%                                  Or error code (negative number).

if (h__km < 0)
    p__hPa = Const.ERROR_HEIGHT_TOO_SMALL;
    return
end
if (h__km > 100)
    p__hPa = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end

if (h__km < 86)
    
    h_prime__km = ConvertToGeopotentialHeight(h__km);
    p__hPa = GlobalPressure_Regime1(h_prime__km, Const);
    return
    
else
    p__hPa = GlobalPressure_Regime2(h__km, Const);
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p__hPa = GlobalPressure_Regime1(h_prime__km, Const)
%   Description:  The mean annual global reference atmospheric pressure,
%                  in hPa, for the first height regime.  See Equations (3a-g).
%
%          Input:  h_prime__km   - Geopotential height, in km'
%
%        Returns:  p__hPa        - Dry air pressure, in hPa.
%                                  Or error code (negative number).
if (h_prime__km < 0)
    p__hPa = Const.ERROR_HEIGHT_TOO_SMALL;
    return
elseif (h_prime__km <= 11)
    p__hPa = 1013.25 * (288.15 / (288.15 - 6.5 * h_prime__km))^(-34.1632 / 6.5);
    return
elseif (h_prime__km <= 20)
    p__hPa = 226.3226 * exp(-34.1632 * (h_prime__km - 11) / 216.65);
    return
elseif (h_prime__km <= 32)
    p__hPa = 54.74980 * (216.65 / (216.65 + (h_prime__km - 20)))^34.1632;
    return
elseif (h_prime__km <= 47)
    p__hPa = 8.680422 * (228.65 / (228.65 + 2.8 * (h_prime__km - 32)))^(34.1632 / 2.8);
    return
elseif (h_prime__km <= 51)
    p__hPa = 1.109106 * exp(-34.1632 * (h_prime__km - 47) / 270.65);
    return
elseif (h_prime__km <= 71)
    p__hPa = 0.6694167 * (270.65 / (270.65 - 2.8 * (h_prime__km - 51)))^(-34.1632 / 2.8);
    return
elseif (h_prime__km <= 84.852)
    p__hPa = 0.03956649 * (214.65 / (214.65 - 2.0 * (h_prime__km - 71)))^(-34.1632 / 2.0);
    return
else
    p__hPa = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p__hPa = GlobalPressure_Regime2(h__km, Const)
%    Description:  The mean annual global reference atmospheric pressure,
%                  in hPa, for the second height regime.  See Equation (5).
%
%          Input:  h__km         - Geometric height, in km
%
%        Returns:  p__hPa        - Dry air pressure, in hPa.
%                                  Or error code (negative number).

if (h__km < 86)
    p__hPa = Const.ERROR_HEIGHT_TOO_SMALL;
    return
end
if (h__km > 100)
    p__hPa = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end

a_0 = 95.571899;
a_1 = -4.011801;
a_2 = 6.424731e-2;
a_3 = -4.789660e-4;
a_4 = 1.340543e-6;

p__hPa = exp(a_0 + a_1 * h__km + a_2 * (h__km ^ 2) + a_3 * (h__km ^ 3) + a_4 * (h__km ^ 4));
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rho = GlobalWaterVapourDensity(h__km, rho_0, Const)
%    Description:  The mean annual global reference atmospheric water vapour
%                  density, in g/m^3.  See Equation (6).
%
%          Input:  h__km         - Geometric height, in km
%                  rho_0         - Ground-level water vapour density, in g/m^3
%
%        Returns:  rho           - Water vapour density, in g/m^3.
%                                  Or error code (negative number).
if (h__km < 0)
    rho = Const.ERROR_HEIGHT_TOO_SMALL;
    return
end
if (h__km > 100)
    rho = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end

h_0__km = 2;     % scale height

rho =  rho_0 * exp(-h__km / h_0__km);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e__hPa = GlobalWaterVapourPressure(h__km, rho_0, Const)
%    Description:  The mean annual global reference atmospheric water vapour
%                  pressure, in hPa.
%
%          Input:  h__km         - Geometric height, in km
%                  rho_0         - Ground-level water vapour density, in g/m^3
%
%        Returns:  e__hPa        - Water vapour pressure, e(h), in hPa
%                                  Or error code (negative number).

if (h__km < 0)
    e__hPa = Const.ERROR_HEIGHT_TOO_SMALL;
    return
end
if (h__km > 100)
    e__hPa = Const.ERROR_HEIGHT_TOO_LARGE;
    return
end

drho = GlobalWaterVapourDensity(h__km, rho_0, Const);

if (h__km < 86)
    
    % convert to geopotential height
    h_prime__km = ConvertToGeopotentialHeight(h__km);
    T__kelvin = GlobalTemperature_Regime1(h_prime__km, Const);
    
else
    T__kelvin = GlobalTemperature_Regime2(h__km, Const);
end

e__hPa = WaterVapourDensityToPressure(rho, T__kelvin);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k_prime__km = ConvertToGeopotentialHeight(h__km)
%    Description:  Converts from geometric height, in km, to geopotential
%                  height, in km'.  See Equation (1a).
%
%          Input:  k__km         - Geometric height, in km
%
%        Returns:  k_prime__km   - Geopotential height, in km'

k_prime__km = (6356.766 * h__km) / (6356.766 + h__km);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k__km = ConvertToGeometricHeight(h_prime__km)
%   Description:  Converts from geopotential height, in km', to geometric
%                  height, in km.  See Equation (1b).
%
%          Input:  k_prime__km   - Geopotential height, in km'
%
%        Returns:  k__km         - Geometric height, in km

k__km =  (6356.766 * h_prime__km) / (6356.766 - h_prime__km);
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e__hPa = WaterVapourDensityToPressure(rho, T__kelvin)
% Description:  Converts water vapour density, in g/m^3, to water vapour
% pressure, in hPa.  See Equation (8).
%
% Input:  rho       - Water vapour density, rho(h), in g/m^3
% T__kelvin - Temperature, T(h), in Kelvin
%
% Returns:  e__hPa    - Water vapour pressure, e(h), in hPa



e__hPa = (rho * T__kelvin) / 216.7;
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    P.676  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = SlantPathAttenuation(f__ghz, h_1__km, h_2__km, beta_1__rad, Const)
% Calculation the slant path attenuation due to atmospheric gases


%     RayTraceConfig config;
%     config.temperature = GlobalTemperature;
%     config.dry_pressure = GlobalPressure;
%     config.wet_pressure = GlobalWetPressure;

if (beta_1__rad > pi / 2)
    
    % negative elevation angle
    % find h_G and then trace in each direction
    % see Section 2.2.2
    
    % compute refractive index at h_1
    p__hPa = GlobalPressure(h_1__km, Const);
    T__kelvin = GlobalTemperature(h_1__km, Const);
    e__hPa = GlobalWetPressure(h_1__km, Const);
    
    n_1 = RefractiveIndex(p__hPa, T__kelvin, e__hPa);
    
    % set initial h_G at mid-point between h_1 and surface of the earth
    % then binary search to converge
    h_G__km = h_1__km;
    delta = h_1__km / 2;
    diff = 100;
    
    
    while(1)
        
        if (diff > 0)
            h_G__km = h_G__km - delta;
        else
            h_G__km = h_G__km + delta;
        end
        delta = delta / 2;
        
        p__hPa = GlobalPressure(h_G__km, Const);
        T__kelvin = GlobalTemperature(h_G__km, Const);
        e__hPa = GlobalWetPressure(h_G__km, Const);
        
        n_G = RefractiveIndex(p__hPa, T__kelvin, e__hPa);
        
        grazing_term = n_G * (Const.a_0__km + h_G__km);
        start_term = n_1 * (Const.a_0__km + h_1__km) * sin(beta_1__rad);
        
        diff = grazing_term - start_term;
        if (abs(diff) <= 0.001)
            break
        end
    end
    
    % converged on h_G.  Now call RayTrace in both directions with grazing angle
    
    beta_graze__rad = pi / 2;
    result_1 = RayTrace(f__ghz, h_G__km, h_1__km, beta_graze__rad, Const);
    result_2 = RayTrace(f__ghz, h_G__km, h_2__km, beta_graze__rad, Const);
    
    result.angle__rad = result_2.angle__rad;
    result.A_gas__db = result_1.A_gas__db + result_2.A_gas__db;
    result.a__km = result_1.a__km + result_2.a__km;
    result.bending__rad = result_1.bending__rad + result_2.bending__rad;
    result.delta_L__km = result_1.delta_L__km + result_2.delta_L__km;
    
else
    
    result = RayTrace(f__ghz, h_1__km, h_2__km, beta_1__rad, Const);
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = RefractiveIndex(p__hPa, T__kelvin, e__hPa)
%    Description:  Compute the refractive index.
%
%          Input:  p__hPa        - Dry pressure, in hPa
%                  T__kelvin     - Temperature, in Kelvin
%                  e__hPa        - Water vapour pressure, in hPa
%
%        Returns:  n             - Refractive index

% dry term of refractivity
N_dry = 77.6 * p__hPa / T__kelvin;

% wet term of refractivity
N_wet = 72 * e__hPa / T__kelvin + 3.75e5 * e__hPa / (T__kelvin ^ 2);

N = N_dry + N_wet;

n = 1 + N * (10^(-6));

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = RayTrace(f__ghz, h_1__km, h_2__km, beta_1__rad, Const)
%
%    Description:  Traces the ray from terminal h_1 to terminal h_2 and
%                  computes results such as atmospheric absorption loss and
%                  ray path length.
%
%          Input:  f__ghz        - Frequency, in GHz
%                  h_1__km       - Height of the low terminal, in km
%                  h_2__km       - Height of the high terminal, in km
%                  beta_1__rad   - Elevation angle (from zenith), in rad
%
%         Output:  result        - Ray trace result structure
%
%         Returns:  [void]

% Equations 16(a)-(c)
i_lower = floor(100 * log(1e4 * h_1__km * (exp(1. / 100.) - 1) + 1) + 1);
i_upper = ceil(100 * log(1e4 * h_2__km * (exp(1. / 100.) - 1) + 1) + 1);
m = ((exp(2. / 100.) - exp(1. / 100.)) / (exp(i_upper / 100.) - exp(i_lower / 100.))) * (h_2__km - h_1__km);

alpha_i__rad = beta_1__rad;
beta_ii__rad = beta_1__rad;

% initialize results
result.A_gas__db = 0;
result.bending__rad = 0;
result.a__km = 0;
result.delta_L__km = 0;

% initialize starting layer
delta_i__km = LayerThickness(m, i_lower);
h_i__km = h_1__km + m * ((exp((i_lower - 1) / 100.) - exp((i_lower - 1) / 100.)) / (exp(1 / 100.) - 1));
[n_i, gamma_i] = GetLayerProperties(f__ghz, h_i__km + delta_i__km / 2, Const);
r_i__km = Const.a_0__km + h_i__km;

% record bottom layer properties for alpha and beta calculations
r_1__km = r_i__km;
n_1 = n_i;

% summation from Equation 13
for i = i_lower:i_upper - 1
    
    delta_ii__km = LayerThickness(m, i + 1);
    h_ii__km = h_1__km + m * ((exp((i + 1 - 1) / 100.) - exp((i_lower - 1) / 100.)) / (exp(1 / 100.) - 1));
    
    [n_ii, gamma_ii] = GetLayerProperties(f__ghz, h_ii__km + delta_ii__km / 2, Const);
    
    r_ii__km = Const.a_0__km + h_ii__km;
    
    delta_i__km = LayerThickness(m, i);
    
    % Equation 19b
    beta_i__rad = asin(min(1, (n_1 * r_1__km) / (n_i * r_i__km) * sin(beta_1__rad)));
    
    % entry angle into the layer interface, Equation 18a
    alpha_i__rad = asin(min(1, (n_1 * r_1__km) / (n_i * r_ii__km) * sin(beta_1__rad)));
    
    % path length through ith layer, Equation 17
    a_i__km = -r_i__km * cos(beta_i__rad) + sqrt( (r_i__km^2) * (cos(beta_i__rad))^2 + 2 * r_i__km * delta_i__km + (delta_i__km^2));
    
    result.a__km = result.a__km+ a_i__km;
    result.A_gas__db = result.A_gas__db + a_i__km * gamma_i;
    result.delta_L__km = result.delta_L__km + a_i__km * (n_i - 1);     % summation, Equation 23
    
    beta_ii__rad = asin(n_i / n_ii * sin(alpha_i__rad));
    
    % summation of the bending angle, Equation 22a
    % the summation only goes to i_max - 1
    if (i ~= i_upper - 1)
        result.bending__rad = result.bending__rad + beta_ii__rad - alpha_i__rad;
    end
    % shift for next loop
    h_i__km = h_ii__km;
    n_i = n_ii;
    gamma_i = gamma_ii;
    r_i__km = r_ii__km;
end

result.angle__rad = alpha_i__rad;
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta_i__km = LayerThickness(m, i)
%
%    Description:  Thickness of the ith layer.
%
%          Input:  m             - Internal parameter
%                  i             - Layer of interest
%
%        Returns:  delta_i__km   - Layer thickness, in km
%

% Equation 14
delta_i__km = m * exp((i - 1) / 100.);

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n, gamma]= GetLayerProperties(f__ghz, h_i__km, Const)
%
%    Description:  Determine the parameters for the ith layer
%
%          Input:  f__ghz        - Frequency, in GHz
%                  h_i__km       - Height of the ith layer, in km
%                  config        - Structure containing atmospheric params
%
%         Output:  n             - Refractive index
%                  gamma         - Specific attenuation, in dB/km
%

% use function pointers to get atmospheric parameters
T__kelvin = GlobalTemperature(h_i__km, Const);
p__hPa = GlobalPressure(h_i__km, Const);
e__hPa = GlobalWetPressure(h_i__km, Const);

% compute the refractive index for the current layer
n = RefractiveIndex(p__hPa, T__kelvin, e__hPa);

% specific attenuation of layer
gamma = SpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa);

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e__hPa = GlobalWetPressure(h__km, Const)

T__kelvin = GlobalTemperature(h__km, Const);
P__hPa = GlobalPressure(h__km, Const);
rho__g_m3 = max(GlobalWaterVapourDensity(h__km, Const.RHO_0__M_KG, Const), 2 * (10^(-6)) * 216.7 * P__hPa / T__kelvin);

e__hPa = WaterVapourDensityToPressure(rho__g_m3, T__kelvin);

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_i = LineShapeFactor(f__ghz, f_i__ghz, delta_f__ghz, delta)

term1 = f__ghz ./ f_i__ghz;
term2 = (delta_f__ghz - delta .* (f_i__ghz - f__ghz)) ./ ( (f_i__ghz - f__ghz).^2  + (delta_f__ghz.^2));
term3 = (delta_f__ghz - delta .* (f_i__ghz + f__ghz)) ./ ( (f_i__ghz + f__ghz).^2  + (delta_f__ghz.^2));

F_i = term1 .* (term2 + term3);

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N_D = NonresonantDebyeAttenuation(f__ghz, e__hPa, p__hPa, theta)

% width parameter for the Debye spectrum, Equation 9
d = 5.6e-4 * (p__hPa + e__hPa) * theta^(0.8);

% Equation 8
frac_1 = 6.14e-5 / (d * (1 + (f__ghz / d)^2));
frac_2 = (1.4e-12 * p__hPa * theta^(1.5)) / (1 + 1.9e-5 * (f__ghz^(1.5)));
N_D = f__ghz * p__hPa * (theta^2) * (frac_1 + frac_2);

return
end

function result = getOxygenData()
% Spectroscopic data for oxygen attenuation (Table 1)

%
result.f_0 = [50.474214,  50.987745,  51.503360,  52.021429,  52.542418,  53.066934,  53.595775, ...
    54.130025,  54.671180,  55.221384,  55.783815,  56.264774,  56.363399,  56.968211, ...
    57.612486,  58.323877,  58.446588,  59.164204,  59.590983,  60.306056,  60.434778, ...
    61.150562,  61.800158,  62.411220,  62.486253,  62.997984,  63.568526,  64.127775, ...
    64.678910,  65.224078,  65.764779,  66.302096,  66.836834,  67.369601,  67.900868, ...
    68.431006,  68.960312, 118.750334, 368.498246, 424.763020, 487.249273, ...
    715.392902, 773.839490, 834.145546 ...
    ];

result.a_1 = [ ...
    0.975,    2.529,    6.193,   14.320,   31.240,   64.290,  124.600,  227.300, ...
    389.700,  627.100,  945.300,  543.400, 1331.800, 1746.600, 2120.100, 2363.700, ...
    1442.100, 2379.900, 2090.700, 2103.400, 2438.000, 2479.500, 2275.900, 1915.400, ...
    1503.000, 1490.200, 1078.000,  728.700,  461.300,  274.000,  153.000,   80.400, ...
    39.800,   18.560,    8.172,    3.397,    1.334,  940.300,   67.400,  637.700, ...
    237.400,   98.100,  572.300,  183.100 ...
    ];

result.a_2 = ...
    [ ...
    9.651, 8.653, 7.709, 6.819, 5.983, 5.201, 4.474, 3.800, 3.182, 2.618, 2.109, ...
    0.014, 1.654, 1.255, 0.910, 0.621, 0.083, 0.387, 0.207, 0.207, 0.386, 0.621, ...
    0.910, 1.255, 0.083, 1.654, 2.108, 2.617, 3.181, 3.800, 4.473, 5.200, 5.982, ...
    6.818, 7.708, 8.652, 9.650, 0.010, 0.048, 0.044, 0.049, 0.145, 0.141, 0.145 ...
    ];

result.a_3 = ...
    [ ...
    6.690,  7.170,  7.640,  8.110,  8.580,  9.060,  9.550,  9.960, 10.370, ...
    10.890, 11.340, 17.030, 11.890, 12.230, 12.620, 12.950, 14.910, 13.530, ...
    14.080, 14.150, 13.390, 12.920, 12.630, 12.170, 15.130, 11.740, 11.340, ...
    10.880, 10.380,  9.960,  9.550,  9.060,  8.580,  8.110,  7.640,  7.170, ...
    6.690, 16.640, 16.400, 16.400, 16.000, 16.000, 16.200, 14.700 ...
    ];

result.a_4 = ...
    [ ...
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
    0.0, 0.0 ...
    ];

result.a_5 = ...
    [ ...
    2.566,  2.246,  1.947,  1.667,  1.388,  1.349,  2.227,  3.170,  3.558,  2.560, ...
    -1.172,  3.525, -2.378, -3.545, -5.416, -1.932,  6.768, -6.561,  6.957, -6.395, ...
    6.342,  1.014,  5.014,  3.029, -4.499,  1.856,  0.658, -3.036, -3.968, -3.528, ...
    -2.548, -1.660, -1.680, -1.956, -2.216, -2.492, -2.773, -0.439,  0.000,  0.000, ...
    0.000,  0.000,  0.000,  0.000 ...
    ];

result.a_6 = ...
    [ ...
    6.850,  6.800,  6.729,  6.640,  6.526,  6.206,  5.085,  3.750,  2.654,  2.952, ...
    6.135, -0.978,  6.547,  6.451,  6.056,  0.436, -1.273,  2.309, -0.776,  0.699, ...
    -2.825, -0.584, -6.619, -6.759,  0.844, -6.675, -6.139, -2.895, -2.590, -3.680, ...
    -5.002, -6.091, -6.393, -6.475, -6.545, -6.600, -6.650,  0.079,  0.000,  0.000, ...
    0.000,  0.000,  0.000,  0.000 ...
    ];

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = getWaterVapourData()
% Spectroscopic data for water vapor attenuation (Table 2)

%
result.f_0 = ...
    [ ...
    22.235080,  67.803960, 119.995940, 183.310087, 321.225630, 325.152888,  336.227764, ...
    380.197353, 390.134508, 437.346667, 439.150807, 443.018343, 448.001085,  470.888999, ...
    474.689092, 488.490108, 503.568532, 504.482692, 547.676440, 552.020960,  556.935985, ...
    620.700807, 645.766085, 658.005280, 752.033113, 841.051732, 859.965698,  899.303175, ...
    902.611085, 906.205957, 916.171582, 923.112692, 970.315022, 987.926764, 1780.000000 ...
    ];

result.b_1 = ...
    [ ...
    0.1079, 0.0011,   0.0007,  2.273, 0.0470, 1.514,    0.0010, 11.67,   0.0045, ...
    0.0632, 0.9098,   0.1920, 10.41,  0.3254, 1.260,    0.2529,  0.0372, 0.0124, ...
    0.9785, 0.1840, 497.0,     5.015, 0.0067, 0.2732, 243.4,     0.0134, 0.1325, ...
    0.0547, 0.0386,   0.1836,  8.400, 0.0079, 9.009,  134.6,     17506.0 ...
    ];

result.b_2 = ...
    [ ...
    2.144, 8.732, 8.353, .668, 6.179, 1.541, 9.825, 1.048, 7.347, 5.048, ...
    3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 6.731, 6.731, .158, .158, ...
    .159, 2.391, 8.633, 7.816, .396, 8.177, 8.055, 7.914, 8.429, 5.110, ...
    1.441, 10.293, 1.919, .257, .952 ...
    ];

result.b_3 = ...
    [ ...
    26.38, 28.58, 29.48, 29.06, 24.04, 28.23, 26.93, 28.11, 21.52, 18.45, 20.07, ...
    15.55, 25.64, 21.34, 23.20, 25.86, 16.12, 16.12, 26.00, 26.00, 30.86, 24.38, ...
    18.00, 32.10, 30.86, 15.90, 30.60, 29.85, 28.65, 24.08, 26.73, 29.00, 25.50, ...
    29.85, 196.3 ...
    ];

result.b_4 = ...
    [ ...
    .76, .69, .70, .77, .67, .64, .69, .54, .63, .60, .63, .60, .66, .66, ...
    .65, .69, .61, .61, .70, .70, .69, .71, .60, .69, .68, .33, .68, .68, ...
    .70, .70, .70, .70, .64, .68, 2.00 ...
    ];

result.b_5 = ...
    [ ...
    5.087, 4.930, 4.780, 5.022, 4.398, 4.893, 4.740, 5.063, 4.810, 4.230, 4.483, ...
    5.083, 5.028, 4.506, 4.804, 5.201, 3.980, 4.010, 4.500, 4.500, 4.552, 4.856, ...
    4.000, 4.140, 4.352, 5.760, 4.090, 4.530, 5.100, 4.700, 5.150, 5.000, 4.940, ...
    4.550, 24.15 ...
    ];

result.b_6 = ...
    [ ...
    1.00, .82, .79, .85, .54, .74, .61, .89, .55, .48, .52, .50, .67, .65, ...
    .64, .72, .43, .45, 1.00, 1.00, 1.00, .68, .50, 1.00, .84, .45, .84, ...
    .90, .95, .53, .78, .80, .67, .90, 5.00 ...
    ];
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N_o = OxygenRefractivity(f__ghz, T__kelvin, e__hPa, p__hPa, OxygenData)

theta = 300 / T__kelvin;

N = 0;

%     for i = 1 : length(OxygenData.f_0)
%
%         % Equation 3, for oxygen
%         S_i = OxygenData.a_1(i) * 1e-7 * p__hPa * (theta^3) * exp(OxygenData.a_2(i) * (1 - theta));
%
%         % compute the width of the line, Equation 6a, for oxygen
%         delta_f__ghz = OxygenData.a_3(i) * 1e-4 * (p__hPa * (theta^(0.8 - OxygenData.a_4(i))) + 1.1  * e__hPa * theta);
%
%         % modify the line width to account for Zeeman splitting of the oxygen lines
%         % Equation 6b, for oxygen
%         delta_f__ghz = sqrt((delta_f__ghz^2) + 2.25e-6);
%
%         % correction factor due to interference effects in oxygen lines
%         % Equation 7, for oxygen
%         delta = (OxygenData.a_5(i) + OxygenData.a_6(i) * theta) * 1e-4 * (p__hPa + e__hPa) * (theta^0.8);
%
%         F_i = LineShapeFactor(f__ghz, OxygenData.f_0(i), delta_f__ghz, delta);
%
%         % summation of terms...from Equation 2a, for oxygen
%         N = N +  S_i * F_i;
%     end

S_i = OxygenData.a_1 .* 1e-7 .* p__hPa .* (theta^3) .* exp(OxygenData.a_2 .* (1 - theta));
delta_f__ghz = OxygenData.a_3 .* 1e-4 .* (p__hPa .* (theta.^(0.8 - OxygenData.a_4)) + 1.1  .* e__hPa .* theta);
delta_f__ghz = sqrt((delta_f__ghz.^2) + 2.25e-6);
delta = (OxygenData.a_5 + OxygenData.a_6 .* theta) .* 1e-4 .* (p__hPa + e__hPa) .* (theta^0.8);
F_i = LineShapeFactor(f__ghz, OxygenData.f_0, delta_f__ghz, delta);

N_D = NonresonantDebyeAttenuation(f__ghz, e__hPa, p__hPa, theta);
N = sum(S_i .* F_i);

N_o = N + N_D;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N_w = WaterVapourRefractivity(f__ghz, T__kelvin, e__hPa, P__hPa, WaterVapourData)

theta = 300 / T__kelvin;

N_w = 0;



%     for i = 1 : length(WaterVapourData.f_0)
%
%         % Equation 3, for water vapour
%         S_i = 0.1 * WaterVapourData.b_1(i) * e__hPa * (theta^3.5) * exp(WaterVapourData.b_2(i) * (1 - theta));
%
%         % compute the width of the line, Equation 6a, for water vapour
%         delta_f__ghz = 1e-4 * WaterVapourData.b_3(i) * (P__hPa * (theta^(WaterVapourData.b_4(i))) + WaterVapourData.b_5(i) * e__hPa * (theta^(WaterVapourData.b_6(i))));
%
%         % modify the line width to account for Doppler broadening of water vapour lines
%         % Equation 6b, for water vapour
%         term1 = 0.217 * (delta_f__ghz^2) + (2.1316e-12 * (WaterVapourData.f_0(i)^2) / theta);
%         delta_f__ghz = 0.535 * delta_f__ghz + sqrt(term1);
%
%         % Equation 7, for water vapour
%         delta = 0;
%
%         F_i = LineShapeFactor(f__ghz, WaterVapourData.f_0(i), delta_f__ghz, delta);
%
%         % summation of terms...from Equation 2b, for water vapour
%         N_w = N_w  + S_i * F_i;
%     end

S_i = 0.1 * WaterVapourData.b_1 .* e__hPa .* (theta.^3.5) .* exp(WaterVapourData.b_2 .* (1 - theta));
delta_f__ghz = 1e-4 * WaterVapourData.b_3 .* (P__hPa .* (theta.^(WaterVapourData.b_4)) + WaterVapourData.b_5 .* e__hPa .* (theta.^(WaterVapourData.b_6)));
term1 = 0.217 .* (delta_f__ghz.^2) + (2.1316e-12 .* (WaterVapourData.f_0.^2) ./ theta);
delta_f__ghz = 0.535 .* delta_f__ghz + sqrt(term1);
delta = 0;
F_i = LineShapeFactor(f__ghz, WaterVapourData.f_0, delta_f__ghz, delta);
N_w = sum(S_i .* F_i);


return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = SpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa)

OxygenData = getOxygenData();
WaterVapourData = getWaterVapourData();

gamma_o = OxygenSpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa, OxygenData);
gamma_w = WaterVapourSpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa, WaterVapourData);

gamma = gamma_o + gamma_w;   % [Eqn 1]

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma_o = OxygenSpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa, OxygenData)

% partial Eqn 1
N_o = OxygenRefractivity(f__ghz, T__kelvin, e__hPa, p__hPa, OxygenData);
gamma_o = 0.1820 * f__ghz * N_o;

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma_w = WaterVapourSpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa, WaterVapourData)

% partial Eqn 1
N_w = WaterVapourRefractivity(f__ghz, T__kelvin, e__hPa, p__hPa, WaterVapourData);
gamma_w = 0.1820 * f__ghz * N_w;

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A__db = TerrestrialPath(f__ghz, T__kelvin, e__hPa, p__hPa, r_0__km)

gamma = SpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa);

% Equation 10
A__db = gamma * r_0__km;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e__hPa = WaterVapourDensityToPartialPressure(rho__g_m3, T__kelvin)

% Equation 4
e__hPa = (rho__g_m3 * T__kelvin) / 216.7;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = LineOfSightParams()
% Initilizes the structure LineOfSightParams
% Heights
result.z__km = [0 0 ];

% Distances
result.d__km = 0;           % Path distance between terminals
result.r_0__km = 0;         % Direct ray length
result.r_12__km = 0;        % Indirect ray length
result.D__km = [0 0 ];

% Angles
result.theta_h1__rad = 0;      %Take-off angle from low terminal to high terminal, in rad
result.theta_h2__rad = 0;      % Take-off angle from high terminal to low terminal, in rad
result.theta = [0 0];

% Misc
result.a_a__km = 0;             % Adjusted earth radius
result.delta_r__km = 0;         % Ray length path difference
result.A_LOS__db = 0;           % Loss due to LOS path
return
end