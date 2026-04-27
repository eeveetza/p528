clear all
close all
fclose all;
clc
use_reflection = false;

% immediate printing to command window in octave
% if (isOctave)
%     page_screen_output(0);
%     page_output_immediately(1);
% end

path{1} = './Data Tables';

%% begin code
% Collect all the filenames .csv in the folder pathname that contain the profile data
ff = [1000, 5000, 10000];
ffstr ={'1,000'; '5,000'; '10,000'};
pp = [0.01, 0.05, 0.1, 0.5, 0.95];
h1 = 1.5;
h2 = 300;



cnt_fail = 0;
cnt_pass = 0;

for kk=1:length(ff)

    for ii=1:length(pp)
        %for ii = 1:1
        kindex=1;
        filename1 = [ffstr{kk} ' MHz - Lb(' num2str(pp(ii)) ')_P528_wo_reflection.csv'];
        fprintf(1,'***********************************************\n');
        fprintf(1,' Processing file %s ...\n', filename1);
        fprintf(1,'***********************************************\n');



        fid=fopen(filename1,'w');
        if (fid==-1)
            return;
        end


        % First line is of the following format
        % 1200MHz / Lb(0.01) dB
        firstline = [num2str(ff(kk)) 'MHz / Lb(' num2str(pp(ii)) ') dB\n'];

        fprintf(fid, firstline);
        fprintf(fid, ',h2(m),300\n');
        fprintf(fid, ',h1(m),1.5\n');
        fprintf(fid, 'D(km),FSL\n');

        d = 10:1:100;

        for dd = 1:length(d)

            result = tl_p528(d(dd),h1, h2, ff(kk), 0, pp(ii)*100, use_reflection);


            fprintf(fid,'%f,%.1f,%.1f\n', d(dd), result.A_fs__db, result.A__db );

        end

    end

fclose(fid);
end

fprintf(1,'Successfully passed %d out of %d tests\n', cnt_pass, cnt_pass+cnt_fail);