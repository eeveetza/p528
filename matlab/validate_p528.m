clear all
close all
fclose all;
clc


% immediate printing to command window in octave
% if (isOctave)
%     page_screen_output(0);
%     page_output_immediately(1);
% end

path{1} = './Data Tables';

%% begin code
% Collect all the filenames .csv in the folder pathname that contain the profile data
for i = 1:length(path)
    filenames{i} = dir(fullfile(path{i}, '*.csv')); % filenames{1}(i).name is the filename
end

cnt_fail = 0;
cnt_pass = 0;

for pp=1:length(path)
    datanumber=length(filenames{pp});
    
    for ii=1:datanumber
        %for ii = 1:1
        kindex=1;
        filename1 = filenames{pp}(ii).name;
        fprintf(1,'***********************************************\n');
        fprintf(1,' Processing file %d/%d: %s%s ...\n', ii, datanumber, path{pp}, filename1);
        fprintf(1,'***********************************************\n');
        filename=[path{pp} '/' filenames{pp}(ii).name];
        fid=fopen(filename,'r');
        if (fid==-1)
            return;
        end
        
        
        % First line is of the following format
        % 1200MHz / Lb(0.01) dB
        readLine = fgetl(fid);
        rgx = '[-+]?\d+\.?\d*([eE][-+]?\d+)?';
        dummy=regexp(readLine,rgx,'match');
        f = str2double(dummy{1});
        p = str2double(dummy{2});
        % read the h2 values
        
        readLine = fgetl(fid);
        dummy = regexp(readLine,',','split');
        
        for i = 1:length(dummy(3:end))-3
            h2(i) = str2double(dummy(i+2));
        end
        
        % read the h2 values
        readLine = fgetl(fid);
        dummy = regexp(readLine,',','split');
        
        for i = 1:length(dummy(3:end))-3
            h1(i) = str2double(dummy(i+2));
        end
        
        fgetl(fid);
        
        
        
        fprintf(1,'%20s %20s %20s \n', 'MATLAB', 'REF TABLE', 'DELTA');
        count = 1;
        while(1)
            readLine=fgetl(fid);
            if (readLine==-1)
                break
            end
            
            dummy=regexp(readLine,',','split');
            D(count) = str2double(dummy(1));
            FSL(count) = str2double(dummy(2));
            for i = 1:length(dummy(3:end))-3
                tl_ref(count,i) = str2double(dummy(i+2));
            end
            
            count = count + 1;
        end
        fclose(fid);
        delta = 0;
        
        
        
        for i = 1:1000:length(D)
            for j = 1:length(h1)
                
                result = tl_p528(D(i),h1(j), h2(j), f, 0, p*100);
                
                delta = round(10.0 * (result.A__db-tl_ref(i,j)) ) / 10.0;
                
                if (abs(delta)>0.1)
                    cnt_fail = cnt_fail + 1;
                else
                    cnt_pass = cnt_pass + 1;
                end
               
                fprintf(1,'%20.1f %20.1f %20.1f \n', result.A__db, tl_ref(i,j) , delta );

            end
            
        end
        
        
        count = count + 1;
    end
    
    
end

fprintf(1,'Successfully passed %d out of %d tests\n', cnt_pass, cnt_pass+cnt_fail);