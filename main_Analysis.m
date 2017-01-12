
clear all;
clc;
close all;

%% Select Movies

% load movie file
[files,filepath,numFiles] = SelectFilesToAnalyze;

%%

fullAnalysis = 1; %

for z = 1:numFiles
    close all;
    
    % load movie
    disp(files{z});
    if fullAnalysis == 1
        [mov numFrames dirpath frameDuration fname] = LoadMovie(filepath{z},files{z});
    else
        [tmp fname] = fileparts(files{z});
        dirpath = fullfile(filepath{z},fname);
        frameDuration = 0.032; 
    end
    
    % load gaussian fit positions
    try
        pixelPos = dlmread(fullfile(dirpath,'pixelPositions.txt'));
        status = 1;
    catch
        disp('Pixel position does not exist');
        status = 0;
    end
    
    if status == 1
        pixelPos(isnan(pixelPos(:,end)),:) = [];
        if size(pixelPos,2) == 7
            pixelPos = pixelPos(:,[1:5 7]); %excluded column 6 (snr)
        end
        % [ x y intensity error frameIndex SNR]

        
        % Filter positions!!!!!!!!!!
        

        % Track particles
        minFrames = 10;
        max_dsp = 10;
        trackParam.mem = 15;      % number of frames that the fluorophore can disappear %0
        trackParam.good = minFrames;     % minimum number of frames that is considered a good track
        trackParam.dim = 2;
        trackParam.quiet = 1;
        trackParam.pixels2um = .086*2;
        
        set_parameters = [minFrames max_dsp trackParam.mem];
        try
            %[pixelTracks, posTracks] = FindTracks(filteredPixelPos,max_dsp,trackParam);
            [pixelTracks, posTracks] = FindTracks(pixelPos,max_dsp,trackParam); %
            trackStatus = 1;
        catch
            trackStatus = 0;
        end

        if trackStatus == 1
            % save tracks
            dlmwrite(fullfile(dirpath,'PixelTracks.txt'),pixelTracks,'delimiter','\t','newline','pc');
            dlmwrite(fullfile(dirpath,'PosTracks.txt'),posTracks,'delimiter','\t','newline','pc');
            dlmwrite(fullfile(dirpath,'analysis_parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
        else
            disp(['Errors tracking in ' fname]);
        end
    end
end






