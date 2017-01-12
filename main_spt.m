%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script prompts the user to select multiple .tif files movie for single
% molecule particle tracking analysis.  All of the results are saved within
% their respective folder where the selected movie is located. This one is 
% fully automated.  The user will only have to input parameters at the
% beginning.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

%% Select Movies

%addpath('./bioformats');

% load tif movie file
[files,filepath,numFiles] = SelectFilesToAnalyze;

%% Extract Positions

particleSiz = 7; %3 original, 5 good
brightThreshold = 7000; %200 original, 5000 good

showplot = 1;
showFit = 0; %0
showcnt = 0; %
set_parameters = [particleSiz brightThreshold];
for z = 1:numFiles
    tic
    close all
    disp(files{z})
    [mov,numFrames,dirpath,frameRate,fname] = LoadMovie(filepath{z},files{z});
    fits = []; %array containing number of fitted centroids for each frame
    % get flurophore positions
    pixelPos = [];
    for frameIndex = 1:numFrames 
        %for frameIndex = 1:1
        img = mov(:,:,frameIndex);
        
        % bandpass filter image
        imgb = bpass(double(img),1,particleSiz);

        % find peaks
        pk = pkfnd(double(imgb),brightThreshold,particleSiz);
        
        % if peaks are found, find centroids
        if ~isempty(pk)
            
            % calculate centroid moments
            cnt = cntrd(double(imgb),pk,particleSiz+1);
        
            
             % show movie with initial centroids
%             if showcnt == 1
%                 %img = mov(:,:,frameIndex); %
%                 figure(1); cla; hold on;
%                 imagesc(img); colormap gray;
%                 scatter(cnt(:,1),cnt(:,2),200,'gs');
%                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                 axis tight
%                 title(['centroids ' num2str(frameIndex) ' out of ' num2str(numFrames)]);
%             end
            
            % 2D Gaussian fit around each centroid
            scale = 5; %5
            p0(3) = ceil(particleSiz/2);
            p0(4) = max(img(:));
            p0(5) = min(img(:));
            centroid = []; badIndex = []; error = []; sigma = []; 
            intensity = []; SNR = []; ratio = []; 
            
            
            
            for i = 1:size(cnt,1) %loop through every centroid
            %for i = 1:5 %loop through every centroid
                 try 
                    p0(1:2) = cnt(i,1:2); %get x y positions of centroid 
                    %[parameters,img,e1] = Gaussian2DFiterror(img,cnt(i,1:2),p0,scale,0);
                    %[parameters,img] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                    [parameters,imgbb] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                    centroid = [centroid; parameters(1:2)]; % position of fit
                    meanSigma = parameters(3); % standard deviation of fit
                    N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
                    sigma = [sigma; meanSigma];
                    intensity = [intensity; N];
                    error = [error; meanSigma/sqrt(N)];
                    SNR = [SNR; parameters(4)/parameters(5)];
                    
                    fits1 = [frameIndex length(centroid)];

                    % show 2D Gaussian Fit
                    if showFit == 1
                        disp(['SNR = ' num2str(SNR(i)) '; Error = ' num2str(error(i))]);
                        [sizey sizex] = size(img);
                        [X,Y] = meshgrid(1:sizey,1:sizex);
                        figure(1); cla; 
                        subplot(2,2,1); imagesc(img); colormap gray;
                        subplot(2,2,3); surf(img); 
                        subplot(2,2,2);
                        contour3(img,'--');
                        hold on;
                        contour3(Gaussian2D(parameters,X,Y),'-');
                    end
                    
               

                catch
                    badIndex = [badIndex; i];
                    disp('Error fitting Gaussian');
                end

            end
            
            fits = vertcat(fits,fits1); %add frame number v. # of centroids fitted to fits array
            % store fit info
            pixelPos = [pixelPos; centroid intensity sigma error SNR ones(size(centroid,1),1)*frameIndex];
        
            % show movie with fitted centroids
            if showplot == 1
                %img = mov(:,:,frameIndex); %
                figure(2); cla; hold on;
                %figure; cla; hold on;
                imagesc(img); colormap gray;
                scatter(centroid(:,1),centroid(:,2),200,'gs');
                %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
                axis tight
                title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
            end
            
            % show specific frames
%             if frameIndex == 1 || 1000
%                 %img = mov(:,:,frameIndex); %
%                 figure; cla; hold on;
%                 imagesc(img); colormap gray;
%                 scatter(centroid(:,1),centroid(:,2),200,'gs');
%                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                 axis tight
%                 title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
%             end
        end        
    end
    
    % save centroid positions
    dlmwrite(fullfile(dirpath,'pixelPositions.txt'),pixelPos,'delimiter','\t','newline','pc');
    dlmwrite(fullfile(dirpath,'parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
    dlmwrite(fullfile(dirpath,'fits.txt'),fits,'delimiter','\t','newline','pc');
    toc
end

