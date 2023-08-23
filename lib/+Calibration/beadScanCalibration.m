clear all
close all
clc

directory = 'C:\Users\ezrab\Desktop\githubRepos\smlfm_caic\smlfm\lib\+Calibration\testCalibrationData';
filename = 'calib_bead_2deg.tif';
filepath = fullfile(directory,filename);

% knowledge of focus position each frame
zRange = (-2000:50:2000)*1e-9; % z positions z stack calibration
reps    = 1; % repeat acquisitions per z-position
% generate coordinate list
z = zeros(length(zRange)*reps,1);
frame = 1;
for i=1:length(zRange)
    for j=1:reps
        z(frame) = zRange(i);
        frame = frame + 1;
    end
end

% localisation parameters
method = 'centroid'; % 'centroid', 'symmetric gaussian'
sigma1 = 2;
sigma2 = 3;
DoGmin = 20;
w = 3;


%% Background estimation

fprintf('Estimating background...')
N = 100; % number of frames used to estimate the background
k_med = 100; % kernelsize of median filter used in background estimation
bkgnd = Calibration.estimateBackground(filepath,N,k_med);
fprintf(' Done!\n')


%% Localisation

% set up pointer for reading frames one by one
t = Tiff(filepath,'r');

% initialize
frame = 0; id = 0;
localisations = {};

while true % read frames
    if ~mod(frame,10); fprintf('Processing frame %d\n',frame); end
    frame = frame + 1;
    
    img = double(t.read()) - bkgnd;
    
    % find local maxima (localisation candidates)
    [x,y] = Calibration.find2DlocalisationCandidates(img,sigma1,sigma2,DoGmin);
    
    % get better localisation for all candidates in this frame
    localisationsCurrentFrame = Calibration.fit2Dlocalisations(x,y,img,w,method);
    
    % add frame and z column
    for i=1:length(localisationsCurrentFrame)
        localisationsCurrentFrame(i).frame = frame;
        localisationsCurrentFrame(i).z = z(frame);
    end
    
    % append results to localisations
    if frame == 1; localisations = localisationsCurrentFrame;
    else; localisations = [localisations localisationsCurrentFrame];
    end
    
    % move pointer to next frame unless already at the last frame
    if t.lastDirectory(); break; else; t.nextDirectory(); end

end
close(t)


frame = [localisations.frame]';
x = [localisations.x]';
y = [localisations.y]';
int = [localisations.intensity]';

figure;
scatter(x,y,100,int,'.'); colorbar;
axis equal; grid on


%% Group localisation for each bead using first frame as seed

uniqueFrames = unique(frame);
seed_x = x(frame==1);
seed_y = y(frame==1);
numSeeds = length(seed_x);

% initialise arrays to contain temporally grouped localisations
x_grouped = nan(length(uniqueFrames),numSeeds);
y_grouped = nan(length(uniqueFrames),numSeeds);

% loop through all frames and link localisations to a seed
for i=1:length(uniqueFrames)
    x_i = x(frame==i);
    y_i = y(frame==i);
    for j=1:numSeeds
        dist = (seed_x(j) - x_i).^2 + (seed_y(j) - y_i).^2;
        [minDist,id] = min(dist);
        x_grouped(i,j) = x_i(id);
        y_grouped(i,j) = y_i(id);
        seed_x(j) = x_i(id); % update seed
        seed_y(j) = y_i(id); % update seed
    end
end

% Check successful temporal grouping
figure;
for i=1:numSeeds
    scatter(x_grouped(:,i),y_grouped(:,i),'.'); hold on
end
axis equal; grid on
title('Temporally grouped localisations')


%% Regression on each group

%% Group by clustering on slope

%% Get an estimate of the position of the central lens by crossing lines

%% Rotation estimate by grouping localisations