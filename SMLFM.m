%% SMLFM.m
% Authors: Ruth R. Sims (1), Kevin O'Holleran (1), Sohaib Abdul Rehman (1), 
% Ezra Bruggeman (2) and Sam Daly (2)
% 1. Cambridge Advanced Imaging Centre, Downing Site, Cambridge, CB2 3DY, UK
% 2. Yusuf Hamied Department of Chemistry, Lensfield Road, University of 
%    Cambridge, Cambridge, CB2 1EW, UK

% Purpose: 
% It takes in 2D localisation data (x,y) captured on a fourier light field 
% microscope and turns it into 3D localisations (x,y,z). It does this by 
% first assigning (x,y) to (x,y,u,v) space and using microscope parameters 
% to calculate z position via parallax. 

% For more information, see: https://doi.org/10.1364/OPTICA.397172
% R. R. Sims, S. Abdul Rehman, M. O. Lenz, S. I. Benaissa, E. Bruggeman, 
% A. Clark, E. W. Sanders, A. Ponjavic, L. Muresan, S. F. Lee and 
% K. O'Holleran, Optica, 2020, 7, 1065. 

% Inputs: 
% Ensure parameters are correctly set in section 1 and use the gui to 
% select a 2D localizations file. An output folder will be generated here.
% The MLA can be reorientated on line 78 and 

% Outputs: 
% 3D localisations [x y z error(lateral) error(axial) views
% intensity frame] as a .csv, and a VISP-compatible .3d file 
% (https://www.nature.com/articles/nmeth.2566).

% Note: SMLFM.m has been tested on MATLAB R2022b and requires the 
% Statistics and Machine Learning Toolbox.

clear vars; close all; clc;

tic
addpath('lib')

%% 1. Set parameters

type_mla     = 'hexagonal'; % 'hexagonal' or 'square'
locs_format  = 'Peakfit'; % 'Peakfit', 'Thunderstorm', 'Picasso'
NA           = 1.27; % numerical aperture of objective
nImmersion   = 1.33; % immersion refractive index
nMedium      = 1.33; % specimen/medium refractive index
f_obj        = 200/60; % in mm
f_FourierLens= 175; % in mm
f_TubeLens   = 200; % in mm
f_MLA        = 175; % in mm
lens_pitch   = 2390; % in microns
pixel_size   = 16; % camera pixel size in microns
pixel_size_sample = pixel_size/(f_TubeLens/f_obj*f_MLA/f_FourierLens); % pixel size in sample space (microns)
sizeOptic = 10000; % size of MLA optic (in microns)
mlaRotation = 0;
mlaCentrePos = [0 0]*(f_TubeLens/f_obj)*(f_MLA/f_FourierLens); % in microns (adjust to suit data plotted on line 100)

z_calib = 1.534; % calibration between optical and physical z

save = 'no'; % save output files: 'yes' or 'no'

%% 2. Read localisation file

% select 2D localisation file (must be a .csv file)
[file,path] = uigetfile({'*.csv'}, 'Select 2D Localisations', 'C:\');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(file)]);
end
filepath = [path,file];

% if save = 'yes' on line 42 then a '3D Fitting' folder will be made and results saved here
outputFolder = [path '3D Fitting\'];

% read-in 2D localisation file
locs_2d = File.readLocalisationFile(filepath,locs_format,pixel_size_sample);

%% 3. Rotate x and y

theta = (30.8)*pi/180; % change to match the orientation of the MLA
x = locs_2d(:,2);
x = x-mean(x);
y = locs_2d(:,3);
y = y-mean(y);
locs_2d(:,2) = (x.*cos(theta) - y.*sin(theta));
locs_2d(:,3) = (x.*sin(theta) + y.*cos(theta));
scatter(locs_2d(:,2),locs_2d(:,3),[],'.');

%% 4. Initialise MLA, Microscope and LightFieldLocalisation objects

microLensArray = Classes.MicroLensArray(type_mla,f_MLA,lens_pitch,mlaCentrePos,sizeOptic);
lfm = Classes.FourierLFM(NA,f_obj,f_TubeLens,f_FourierLens,pixel_size,nImmersion,nMedium,microLensArray);
lfLocs = Classes.LightFieldLocalisations(locs_2d,microLensArray,lfm);

%% 5. Filter and rotate localisations

lfLocs = lfLocs.resetFilteredLocs;
lfLocs = lfLocs.filterRho([0 0.8]); % 0 0.6
%lfLocs = lfLocs.filterSpotSize([0.1 1]);

lfLocs.plotXYs(3) % break code here to check if microlenses are correctly assigned

lfLocs.plotUVs

% estimate MLA rotation
dTheta = 5; % degrees, will evaluate [-dTheta/2, dTheta/2]
% mla_rotation = -Calibration.estimateMLArotation(locs_2d(:,2:3),dTheta,'radon',0);
% fprintf('Estimated MLA rotation: %.3f degrees\n',mla_rotation*180/pi )


%% 6. Find system aberrations

fit_params = {};
fit_params.frame_range = [lfLocs.minFrame 1000];
fit_params.max_disparity = 5; % find locs from -5 to 5 um
fit_params.dist_search = 0.5;
fit_params.angle_tol = 2*pi/180;
fit_params.threshold = 1;
fit_params.min_views = 3; %5

mla_rotation = -0*pi/180;
lfLocs = lfLocs.resetFilteredLocs;
%lfLocs = lfLocs.filterRho([0 0.6]);
%lfLocs = lfLocs.filterSpotSize([0.1 1]);
lfLocs.plotXYs(3)

lfLocs = lfLocs.rotateUV(mla_rotation);
lfLocs = lfLocs.setAlpha('integrateSphere');

abberation_params = {};
abberation_params.axial_window = 1;
abberation_params.min_views = 3; 
abberation_params.photon_thresh = 1;

fprintf('Fitting first 1000 frames for aberration correction')
numWorkers = 8; 
[locs3D,fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs,lfm,fit_params,numWorkers);
correction = Fitting.calculateViewError(fit_data,lfm,lfLocs.filteredLocs,abberation_params);

lfLocs = lfLocs.correctUV(correction);
fprintf('Global abberation calculated (um): views listed below (u,v,dx,dy)')
correction(:,1:4)

%% 7. Fit full data set on corrected localisations

fit_params.frame_range = [lfLocs.minFrame lfLocs.maxFrame];
fit_params.max_disparity = 8;
fit_params.dz = 0.5;
fit_params.angle_tol = 1*pi/180;
fit_params.threshold = 0.3;
fit_params.min_views = 2;
fprintf('Fitting whole data set')
[locs3D,fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs,lfm,fit_params,numWorkers);
fprintf('Done!\n\n')
fprintf('Total number of 2d localisations used for fitting: %d\n',sum(locs3D(:,end-2)))
fprintf('Final number of 3d localisations: %d\n',size(locs3D,1))
toc

%% 7. Plotting

locs3D(:,3) = locs3D(:,3) * z_calib; % apply z calibration to z coordinates

x = locs3D(:,1);
y = locs3D(:,2);
z = locs3D(:,3);
n_views     = locs3D(:,end-2); % number of views used to fit the localization
lateral_err = 1000*locs3D(:,4); % fitting error in x and y
axial_err   = 1000*locs3D(:,5); % fitting error in z

keep = lateral_err<200 & n_views>3; % thresholds for plotting data

figure(1)
scatter3(x(keep),y(keep),z(keep),80,z(keep),'.')
xlabel('x \mum'); ylabel('y \mum'); zlabel('z \mum')
c = colorbar; c.Label.String = 'Z (\mu m)'; % colormap hot
axis equal
% set(gca,'XLim',[-10 10],'YLim',[-10 10],'ZLim',[-10 10],'FontSize',18)

figure(2)
subplot(1,3,1); histogram(lateral_err(keep),1:5:200);
xlabel('Lateral fit error (nm)'); ylabel('Occurence')
subplot(1,3,2); histogram(axial_err(keep),1:5:200);
xlabel('Axial fit error (nm)'); ylabel('Occurence')
subplot(1,3,3); histogram(locs3D(keep,end-1));
xlabel('Number of photons'); ylabel('Occurence')
set(gcf,'Position',[100 100 1200 300])

figure(3)
histogram2(locs3D(keep,end-1),axial_err(keep),'DisplayStyle','tile','ShowEmptyBins','on');
xlabel('Photons'); ylabel('Axial precision');

%% 8. Write the result

switch save % if save = 'yes' in line 57
    
    case 'yes'
    fprintf('Writing localisations to output files...\n')
    mkdir(outputFolder);
    VISP = [locs3D(:,1)*1000, locs3D(:,2)*1000, locs3D(:,3)*1000, locs3D(:,7), locs3D(:,8)];
    dlmwrite([outputFolder 'VISP_(' datestr(now,'dd-mm-yyyy_MM-HH') ').3d'],VISP ,'delimiter','\t','precision',7);
    writematrix(locs3D, [outputFolder 'locs3D.csv']);
    %writematrix(phoabvbkg, [outputFolder 'locSNR_(' datestr(now,'dd-mm-yyyy_MM-HH') ').csv']);

    fprintf('Saving plots...\n')
    saveas(figure(1), [outputFolder 'plot3D.fig']);
    saveas(figure(2), [outputFolder 'graphs.fig']);
    saveas(figure(3), [outputFolder 'histogram.fig']);
    
    case 'no'

end
