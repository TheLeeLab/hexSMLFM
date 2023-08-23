% Authors: K O'Holleran & Ruth Sims
% Summary: example of estimating centre of pupil in (u,v) coordinates using
% a bead scan localisation file. Z values are not required - simply a
% single bead (very important it is only one) and a multitude of z values.
% These are projected onty (x,y) so z info is lost - calibration can even
% be turning the focus wheel slowly by hand.

name_file='../testdata/beadScanPeakFit.csv';  % name of file containing loc results
% the results table is a 2 D array in row format:
% (frame_no, u, v, x', y', sigma_x, sigma_y, signal, background, flag)'
% u and v and integers from -n:n. x' and y' are camera pixel localisations
% sigmas are similarly in pixels

% hard-coding the approximate top left corner of the u=0, v=0 lens
pix_per_lens = 1015/6.5; % 1015um pitch by 6.5um camera pixel
camera_gain = 1.9;
ppl = 2*1015/11;
x0 = 450-ppl/2;
y0 = 450-ppl/2;  
peak_fit = true;
if peak_fit
    pfData = csvread([name_file],1,1);
    [nLoc,~] = size(pfData);
    all_results = zeros(nLoc,10);   
    all_results(:,1) = pfData(:,1); %frame  
    all_results(:,2) = floor((pfData(:,2)-x0)/ppl);% u
    all_results(:,3) = floor((pfData(:,3)-y0)/ppl); % v
    all_results(:,4) = pfData(:,11); % x
    all_results(:,5) = pfData(:,12); % y
    all_results(:,6) = pfData(:,13); % dx
    all_results(:,7) = pfData(:,14); % dy
    all_results(:,8) = pfData(:,9); % signal  
    all_results(:,9) = pfData(:,8); % background
    all_results(:,10) = pfData(:,15); % precision
    all_results(abs(all_results(:,2))>2,:) = [];
    all_results(abs(all_results(:,3))>2,:) = [];
end

all_results(abs(all_results(:,2))>2,:) = [];
all_results(abs(all_results(:,3))>2,:) = [];

view_mask2 = [0 1 1 1 0;...
             1 1 1 1 1;...
             1 1 1 1 1;...
             1 1 1 1 1;
             0 1 1 1 0;];

% filter by view and by sigma_x, sigma_y (sugget looking at hist(all_results(:,6),100)
% to establish thresholds (and column 7)
sigma_x_range = [1,5];
sigma_y_range = [1,5];
all_results = filter_localisations(all_results,view_mask2,sigma_x_range,sigma_y_range);
all_results(all_results(:,8)<10000,:) = [];
 
[fv,fp] = find_image_plane(all_results,pix_per_lens);  
[all_results, corr_map] =  remove_aberration(all_results,fv,10,7,pix_per_lens);

figure(1)
scatter(all_results(:,4)-all_results(:,2)*pix_per_lens,all_results(:,5)-all_results(:,3)*pix_per_lens,200,all_results(:,1),'.')

%% Go through different lenses fitting straight lines
% Fit is performed to the group of all localisations in the view from all
% frames - should only be done with a single bead.
a = [];
data2 = [];
fits = [];
figure(3)
hold on
us = unique(all_results(:,2));
vs = unique(all_results(:,3));
for u=1:numel(us)
    for v =  1:numel(vs)
        data = all_results(all_results(:,2)==us(u) & all_results(:,3)==vs(v),:);
        if ~isempty(data)
            x = data(:,4);
            y = data(:,5);
            s = robustfit(x,y);
            plot(x,y,'o')
            pause(0.1)
            fits = [fits, [s;u;v]];        
        end
    end
end
x = 0:900;
for i = 1:size(fits,2)
    plot(x,fits(2,i).*x+fits(1,i),'--r','MarkerSize',1);
end
set(gca,'XLim',[0 900],'YLim',[0 900],'FontSize',18)
fits = fits';
ones_t = ones(size(fits,1),1);
centre = [ones_t -fits(:,2)]\fits(:,1);
scatter(centre(2),centre(1),400,'gx','LineWidth',3)

    
