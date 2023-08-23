function [points,n_views,n_duplicates] = find_valid_LF2_points(seed,in_points,pix_per_lens, dis_tol, max_disparity, angle_tol,u_scaling,pix_dx,model_num,NA,n)
%FIND_VALID_LF2_POINTS Returns a subset of in_points which are candidates
%for being related to a single point - seed
% seed - single row if LF2 data
% in_points - array of LF2 data
% pix_per_lens - pixels per lens
% dis_tol = the window to accept points as belonging to a given z
% max_disparity - limit of possible disparity between adjacent views
% angle tol - tolerance for selecting with angles around relative uv

dz = 0.1;
sx = seed(1,4);
sy = seed(1,5);
su = seed(1,2);
sv = seed(1,3);

% first remove self view points
in_points(in_points(:,2)==su & in_points(:,3)==sv,:) = [];

% disparity from seed and relative u and v
dx = in_points(:,4) - sx;
dy = in_points(:,5) - sy;
du = in_points(:,2) - su;
dv = in_points(:,3) - sv;
r = sqrt(dx.^2+dy.^2);

% angle selection is relative to the uv index of the seed and dx dy on
% camera pixels
angles_xy = atan2(dy,dx);
angles_uv = atan2(dv,du);
% angle selection
angle_diff=min(abs(angles_xy-angles_uv),pi-abs(angles_xy-angles_uv));
angle_selection = abs(angles_xy-angles_uv)<angle_tol| abs(angles_xy-angles_uv)>(pi-angle_tol);
%in_points(~angle_selection,:) = [];

% rect selection
ortho_dist=r.*sin(abs(angles_xy-angles_uv));
rect_selection = abs(ortho_dist) < 0.22;
in_points(~rect_selection,:) = [];

% select against relative disparity
du = in_points(:,2) - su;
dv = in_points(:,3) - sv;
dx = in_points(:,4) - sx;
dy = in_points(:,5) - sy;
r = sqrt(dx.^2+dy.^2);
ruv = pix_per_lens*pix_dx*sqrt(du.^2+dv.^2); % lens pitch * uv length
rho = sqrt(du.^2+dv.^2);
disparity = (r-ruv)./rho;

zrange = -max_disparity:dz:max_disparity; 
num_points = zeros(numel(zrange),1);
% num_points will be in integration of the number of points in
% each disparity gradient range
for z = 1:numel(zrange)
    num_points(z) = sum(disparity>(zrange(z)-dis_tol) & disparity<=(zrange(z)+dis_tol));
end
[~,best_z_ind] = max(num_points);
best_z = zrange(best_z_ind);
% apply z-selection
z_selection = disparity>(best_z-dis_tol) & disparity<=(best_z+dis_tol);
in_points = in_points(z_selection,:);

% get rid of multiple points in single views
% This section looks for multiple localisation within a single view
% that made the selection - we only need one so pick one... this
% could be done better I'm sure! I take the one which has best fit. (not
% all possible combinations checked though, just as they are found)

% initialise set of localisatins to fit with the seed being
% used first
points = seed;
n_duplicates = 0;
n_views = 0;
if ~isempty(in_points)
    [uv, ~, iuv] = unique(in_points(:,2:3), 'rows');
    % generates indices for localisations belonging to a view
    indices = accumarray(iuv, find(iuv), [], @(rows){rows});
    n_views = length(indices);
    view_count = cellfun(@(x) length(x),indices);
    % extract views with only one localisation
    points = [points; in_points(ismember(iuv,find(view_count==1)),:)];
    
    indices(view_count==1) = [];
    nd = numel(indices);
    n_duplicates = nd;
    % If there is more than one point in a view take the one which fits the unique
    % points best
    if nd~=0
        for i = 1:nd
            view_ind = indices{i};
            n_view = size(view_ind,1);
            d = zeros(n_view,1);
            for j = 1:n_view
                temp1 = points;
                temp2 = in_points;
                model = fit_function_LF2([temp1;temp2(view_ind(j),:)],model_num,pix_per_lens,NA,n,u_scaling,pix_dx);
                dum = dist_fcn_LF2(model,[temp1;temp2(view_ind(j),:)],model_num,pix_per_lens,NA,n,u_scaling,pix_dx);
                d(j) = sqrt(sum(dum(:).^2));
            end
            [~,choose] = min(d);
            points = [points; in_points(view_ind(choose),:)];
        end
    end
    
    % if there are two views that share common u or v then reject the grouping
    uvar = var(points(:,2));
    vvar = var(points(:,3));
    if uvar<0.1 || vvar <0.1 && n_views == 1
        points = [];
        n_views = 0;
    end
    
    
end

end

