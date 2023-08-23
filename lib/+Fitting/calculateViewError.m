function [corr] = calculateViewError(fit_data,lfm,locs_2d,abberation_params)
%FIND_AVERAGE_FIT_ERROR returns average fit error (sort of average aberration) across all fitted points
nLoc = size(fit_data,2);
min_view = abberation_params.min_views;
photon_thresh = abberation_params.photon_thresh;
defocus_limit = abberation_params.axial_window;
rhoScaling = lfm.rhoScaling;
views = unique(locs_2d(:,2:3),'rows');
n_views = size(views,1);
corr = zeros(n_views,5);
corr(:,1:2) = views;
for i = 1:nLoc
    inliers = fit_data(i).inliers;
    temp = inliers;
    error_fit = Fitting.forwardModelError(fit_data(i).model,temp,rhoScaling);
    nin = size(inliers,1);
    photons = fit_data(i).photons;
    if nin>min_view && abs(fit_data(i).model(3))<defocus_limit && photons>photon_thresh
        for j = 1:n_views
            dum =  error_fit(inliers(:,2)==views(j,1) & inliers(:,3)==views(j,2),:);
            if ~isempty(dum)
                corr(j,3:4) =  corr(j,3:4)+dum;
                corr(j,5) =  corr(j,5)+1;
            end
        end
    end
end
corr(:,3:4) = corr(:,3:4)./[corr(:,5) corr(:,5)];

end

