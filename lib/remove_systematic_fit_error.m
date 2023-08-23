function [corr] = remove_systematic_fit_error(total_fit,all_results,fit_terms,pix_per_lens,u_scale,min_view,defocus_limit,photon_thresh,NA,n,dx)
%FIND_AVERAGE_FIT_ERROR returns average fit error (sort of average aberration) across all fitted points
nLoc = size(total_fit,2);

views = unique(all_results(:,2:3),'rows');
n_views = size(views,1);
corr = zeros(n_views,5);
corr(:,1:2) = views;
for i = 1:nLoc
    inliers = total_fit(i).inliers;
    temp = inliers;
    error_fit = dist_fcn_LF2(total_fit(i).model,temp,fit_terms,pix_per_lens,NA,n,u_scale,dx);
    nin = size(inliers,1);
    photons = total_fit(i).photons;
    if nin>min_view && abs(total_fit(i).model(3))<defocus_limit && photons>photon_thresh
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

