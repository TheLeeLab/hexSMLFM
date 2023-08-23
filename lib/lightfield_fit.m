function [fitted_points,total_fit] = lightfield_fit(all_results,frame_range,pix_per_lens, dx, u_scaling,max_disparity,dz,angle_tol,fit_terms,fit_threshold,min_views,NA,n,numWorkers)
%LIGHTFIELD_FIT Summary of this function goes here

% total_fit holds the fit data in a struct so that individually chosen
% points can be inspected
min_frame = frame_range(1);
max_frame = frame_range(2);
total_fit=[];
% data0 just holds the fitted model in 2D array with row format
% (model,error_fit_model(1,2), error_fit_model(3), no_views_in_fit,no_photons,frame_no.) for fit terms = 0
% for higher fit_terms the model is extended to (x,y,a,b,c..., frame_o)
fitted_points = [];

parfor (frame=min_frame:max_frame,numWorkers)
    if ~mod(frame,1000); fprintf('Processing frame %d/%d\n',frame,max_frame); end
    
    % sort potential seeds by central view first and number of photons
    loc = all_results(all_results(:,1)==frame,:); 
    loc_cen = loc(abs(loc(:,2))<0.1 & abs(loc(:,3))<0.1,:);
    loc_out = loc(~(loc(:,2)==0 & loc(:,3)==0),:);
    [~,indSort] = sort(loc_cen(:,8),'descend');
    loc_cen = loc_cen(indSort,:);
    [~,indSort] = sort(loc_out(:,8),'descend');
    loc_out = loc_out(indSort,:);  
    loc = [loc_cen;loc_out];

    %loop over localisations until few left
    iter = 0;
    while size(loc,1)>(min_views) && iter<100
        iter = iter+1;
        seed = loc(1,:);
        %loc_near = loc_out;        
        [loc_fit,no_views,~] = find_valid_LF2_points(seed,loc,pix_per_lens,dz,max_disparity,angle_tol,u_scaling,dx,fit_terms,NA,n);
        % if there are now fewer than a certain no. views don't bother fitting, move on
        if no_views<min_views
                loc(1,:) = [];         
            continue
        end
        
       % fit to points
        loc_temp = loc_fit;
        [model, stdx, ~] = fit_function_LF2(loc_temp,fit_terms,pix_per_lens,NA,n,u_scaling,dx);
        % calculate error in model.
        dist = sqrt(stdx(1).^2+stdx(2).^2+stdx(3).^2);
        % if final fit distance is greater than n then don't fit, move on
        if dist>fit_threshold
            loc(1,:) = [];       
            continue
        end
        
        % if dist threshold met then remove fitted localisations from loc
        nFitted = size(loc_fit,1);
        for i = 1:nFitted
            [q,qidx] = ismember(loc_fit(i,:), loc, 'rows');
            if q
                loc(qidx,:) = [];
            end
        end
        
        %add to modelled and fitted points
        temp = {};
        if(~isempty(loc_fit))
            temp.frame=frame;
            sigx = loc_fit(:,6);
            sigy = loc_fit(:,7);
            amp = loc_fit(:,8); % number of photons
            intCounts = amp;
            %fit Model: (x,y,a,fit_error lateral, fit error axial, no views, no photon counts in fit, frame no.
            fitModel = [model;mean(stdx(1:2));stdx(3);no_views+1; sum(intCounts(:)); frame];
            temp.model= model;
            temp.inliers=loc_fit;
            temp.photons = sum(intCounts(:));
            temp.error = stdx;
            temp.u_scale = u_scaling;
            temp.ppl = pix_per_lens;
            temp.dx = dx;
            total_fit = [total_fit temp];
            fitted_points = [fitted_points;fitModel'];
        end      
    end
end


end