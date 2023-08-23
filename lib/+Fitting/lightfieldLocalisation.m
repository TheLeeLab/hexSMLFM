function [fitted_points,total_fit] = lightfieldLocalisation(locs_2d,fourierLFM,fit_params,numWorkers)
%LIGHTFIELDLOCALISATION Based on lightfield_fit.m

% total_fit holds the fit data in a struct so that individually chosen
% points can be inspected
frame_range = fit_params.frame_range;
min_frame = frame_range(1);
max_frame = frame_range(2);
min_views = fit_params.min_views;
fit_threshold = fit_params.threshold;
rhoScaling = fourierLFM.rhoScaling;
total_fit=[];
% fitted_points just holds the fitted model in 2D array with row format
% (model,error_fit_model(1,2), error_fit_model(3), no_views_in_fit,no_photons,frame_no.) for fit terms = 0
fitted_points = [];

parfor (frame=min_frame:max_frame,numWorkers)
    if ~mod(frame,1000); fprintf('Processing frame %d/%d\n',frame,max_frame); end
    
    % sort potential seeds by central view first and number of photons
    candidates = locs_2d(locs_2d(:,1)==frame,:); 
    loc_cen = candidates(abs(candidates(:,2))<0.1 & abs(candidates(:,3))<0.1,:);
    loc_out = candidates(~(candidates(:,2)==0 & candidates(:,3)==0),:);
    [~,indSort] = sort(loc_cen(:,8),'descend');
    loc_cen = loc_cen(indSort,:);
    [~,indSort] = sort(loc_out(:,8),'descend');
    loc_out = loc_out(indSort,:);  
    candidates = [loc_cen;loc_out];

    %loop over localisations until few left
    iter = 0;
    while size(candidates,1)>(min_views) && iter<100
        iter = iter+1;
        seed = candidates(1,:);
        %loc_near = loc_out;        
        [loc_fit,no_views,~] = Fitting.groupLocalisations(seed,candidates,fit_params,fourierLFM);
        % if there are now fewer than a certain no. views don't bother fitting, move on
        if no_views<min_views
                candidates(1,:) = [];         
            continue
        end
        
       % fit to points
        loc_temp = loc_fit;
        [model, stdx, ~] = Fitting.backwardModel(loc_temp,rhoScaling);
        % calculate error in model.
        dist = sqrt(stdx(1).^2+stdx(2).^2+stdx(3).^2);
        % if final fit distance is greater than n then don't fit, move on
        if dist>fit_threshold
            candidates(1,:) = [];       
            continue
        end
        
        % if dist threshold met then remove fitted localisations from loc
        nFitted = size(loc_fit,1);
        for i = 1:nFitted
            [q,qidx] = ismember(loc_fit(i,:), candidates, 'rows');
            if q
                candidates(qidx,:) = [];
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
            total_fit = [total_fit temp];
            fitted_points = [fitted_points;fitModel'];
        end      
    end
end


end