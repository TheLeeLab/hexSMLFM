function [corrected, corr_map] = remove_aberration(localisations, frame, repeats,nlens,ppl)
%REMOVE_ABERRATION Takes input:
%   localisations: 2D array with rows - (frame_no,u,v, x, y, sigma_x, sigma_y, amplitude, bakground, flag?)
%   frame: first frame (of a set of repeats) to calculate aberration
%   repeats: how mamy frames were taken at each z position
%   nlens: number lenses across width of data set (must be odd)
%
% output:
%   corrected: 2D array in same format as localisations but with aberration removed
%   corr_map: 2D array in space of (u,v) with pixel magnitudes of
%   correction

% pull out subset of localisations for plane of interest denoted by frame
loc=localisations(localisations(:,1)>=frame & localisations(:,1)<=(frame+repeats-1),:);
 t = (nlens-1)/2;
% initialise and calculate correction map
corr_map = zeros(nlens,nlens,2);
for u=1:nlens
    for v=1:nlens
        temp=loc(loc(:,2)==u-3 & loc(:,3)==v-3,:);
        if ~isempty(temp)
            corr_map(u,v,1)=mean(temp(:,4))-ppl*temp(1,2);
            corr_map(u,v,2)=mean(temp(:,5))-ppl*temp(1,3);
        end
    end
end
mask = corr_map(:,:,1)~=0;
corr_map(:,:,1) = mask.*(corr_map(:,:,1)-corr_map(t,t,1));
corr_map(:,:,2) = mask.*(corr_map(:,:,2)-corr_map(t,t,2));


%correct full localisations list with corr_map
corrected = localisations;
for i = 1:size(localisations(:,1),1)
    corrected(i,4) = localisations(i,4)-corr_map(localisations(i,2)+t,localisations(i,3)+t,1);
    corrected(i,5) = localisations(i,5)-corr_map(localisations(i,2)+t,localisations(i,3)+t,2);
end


end

