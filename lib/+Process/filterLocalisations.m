% Filter localisation by view by using an array mask to remove views from 
% Fourier light field data.
% 
% localisations: array, list in form (frame, x', y', u, v, dx, dy, ampl, flag)
% view_mask: 2D binary array indicating which views to keep. Must be odd in size
% dx_range: [dx_min, dx_max], min and max pixel values for sigma x
% dy_range: [dy_min, dy_max], min and max pixel values for sigma y
%
% notes: x and y should be in camera pixels, u and v should be
% -num_views:+num_views

function localisations = filterLocalisations(localisations,...
    view_mask, dx_range, dy_range)

[m,n] = size(view_mask);

[iu, iv] = find(~view_mask);
iu = iu - (m+1)/2;
iv = iv - (n+1)/2;
for i = 1 :numel(iu)
   localisations(localisations(:,2)==iu(i) & localisations(:,3)==iv(i),:) = [];
end

localisations(localisations(:,6)<dx_range(1),:) = [];
localisations(localisations(:,6)>dx_range(2),:) = [];
localisations(localisations(:,7)<dy_range(1),:) = [];
localisations(localisations(:,7)>dy_range(2),:) = [];

end