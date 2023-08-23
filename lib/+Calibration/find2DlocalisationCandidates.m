function [x,y] = find2DlocalisationCandidates(img,sigma1,sigma2,DoGmin)
% FIND2DLOCALISATIONCANDIDATES find local maxima in an image after a
% difference-of-gaussian (DoG) filter.
% 
% DoG = (img guassian blurred with sigma1) - (img guassian blurred with sigma2)
% 
% input:
%   img    - array, the image
%   sigma1 - double, standard deviation of the first gaussian filter
%   sigma1 - double, standard deviation of the second gaussian filter
%   DoGmin - double, ignore all maxima in regions where DoG < DoGmin
% output:
%   x - nx1 array, the x-coordinates of identified local maxima (pixel resolution)
%   y - nx1 array, the y-coordinates of identified local maxima (pixel resolution)

% pre-processing filter
imgMed = medfilt2(img,[3,3],'symmetric');
G1 = imgaussfilt(imgMed,sigma1);
G2 = imgaussfilt(imgMed,sigma2);
DoG = G1 - G2;
DoG(DoG < DoGmin) = 0;

% find local maxima (localisation candidates)
BW = imregionalmax(DoG);
[x,y] = find(BW); % get coordinates of local maxima

end