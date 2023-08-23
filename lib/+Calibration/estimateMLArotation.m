function mla_rotation = estimateMLArotation(xy,rotationWindow,method,flagshow)
% ESTIMATEMLAROTATION estimates the rotation that needs to be applied to the
% localisations to make a horizontal row of microlenses be aligned with
% the x-axis.
% 
% xy - array, two column localisation coordinates in units of pixels
% rotationWindow - double, a window of possible rotations the function will
% evaluate. It will evaluate [-rotationWindow/2, rotationWindow/2]
% method - string, 'radon' or 'autocorr'
% flagshow - bool, 1 to show some figures for troubleshooting
% 
% mla_rotation - the angle in radians that needs to be applied to correct
% for rotational misalignment of the microlens array using:
%    x_rot = cos(mla_rotation)*x + sin(mla_rotation)*y;
%    y_rot = cos(mla_rotation)*y - sin(mla_rotation)*x;
% assuming that x and y are defined using the center of the pupil as the
% origin

x = xy(:,1);
y = xy(:,2);

switch method
    case 'radon'
        mla_rotation = estimateMLArotation_radon(x,y,rotationWindow,flagshow);
    case 'autocorr'
        fprintf('"%s" is not yet implemented.\n',method)
    otherwise
        fprintf('"%s" is not a valid method for estimateMLArotation.m\n',method)
end


end

function mla_rotation = estimateMLArotation_radon(x,y,rotationWindow,flagshow)
% ESTIMATEMLAROTATION_RADON estimates the rotation that needs to applied to
% microlens array to make a horizontal row of microlenses be aligned with
% the x-axis.
% 
% x and y - arrays, localisation coordinates in units of pixels
% rotationWindow - double, a window of possible rotations the function will
% evaluate. It will evaluate [-rotationWindow/2, rotationWindow/2]
% flagshow - bool, 1 to show some figures for troubleshooting
% 
% mla_rotation - the angle in radians that needs to be applied to correct
% for rotational misalignment of the microlens array using:
%    x_rot = cos(mla_rotation)*x + sin(mla_rotation)*y;
%    y_rot = cos(mla_rotation)*y - sin(mla_rotation)*x;
% assuming that x and y are defined using the center of the pupil as the
% origin

% hardcoded parameters to limit number of inputs
pixsize_2dhist  = 0.5; % use bins of half a pixel size
gaussfilt_sigma = 1; % sigma of gaussian blur to image before radon transform
stepsize_radon  = 0.01; % degrees, stepsize evaluated orientations radon transform

% generate an image form localisations using a 2d histogram
Xedges = min(x):pixsize_2dhist:max(x);
Yedges = min(y):pixsize_2dhist:max(y);
[img,~,~] = histcounts2(x,y,Xedges,Yedges);

% apply gaussian filter with large kernel
img = imgaussfilt(img,gaussfilt_sigma);

% perform a radon transform for the provided rotationwindow
theta = -rotationWindow/2:stepsize_radon:rotationWindow/2;
[R,xp] = radon(img,theta);

% find maximum in radon transform to get mla rotation
[~,snd] = max(R(:));
[ij,ji] = ind2sub(size(R),snd);
mla_rotation = theta(ji) - theta(round(length(theta)/2));
mla_rotation = (pi/180)*mla_rotation;

if flagshow
    
    figure; 
    subplot(4,4,[1,2,5,6]); imshow(img,[]);
    title('Blurred image generated from localisations'); colorbar
    
    subplot(4,4,[9,10,13,14]);
    imshow(R,[],'XData',theta,'YData',xp)
    hold on; plot(theta(ji),xp(ij),'o')
    hold on; plot([theta(ji) theta(ji)],[min(xp) max(xp)])
    xlabel('\theta (degrees)')
    ylabel('x'''); axis square; axis on
    title('Radon transform')
    legend('Maximum','Estimated rotation')
    
    subplot(4,4,[3,4,7,8]);
    plot(theta,max(R,[],1)); xlabel('\theta'); ylabel('Maximum projection radon transform')
    hold on; plot([theta(ji) theta(ji)],[0 1.1*max(max(R,[],1))]); axis tight
    legend('Maximum proj radon trans','Maximum','location','southeast')
    
    % apply the rotation and plot points before and after rotation
    x_temp = x - mean(x);
    y_temp = y - mean(y);
    x_rot = cos(mla_rotation)*x_temp + sin(mla_rotation)*y_temp;
    y_rot = cos(mla_rotation)*y_temp - sin(mla_rotation)*x_temp;
    subplot(4,4,[11,12,15,16]);
    scatter(x_temp,y_temp,'.'); hold on
    scatter(x_rot,y_rot,'.')
    grid on; axis equal
    legend('before rotation','after rotation')
    
end

end
