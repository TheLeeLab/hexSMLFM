clc
clear all
close all

w = 15;

x0 = 6.6;
y0 = 4;
amplitude = 500;
sigma = 2;

% generate a fake dataset
x = -w:w;
[x,y] = meshgrid(x,x);
xdata = cat(3,x,y);
img = symmetric2DGaussian([amplitude,x0,y0,sigma],xdata);
img = poissrnd(img);

imshow(img,[])

[amplitude,x0,y0,sigma,resnorm,residual] = fitSymmetric2DGaussian(img,w);
hold on
plot(x0+w+1,y0+w+1,'rx')

function [amplitude,x0,y0,sigma,resnorm,residual] = fitSymmetric2DGaussian(img,w)
% FITSYMMETRIC2DGAUSSIAN fits a symmetric 2d gaussian to a square image.
% input:
%   img - square array of size [2*w+1, 2*w+1], the image
%   w   - int, as defined above (could be calculated each time, but as it usually is
%              a constant variable it can just be supplied to the function)
% output:
%   amplitude - double, the amplitude of the fitted gaussian
%   x0 and y0 - double, the x and y coordinate of the center of the
%               gaussian, relative to the center of the middle pixel which
%               is at (0,0)
%   sigma     - double, the standard deviation of the fitted gaussian
%   resnorm   - double, sum of the square of the residuals, i.e. sum(residual(:).^2)
%   residual  - array, the residual of the fit, same size as img

% get (x,y) meshgrid
x = -w:w;
[x,y] = meshgrid(x,x);
xdata = cat(3,x,y);

% define lower and upper bounds [amplitude,xo,yo,sigma]
lb = [0,-w,-w,0];
ub = [realmax('double'),w,w,w^2];
x0 = [100,0,0,1]; % initial values

% least squares fitting
options = optimset('Display','off'); % don't flood command window with messages
[x,resnorm,residual,~] = lsqcurvefit(@symmetric2DGaussian,x0,xdata,img,lb,ub,options);

amplitude = x(1);
x0        = x(2);
y0        = x(3);
sigma     = x(4);

end

function F = symmetric2DGaussian(x,xdata)
% SYMMETRIC2DGAUSSIAN function describing a symmetric 2d gaussian used as
% input for lsqcurvefit during fitting.
% input:
%   x     - 1x4 array containing parameters [amplitude, x0, y0, sigma]
%   xdata - NxNx2 matrix, containing spatial x and y values at which to
%           calculate the function
% output:
%   F - NxN matrix, the calcualted 2d gaussian
F = x(1)*exp(-(  (xdata(:,:,1)-x(2)).^2  + ...
                 (xdata(:,:,2)-x(3)).^2  )/(2*x(4).^2));
end

% function F = asymmetric2DGaussian(x,xdata)
% F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + ...
%                (xdata(:,:,2)-x(4)).^2/(2*x(5)^2)));
% end
% 
% function F = rotatedAsymmetric2DGaussian(x,xdata)
% F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + ...
%                (xdata(:,:,2)-x(4)).^2/(2*x(5)^2)));
% end
