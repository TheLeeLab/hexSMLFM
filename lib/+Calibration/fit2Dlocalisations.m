function locs = fit2Dlocalisations(x,y,img,w,method)
% FIT2DLOCALISATIONS  

% create meshgrid for centroid calculation in small 2*w+1 square roi
x_mesh = -w:w;
[x_mesh,y_mesh] = meshgrid(x_mesh,x_mesh);
y_mesh = -y_mesh;

locs = {};
switch method
    
    case 'centroid'
        for i=1:length(x)
            try
                roi = img(x(i)-w:x(i)+w,y(i)-w:y(i)+w); % crop square roi around localisation
                intensity = sum(roi(:));
                locs(i).intensity = intensity;
                locs(i).x = x(i) + sum(x_mesh.*roi,'all')/intensity;
                locs(i).y = y(i) + sum(y_mesh.*roi,'all')/intensity;
            catch % localisation was too close to edge to crop an roi, skip it
            end
        end
    
    case 'symmetric gaussian'
        for i=1:length(x)
            try
                roi = img(x(i)-w:x(i)+w,y(i)-w:y(i)+w); % crop square roi around localisation
                intensity = sum(roi(:));
                [params,~,~] = fit2DGaussian(roi,w,'symmetric');
                locs(i).intensity = intensity;
                locs(i).amplitude = params.amplitude;
                locs(i).x = x(i) + params.x0;
                locs(i).y = y(i) + params.y0;
                locs(i).sigma = params.sigma;
            catch % localisation was too close to edge to crop an roi, skip it
            end
        end
    
    case 'asymmetric gaussian'
        for i=1:length(x)
            try
                roi = img(x(i)-w:x(i)+w,y(i)-w:y(i)+w); % crop square roi around localisation
                intensity = sum(roi(:));
                [params,~,~] = fit2DGaussian(roi,w,'asymmetric');
                locs(i).intensity = intensity;
                locs(i).amplitude = params.amplitude;
                locs(i).x = x(i) + params.x0;
                locs(i).y = y(i) + params.y0;
                locs(i).sigmax = params.sigmax;
                locs(i).sigmay = params.sigmay;
            catch % localisation was too close to edge to crop an roi, skip it
            end
        end
        
    case 'rotated asymmetric gaussian'
        for i=1:length(x)
            try
                roi = img(x(i)-w:x(i)+w,y(i)-w:y(i)+w); % crop square roi around localisation
                intensity = sum(roi(:));
                [params,~,~] = fit2DGaussian(roi,w,'rotated asymmetric');
                locs(i).intensity = intensity;
                locs(i).amplitude = params.amplitude;
                locs(i).x = x(i) + params.x0;
                locs(i).y = y(i) + params.y0;
                locs(i).sigmax = params.sigmax;
                locs(i).sigmay = params.sigmay;
                locs(i).theta = params.theta;
            catch % localisation was too close to edge to crop an roi, skip it
            end
        end

    otherwise
        fprintf('Unexpected value for "method" in function "fit2Dlocalisations".')
end
end

function [params,resnorm,residual] = fit2DGaussian(img,w,type)
% FITSYMMETRIC2DGAUSSIAN fits a symmetric 2d gaussian to a square image.
% input:
%   img  - square array of size [2*w+1, 2*w+1], the image
%   w    - int, as defined above (could be calculated each time, but as it usually is
%               a constant variable it can just be supplied to the function)
%   type - string, 'symmetric', 'asymmetric' or 'rotated asymmetric' type
%               of 2d gaussian
% output:
%   params - ...
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

switch type
    
    case 'symmetric'
        % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigma]
        lb = [0,-w,-w,0];
        ub = [realmax('double'),w,w,w^2];
        x0 = [100,0,0,1]; % initial values
        options = optimset('Display','off'); % don't flood command window with messages
        [x,resnorm,residual,~] = lsqcurvefit(@symmetric2DGaussian,x0,xdata,img,lb,ub,options);
        params.amplitude = x(1);
        params.x0        = x(2);
        params.y0        = x(3);
        params.sigma     = x(4);
        
    case 'asymmetric'
        % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigmax,sigmay]
        lb = [0,-w,-w,0,0];
        ub = [realmax('double'),w,w,w^2,w^2];
        x0 = [100,0,0,1,1]; % initial values
        options = optimset('Display','off'); % don't flood command window with messages
        [x,resnorm,residual,~] = lsqcurvefit(@asymmetric2DGaussian,x0,xdata,img,lb,ub,options);
        params.amplitude = x(1);
        params.x0        = x(2);
        params.y0        = x(3);
        params.sigmax    = x(4);
        params.sigmay    = x(5);
        
    case 'rotated asymmetric'
        % set lower and upper bounds and initial values of x = [amplitude,xo,yo,sigmax,sigmay,theta]
        lb = [0,-w,-w,0,0,-pi/4];
        ub = [realmax('double'),w,w,w^2,w^2,pi/4];
        x0 = [100,0,0,1,1,0]; % initial values
        options = optimset('Display','off'); % don't flood command window with messages
        [x,resnorm,residual,~] = lsqcurvefit(@rotatedAsymmetric2DGaussian,x0,xdata,img,lb,ub,options);
        params.amplitude = x(1);
        params.x0        = x(2);
        params.y0        = x(3);
        params.sigmax    = x(4);
        params.sigmay    = x(5);
        params.theta     = x(6);
        
    otherwise
        disp('Unexpected value for "type" in function "fit2DGaussian".')
end

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

function F = asymmetric2DGaussian(x,xdata)
% ASYMMETRIC2DGAUSSIAN function describing an asymmetric 2d gaussian used as
% input for lsqcurvefit during fitting.
% input:
%   x     - 1x5 array containing parameters [amplitude, x0, y0, sigmax, sigmay]
%   xdata - NxNx2 matrix, containing spatial x and y values at which to
%           calculate the function
% output:
%   F - NxN matrix, the calcualted 2d gaussian
F = x(1)*exp(-(  ((xdata(:,:,1)-x(2)).^2)/(2*x(4).^2) + ...
                 ((xdata(:,:,2)-x(3)).^2)/(2*x(5).^2)  ));
end

function F = rotatedAsymmetric2DGaussian(x,xdata)
% ROTATEDASYMMETRIC2DGAUSSIAN function describing a rotated asymmetric 2d
% gaussian used as input for lsqcurvefit during fitting.
% input:
%   x     - 1x6 array containing parameters [amplitude, x0, y0, sigmax, sigmay, theta]
%   xdata - NxNx2 matrix, containing spatial x and y values at which to
%           calculate the function
% output:
%   F - NxN matrix, the calcualted 2d gaussian

% rotate x and y coordinate space by angle theta = x(6) around origin
xdatarot(:,:,1) = xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
xdatarot(:,:,2) = xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));

% rotate centroid position by angle theta = x(6) around origin
x0rot = x(2)*cos(x(6)) - x(3)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(3)*cos(x(6));

% asymmetric 2d gaussian but using rotated coordinates
F = x(1)*exp(-((xdatarot(:,:,1)-x0rot).^2/(2*x(4)^2) + ...
               (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2)));
end
