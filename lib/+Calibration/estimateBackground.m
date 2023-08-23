function bkgnd = estimateBackground(filepath,N,k_med)
% ESTIMATEBACKGROUND estimate spatially varrying background from an stack
% using the first N frames and reading one frame at a time into memory.
%
% input:
%   filepath - string, full path to the image stack
%   N        - int, number of frames used for estimating the background
%   k_med    - int, kernel size of the median filter applied to the average
%              of the first N frames of the stack. The kernel should be
%              much larger than non-vbackground features in the images and
%              smaller than any features in the background so they don't
%              get smoothed out.
% output:
%   bkgnd - array, estimate of the background

t = Tiff(filepath,'r'); % set up pointer for reading frames one by one
avgProjection = double(t.read()); % read first frame
counter = 1; % counter to keep track of which frame we are on
while true && (counter <= N)
    % move pointer to next frame unless we were already at the last frame
    if t.lastDirectory(); break; else; t.nextDirectory(); end
    avgProjection = avgProjection + double(t.read());
    counter = counter + 1;
end
close(t)
avgProjection = avgProjection/counter;

% large kernel median filter to get rid of any dense structures that are
% not background
bkgnd = medfilt2(avgProjection,[k_med,k_med],'symmetric');

end