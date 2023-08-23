function [img] = locs2image(locs,pixelSize)
%LOCS2IMAGE This function takes in a list of xy localisations (two columns)
%and outputs an image. Localisations should be same unit as pixelSize
xmin = min(locs(:,1));
xmax = max(locs(:,1));
ymin = min(locs(:,2));
ymax = max(locs(:,2));

xrange = xmin:pixelSize:xmax;
yrange = ymin:pixelSize:ymax;
[img,~,~] = histcounts2(locs(:,1),locs(:,2),xrange,yrange);

end