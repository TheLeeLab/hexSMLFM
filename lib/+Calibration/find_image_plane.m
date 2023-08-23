% Summmary: returns the frame number with:
% (1) minimum spot position variance after converting to local perspective 
% view coordinates (frame_No_var).
% (2) minimum median spot sigma (frame_No_width)
% Inputs: locs - localisation array in SMLFM format
% Inputs: pixPerLens
% NOTE: designed to work for bead scan data with single bead being scanned 
% in z.
%
% What is this useful for? Finding z = 0 (minimum positional variance).
% Finding sharpest imaged frame (minimum spot width). You would think these
% may be the same but they often aren't.

function [frame_No_var,frame_No_width]=find_image_plane(locs,pixPerLens)

std_x=[];
std_y=[];
mwx=[];
mwy=[];

frame_all=[];
for frame=1:max(locs(:,1))    
    data_frame=locs(locs(:,1)==frame,:);
    x=data_frame(:,4)-data_frame(:,2)*pixPerLens;
    y=data_frame(:,5)-data_frame(:,3)*pixPerLens;
    w_x=data_frame(:,6);
    w_y=data_frame(:,7);
    
    std_x=[std_x std(x)];
    std_y=[std_y std(y)];
    
    mwx=[mwx median(w_x)];
    mwy=[mwy median(w_y)];
    
    frame_all=[frame_all frame];
end

var_r=std_x.^2+std_y.^2;
width=sqrt(mwx.^2+mwy.^2);

[~,frame_No_var]=min(var_r);
[~,frame_No_width]=min(width);
end
