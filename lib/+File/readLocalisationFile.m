% Read a localisation file (csv, txt, possible other) and return certain columns.
%
% @param filepath....: string, full path to the localisation file
% @param id_frame....: int, column containing frame in the localisations file
% @param id_x........: int, column containing x-coordinates
% @param id_y........: int, column containing y-coordinates
% @param id_sigma_x..: int, column containing sigma_x
% @param id_sigma_y..: int, column containing sigma_y
% @param id_int......: int, column containing intensity
% @param id_bkgnd....: int, column containing background
% @param id_precision: int, column containing precision (or uncertainty)
%
% @output locs: nx8 array, useful columns from localisation dataset

function locs = readLocalisationFile(filepath,type,pixel_size_sample)
scale = 1;
switch type
    case 'Picasso'
        id_frame = 2;
        id_x     = 3;
        id_y     = 4;
        id_sigx  = 5;
        id_sigy  = 5;
        id_int   = 6;
        id_bkgnd = 7;
        id_prec  = 9;
        scale= 1000;
    case 'Thunderstorm'
        id_frame = 1; % column containing frame
        id_x     = 2; % column containing x coordinates (in nm)
        id_y     = 3; % column containing y coordinates (in nm)
        id_sigx  = 4; % column containing sigma_x (not used for fitting)
        id_sigy  = 4; % column containing sigma_y (not used for fitting)
        id_int   = 5; % column containing intensity (not used for fitting)
        id_bkgnd = 6; % column containing background (not used for fitting)
        id_prec  = 7; % column containing precision or uncertainty (not used for fitting)
        scale = 1000;
    case 'Peakfit'
        id_frame = 1; % column containing frame
        id_x     = 10; % column containing x coordinates (in pixels)
        id_y     = 11; % column containing y coordinates (in pixels)
        id_sigx  = 13; % column containing sigma_x (not used for fitting)
        id_sigy  = 13; % column containing sigma_y (not used for fitting)
        id_int   = 9; % column containing intensity (not used for fitting)
        id_bkgnd = 8; % column containing background (not used for fitting)
        id_prec  = 14; % column containing precision or uncertainty (not used for fitting)
        scale = 1/pixel_size_sample; % adjust for pixel size
end

data = importdata(filepath);
data = data.data;
locs = nan(size(data,1),8);
locs(:,1) = data(:,id_frame);
locs(:,2) = data(:,id_x)/scale; % in microns
locs(:,3) = data(:,id_y)/scale; % in microns
locs(:,4) = data(:,id_sigx)/scale; % in microns
locs(:,5) = data(:,id_sigy)/scale; % in microns
locs(:,6) = data(:,id_int);
locs(:,7) = data(:,id_bkgnd); 
locs(:,8) = data(:,id_prec)/scale; % in microns

end
