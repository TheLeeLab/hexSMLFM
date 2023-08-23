% Assign localisation to microlenses and determine (u,v) for each.
% 
% @param locs........: nx8 array, localisation data where the columns contain
%                     frame,x,y,sigma_x,sigma_y,intensity,background,precision
% @param type_mla....: string, type of array ('square' or 'hexagonal')
% @param num_ulenses.: int, number of microlenses spanning the pupil in one
%                     dimension (assumes symmetry in x and y dimension)
% @param pix_per_lens: float, ratio of microlens pitch over physical camera pixel size
% 
% @output locs: nx10 array, identical to input 'locs' but with two extra
% columns -- u and v -- inserted after the first column. The new column
% order is: frame,u,v,x,y,sigma_x,sigma_y,intensity,background,precision

function locs = microlensAssignment(locs,type_mla,num_ulenses,pix_per_lens)

if strcmp(type_mla,'square')

    % Let user select the center of the (u,v) = (0,0) lens
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    scatter(locs(:,2),locs(:,3),1,'k.')
    axis equal; axis tight; set(gca,'YDir','reverse'); xlabel('x (nm)'); ylabel('y (nm)');
    set(0,'DefaultAxesTitleFontWeight','normal'); title({'Click on the center of the central lens.','Hit enter to continue.'})
    [x0,y0] = ginput;
    x0 = x0 - pix_per_lens/2;
    y0 = y0 - pix_per_lens/2;
    close(fig)
    
    % calculate u and v for each localisation
    u = floor((locs(:,2)-x0)/pix_per_lens);
    v = floor((locs(:,3)-y0)/pix_per_lens);
    
    % insert the u and v columns
    locs = [locs(:,1),u,v, locs(:,2:end)];
    
    % remove u and v that fall outside the microlenses
    locs(abs(locs(:,2))>ceil(num_ulenses/2),:) = [];
    locs(abs(locs(:,3))>ceil(num_ulenses/2),:) = [];
        
elseif strcmp(type_mla,'hexagonal')
    disp('Hexagonal arrays are not yet supported.'); return
else
    disp(['Sorry, "' mla_type '" is not a supported MLA.']); return
end

end