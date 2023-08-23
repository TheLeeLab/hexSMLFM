classdef FourierLFM
    %SMLFM object contains microscope properties and convenience
    %methods for performing SMLFM.
    
    properties
        NA % numerical aperture of the objective
        objF % focal length of the objective lens (f_obj = M/f_tube)
        tubeLensF % focal length of the tube lens
        fourierLensF % focal length of the fourier lens
        mla % microLensArray object
        camPixelSize % size camera pixels
        nImmersion % refractive index of the cover glass
        nMedium % refractive index of the medium in which the molecule has been embedded
        bfpRadius % radius of the conjugate back focal plane
        numLensesInBfp % number of microlenses in the conjugate back focal plane
        pixPerLens % number of camera pixels per microlens
        magnification % overall magnification to camera plane
        imagePixelSize % size of pixel projected back to sample plane
        rhoScaling % scaling factor to convert microns in image plane to rho
    end
    
    methods
        
        function obj = FourierLFM(NA,objF,tubeLensF,fourierLensF,camPixelSize,nImmersion,nMedium,mla)
            %FOURIERLFM Construct an instance of this class
            obj.NA = NA;
            obj.objF = objF; % in mm
            obj.tubeLensF = tubeLensF; %in mm
            obj.fourierLensF = fourierLensF; % in mm
            obj.mla = mla;
            obj.camPixelSize = camPixelSize; % in microns
            obj.nImmersion = nImmersion;
            obj.nMedium = nMedium; 
            bfpRadius = objF*NA*(fourierLensF/tubeLensF)*1000; % radius of the back focal plane
            obj.bfpRadius = bfpRadius; % in micrometres 
            obj.numLensesInBfp =2*bfpRadius/(mla.pitch); % number of microlenses in the bfp
            obj.pixPerLens =mla.pitch/camPixelSize; % number of pixels per microlens
            obj.magnification = (tubeLensF/objF)*(mla.focalLength/fourierLensF);
            obj.imagePixelSize = camPixelSize/obj.magnification;
            obj.rhoScaling = obj.magnification/bfpRadius;
        end
        
        function [lensCenters,arrayPhase] = generatePhase(obj,wavelength)
            %GENERATEPHASE Calculates the microlens array phase
            
            [x,y] = obj.mla.generateLattice;
            lensCenters = [x,y];
            
            % normalize lens coordinates to the bfp radius
            x = x/obj.bfpRadius;
            y = y/obj.bfpRadius;
            
            % Calculate local field radius from lens centers in lattice
            xrange = linspace(-1,1,obj.numLensesInBfp*obj.pixPerLens);
            [lrf,~] = obj.localRadiusField(xrange,xrange,[x,y]);
            localRadius = obj.bfpRadius*lrf;
            
            % Convert local radius to phase
            k0 = 2*pi/wavelength;
            arrayPhase = -k0/(2*obj.mla.focalLength)*(localRadius.^2);
        end
        
        function fig = showPhase(obj,lensCenters,arrayPhase)
            %SHOWPHASE Displays a figure of the lens center coordinates and
            %phase of the microlens array
            
            x = lensCenters(:,1);
            y = lensCenters(:,2);
            
            fig = figure('position',[100 100 1600 500]);
            
            subplot(1,3,1)
            voronoi(x*1e3,y*1e3); hold on; axis equal;
            xlim([-obj.numLensesInBfp obj.numLensesInBfp]*obj.mla.pitch*1e3);
            ylim([-obj.numLensesInBfp obj.numLensesInBfp]*obj.mla.pitch*1e3);
            plot((obj.nMedium/obj.nImmersion)*obj.bfpRadius*1e3*cos(0:0.1:2.1*pi),(obj.nMedium/obj.nImmersion)*obj.bfpRadius*1e3*sin(0:0.1:2.1*pi),'r');
            xlabel('x (mm)'); ylabel('y (mm)'); title('Microlens array')
            legend('microlens centers','microlens edges (voronoi)','conjugate bfp outline','location','northeast')
            set(gca,'fontsize',10); set(0,'DefaultAxesTitleFontWeight','normal');
            
            subplot(1,3,2)
            imshow(flipud(arrayPhase),[]); colorbar; hold on; axis on
            xlabel('x (camera pixels)'); ylabel('y (camera pixels)'); title('Continuous phase')
            x_outline_bfp = (obj.nMedium/obj.nImmersion)*cos(0:0.1:2.1*pi)*size(arrayPhase,1)/2 + size(arrayPhase,1)/2;
            y_outline_bfp = (obj.nMedium/obj.nImmersion)*sin(0:0.1:2.1*pi)*size(arrayPhase,1)/2 + size(arrayPhase,1)/2;
            plot(x_outline_bfp,y_outline_bfp,'r','LineWidth',2);
            set(0,'DefaultAxesTitleFontWeight','normal');

            subplot(1,3,3)
            imshow(flipud(wrapTo2Pi(arrayPhase)),[]); colorbar; hold on; axis on
            xlabel('x (camera pixels)'); ylabel('y (camera pixels)'); title('Wrapped phase')
            plot(x_outline_bfp,y_outline_bfp,'r','LineWidth',2);
            set(0,'DefaultAxesTitleFontWeight','normal');
        end
        
    end
    
    
    methods (Static)
        
        function [lrf,ids] = localRadiusField(xrange,yrange,axes)
            %LOCALRADIUSFIELD Returns a matrix of real numbers representing distance
            %within xrange and yrange to nearest axis coordinate in axes
            %   xrange: vector -lower_x_lim:dx:upper_x_lim
            %   yrange: vector -lower_y_lim:dx:upper_y_lim
            %   axes: m x 2 matrix where m is number of local axis coordinates
            %
            % returns: lrf: matrix with values for radius to nearest axis
            %          ids: matrix with axes index linking to axes coords list
            % 
            % Author: Kevin O'Holleran

            [x,y] = meshgrid(xrange,yrange);
            nx = numel(xrange);
            ny = numel(yrange);
            num_axes = size(axes,1);

            lrf = inf(nx,ny);
            ids = zeros(nx,ny);
            for i = 1:num_axes
                xi = axes(i,1);
                yi = axes(i,2);
                r_temp = sqrt((x-xi).^2+(y-yi).^2);
                id_temp = r_temp<lrf;
                lrf(id_temp)= r_temp(id_temp);
                ids(id_temp) = i;
            end
        end
    end
    
end
