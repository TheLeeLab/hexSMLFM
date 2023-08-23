classdef LightFieldLocalisations
    %LIGHTFIELDLOCALISATIONS A class that contains light field localisation
    %data. The key property is 'locs' and contains light field formated
    %localisation data. This is in format:
    % [frame u v x y dx dy photons background precision alpha]
    % spatial units of micrometres and (u,v) in normalised pupil coords.
    
    properties
        locs
        filteredLocs
        minFrame
        maxFrame
        nLocs
        nFilteredLocs
        microLensArray
        fourierLFM
    end
    
    methods
        function obj = LightFieldLocalisations(locs_2d,microLensArray,fourierLFM)
            obj.nLocs = size(locs_2d,1);
            obj.fourierLFM = fourierLFM;
            obj.nFilteredLocs = size(locs_2d,1);
            obj.microLensArray = microLensArray;
            locs = zeros(obj.nLocs,10);
            locs(:,1) = locs_2d(:,1);
            locs(:,4:10) = locs_2d(:,2:8);
            obj.locs = locs;
            obj.minFrame = min(locs(:,1));
            obj.maxFrame = max(locs(:,1));
            obj.locs(:,2:3) = initialiseUV(obj);
            obj.filteredLocs = obj.locs;
        end
        
        
        function obj = filterSpotSize(obj,filterRange)
            %filterSpotSize: filters on sigma of localisations
            obj = filter(obj,filterRange,6);
        end
        
        function obj = filterPhotons(obj,filterRange)
            %filterSpotSize: filters on number of photons in spot
            obj = filter(obj,filterRange,8);
        end
        
        function obj = filterRho(obj,filterRange)
            %filterSpotSize: filters on number of photons in spot
            rho = sqrt(obj.filteredLocs(:,2).^2+obj.filteredLocs(:,3).^2);
            remove = rho<filterRange(1) | rho>filterRange(2);
            obj.filteredLocs(remove,:) = [];
            obj.nFilteredLocs = size(obj.filteredLocs,1);
        end
        
        function obj = filterLargestRhoLenses(obj,N)
            %remove largest N furthest our microlenses
            rho = sqrt(obj.filteredLocs(:,2).^2+obj.filteredLocs(:,3).^2);
            rhos = sort(unique(rho),'descend');
            rhoMax = rhos(N);
            remove = rho>rhoMax;
            obj.filteredLocs(remove,:) = [];
            obj.nFilteredLocs = size(obj.filteredLocs,1);
        end
        
        function obj = resetFilteredLocs(obj)
            obj.filteredLocs = obj.locs;
            obj.nFilteredLocs = obj.nLocs;
        end
        
        function obj = rotateUV(obj,theta)
            locs_temp = obj.filteredLocs;
            locs_temp(:,2) = (cos(theta)*locs_temp(:,2)+sin(theta)*locs_temp(:,3));
            locs_temp(:,3) = (cos(theta)*locs_temp(:,3)-sin(theta)*locs_temp(:,2));
            obj.filteredLocs = locs_temp;
        end
        
        function obj = setAlpha(obj,model)
            u = obj.filteredLocs(:,2);
            v = obj.filteredLocs(:,3);
            NA = obj.fourierLFM.NA;
            n = obj.fourierLFM.nMedium;
            switch model
                case 'sphere'
                    rhosq = sqrt(u.^2+v.^2);
                    alphau = -NA.*u./(n*sqrt(1-rhosq.*NA^2/n^2));
                    alphav = -NA.*v./(n*sqrt(1-rhosq.*NA^2/n^2));
                    alpha = [alphau alphav];
                    obj.filteredLocs(:,11:12) = alpha;
                case 'integrateSphere'
                    u_scaling = 2/obj.fourierLFM.numLensesInBfp;
                    [alphau, alphav] = Fitting.phase_average_sphere(u,v,u_scaling,NA,n,10);
                    obj.filteredLocs(:,11:12) = [alphau alphav];
                case 'linear'
                    alphau = u;
                    alphav = v;
                    alpha = [alphau alphav];
                    obj.filteredLocs(:,11:12) = alpha;
            end
        end
        
        function obj = correctUV(obj,corr)
            nViews = size(corr,1);
            locs = obj.filteredLocs;
            for i = 1:nViews
                ids = locs(:,2)==corr(i,1) & locs(:,3)==corr(i,2);
                locs(ids,4) = locs(ids,4)-corr(i,3);
                locs(ids,5) = locs(ids,5)-corr(i,4);
            end
            obj.filteredLocs = locs;
        end
        
        function h = plotUVs(obj)
            h = scatter(obj.filteredLocs(:,2),obj.filteredLocs(:,3),[],'.');
        end
        
        function h = plotXYs(obj,ind)
            h = scatter(obj.filteredLocs(:,4),obj.filteredLocs(:,5),10,obj.filteredLocs(:,ind),'.');
        end
        
    end
    
    
    
    methods (Access=private)
        
        function obj = filter(obj,filterRange,ind)
            remove = obj.filteredLocs(:,ind)<filterRange(1) | obj.filteredLocs(:,ind)>filterRange(2);
            obj.filteredLocs(remove,:) = [];
            obj.nFilteredLocs = size(obj.filteredLocs,1);
        end
        
        function uvs = initialiseUV(obj)
            rhoScaling = obj.fourierLFM.rhoScaling;
            mag = obj.fourierLFM.magnification;
            pos = obj.locs(:,4:5);
            pos(:,1) = pos(:,1)-obj.microLensArray.mlaCentre(1)/mag;
            pos(:,2) = pos(:,2)-obj.microLensArray.mlaCentre(2)/mag;
            pos = pos*rhoScaling;
            uvScale = 2/obj.fourierLFM.numLensesInBfp;
            lensCentres = obj.microLensArray.lensCentres;
            lensIndex = knnsearch(lensCentres*uvScale,pos);
            uvs = obj.microLensArray.lensCentres(lensIndex,:)*uvScale;
        end
        
        
        
    end
end

