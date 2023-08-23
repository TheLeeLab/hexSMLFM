classdef MicroLensArray
    %MICROLENSARRAY Object contains parameters for microlens array used in
    %Fourier Light Field Microscopes.
    
    properties
        latticeType % 'square' or 'hexagonal'
        mlaCentre % spatial (x,y) position in terms of camera space
        focalLength % focal length microlenses
        pitch % along x-direction (parallel with optical table)
        rotation % rotation of the MLA in radians
        sizeOptic % pysical size of the array (assumes it's a square optic)
        lensCentres % list of lens centres in lattice spacing units
    end
    
    methods
        
        function obj = MicroLensArray(latticeType,focalLength,pitch,mlaCentre,sizeOptic)
            %MICROLENSARRAY Construct an instance of this class
            obj.latticeType = latticeType;
            obj.focalLength = focalLength;
            obj.pitch = pitch;
            obj.rotation = 0;
            obj.sizeOptic = sizeOptic;
            obj.mlaCentre = mlaCentre;
            obj.lensCentres = generateLattice(obj);
        end
        
        function lensCentres = generateLattice(obj)
            %GENERATELATTICE Generate the coordinates of the microlens
            %centers in the array first in integer lattice spacing then
            %convert to real space camera micrometre units.
            
            switch obj.latticeType
                case 'square'
                    % Generate centers of square lenses, pretending pitch = 1
                    % Generate more than needed, because gets cropped when rotated
                    N = obj.sizeOptic/obj.pitch;
                    x = -floor(N/2):ceil(N/2);
                    [x,y] = meshgrid(x,x);
                                      
                case 'hexagonal'
                    % Generate centers of hexagonal lenses, pretending pitch = 1
                    % Generate more than needed, because gets cropped when rotated
                    N = obj.sizeOptic/obj.pitch;
                    x = -floor(N/2):ceil(N/2);
                    [x,y] = meshgrid(x,x);
                    dx = repmat([0.5 0],[length(x),ceil(length(x)/2)]);
                    dx = dx(1:length(x),1:length(x));
                    x1 = x + dx'+0.5;
                    y2 = sqrt(3)/2*y;
                    x = y2;
                    y = x1;

                otherwise
                    disp(string(obj.latticeType)+' is not a valid MLA type.')
            end
            
            % keep in abstract units
            lensCentres = [x(:) y(:)];
                                                            
        end
        
        function obj = translateLattice(obj,dx,dy)
            %TRANSLATELATTICE Translate the coordinates of the microlens
            %centers by a shift in x and y
            obj.lensCenters(:,1) = obj.lensCenters(:,1)+dx;
            obj.lensCenters(:,2) = obj.lensCenters(:,2)+dy;

        end
        
        function obj = rotateLattice(obj,theta)
            %TRANSLATELATTICE Rotates the coordinates of the microlens
            %centers around the origin
            x = obj.lensCentres(:,1);
            y = obj.lensCentres(:,2);
            x = x.*cos(theta) - y*sin(theta);
            y = x.*sin(theta) + y*cos(theta);
            obj.lensCentres = [x y];
        end
        
    end
end

