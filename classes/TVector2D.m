classdef TVector2D
    %  Egbert
    %   basic vector object, used to represent object on edges or faces of
    %   a regular cartesian product grid
    
    properties
        % for now refer to E field, but can be any
        y  % 2D array of Ey field (Nx+1,Ny,Nz+1)
        z  % 2D array of Ez field (Nx+1,Ny+1,Nz)
        NdY; NdZ;	% number of fields in corresponding direction
        grid  %   pointer to grid
        type = 'edge' % 'edge' only
        Nyz
    end
    
    methods
        
        function obj = TVector2D(grid,~)
            if nargin >= 1
                if ~isa(grid,'TGrid2D')
                    error('First argument to TVector2D must be a compatible TGrid2D object')
                end
                obj.grid = grid;
                obj.NdY = [ grid.Ny grid.Nz+1];
                obj.NdZ = [ grid.Ny+1 grid.Nz];
                
                %   create a zero vector
                obj = obj.zero;
            end
        end
        %%
        function obj = zero(obj)
            % zero vector elements, assuming size is set
            obj.y = zeros(obj.NdY);
            obj.z = zeros(obj.NdZ);
            obj.Nyz = [prod(obj.NdY);prod(obj.NdZ)];
        end
        %%
        function result = length(obj)
            result = obj.Nyz(1)+obj.Nyz(2);
        end
        %%
        function  result = GetVector(obj)
            % creates standard matlab vector (1-D array) of E fields
            % the ordering is exactly the same as for Edges
            
            result = [reshape(obj.y,obj.Nyz(1),1);reshape(obj.z,obj.Nyz(2),1)];
        end
        %%
        function  obj = SetVector(obj,v)
            % Ey
            i1 = 1; i2 = obj.Nyz(1);
            obj.y = reshape(v(i1:i2), obj.NdY);
            %  Ez
            i1 = i2+1; i2 = i2 + obj.Nyz(2);
            obj.z = reshape(v(i1:i2), obj.NdZ);
        end;
        %%
        function  [ind_i,ind_b] = int_bdry_indices(obj)
            %   return indicies (in 1D array) of interior and boundary edges
            v = GetVector(obj.SetAllBdry);
            % interior edges/faces
            ind_i = find(v==0);
            if nargout == 2
                % boundary edges/faces
                ind_b = find(v);
            end
        end
        %%
        function obj1 = SetAllBdry(obj,c)
            %  set all boundary edges (faces) to specified constant c;
            %   defaults to one; all interior edges (faces) are set to 0
            if nargin == 1
                c = 1;
            end
            obj1 = obj.zero;
            obj1.y(1:end,[1 end]) = c;
            obj1.z([1 end],1:end) = c;
        end
        
        %%
        function result = consistent(obj1,obj2)
            result = strcmp(obj1.type,obj2.type);
            result = result && obj1.grid.Nx == obj2.grid.Nx && ...
                obj1.grid.Ny == obj2.grid.Ny && obj1.grid.Nz == obj2.grid.Nz;
        end
        %%  this is a prototype arithmetic object: dot products will also
        %         ulimately be useful
        function objout = plus(obj1,obj2)
            if ~consistent(obj1,obj2)
                error('objects are not same size or type; cannot add')
            end
            objout = obj1;
            objout.x = objout.x+obj2.x;
            objout.y = objout.y+obj2.y;
            objout.z = objout.z+obj2.z;
        end
        %%
        function objout = minus(obj1,obj2)
            if ~consistent(obj1,obj2)
                error('objects are not same size or type; cannot subtract')
            end
            objout = obj1;
            objout.x = objout.x-obj2.x;
            objout.y = objout.y-obj2.y;
            objout.z = objout.z-obj2.z;
        end
        %%
        function obj = interpFunc(obj1,location,yz)
            %  create a TVector object containiing weights needed for
            %  interpolation of yz component of obj1 to location
            obj = obj1.zero;
            switch yz
                case 'y'
                    yC = cumsum([obj1.grid.DDy]);
                    zC = cumsum([0;obj1.grid.Dz]);
                case 'z'
                    yC = cumsum([0;obj1.grid.Dy]);
                    zC = cumsum([obj1.grid.DDz]);
            end
            yC = yC;%-obj1.grid.origin(1);
            zC = zC-sum(obj1.grid.Dz(1:obj1.grid.Nza));%-obj1.grid.origin(2);
            %   find indicies of "mminimal corner"
            iy = find(location(1) > yC,1,'last');
            iz = find(location(2) > zC,1,'last');
            % find weights
            wy = (yC(iy+1)-location(1))/(yC(iy+1)-yC(iy));
            wz = (zC(iz+1)-location(2))/(zC(iz+1)-zC(iz));
            switch yz
                case 'y'
                    obj.y(iy,iz) = wy*wz;
                    obj.y(iy+1,iz) = (1-wy)*wz;
                    obj.y(iy,iz+1) = wy*(1-wz);
                    obj.y(iy+1,iz+1) = (1-wy)*(1-wz);
                case 'z'
                    obj.z(iy,iz) = wy*wz;
                    obj.z(iy+1,iz) = (1-wy)*wz;
                    obj.z(iy,iz+1) = wy*(1-wz);
                    obj.z(iy+1,iz+1) = (1-wy)*(1-wz);
            end            
        end
    end
end