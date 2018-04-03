classdef TScalar2D
    %
    %   basic scalar object for 2D, derived from 3D version
    
    properties
        % for now refer to E field, but can be any
        v    %   2D array representing scalar field
        NdV	%  dimensions of array
        grid  %   pointer to grid
        type = 'node' % node, cell
        Nv  %   total length of array
    end
    
    methods
        
        function obj = TScalar2D(grid,type)
            if nargin >= 1
                if ~isa(grid,'TGrid2D')
                    error('First argument to TScalar2D must be a TGrid2D object')
                end
                obj.grid = grid;
                if nargin ==2
                    obj.type = lower(type);
                    %   defaults to node
                end
                switch obj.type
                    case 'node'
                        obj.NdV = [grid.Ny+1 grid.Nz+1];
                    case 'cell'
                        obj.NdV = [ grid.Ny grid.Nz];
                    otherwise
                        error('only node or cell allowed for property type')
                end
                %   create a zero vector
                obj.v = zeros(obj.NdV);
                obj.Nv = prod(obj.NdV);
            end
        end
        %%
        function obj = zero(obj)
            % zero vector elements, assuming size is set
            
            obj.v = zeros(obj.NdV);  
            obj.Nv = prod(obj.NdV);
        end
        %%
        function result = length(obj)
            result = obj.Nv;
        end
        %%
        function  result = GetVector(obj)
            % creates standard matlab vector (1-D array) of E fields
            % the ordering is exactly the same as for Edges
            
            result = reshape(obj.v,obj.Nv,1);
        end
        %%
        function  obj = SetVector(obj,v)
            obj.v = reshape(v, obj.NdV);
        end;
        %%
        function  [ind_i,ind_b] = int_bdry_indices(obj)
            %   return indicies (in 1D array) of interior and boundary edges
            switch obj.type
                case 'cell'
                    error('no such thing as boundary cells!')
                case 'node'
                    temp = GetVector(obj.SetAllBdry);
                    % inner nodes
                    ind_i = find(temp==0);
                    if nargout == 2
                        % boundary nodes
                        ind_b = find(temp);
                    end
            end
        end
        %%
        function obj1 = SetAllBdry(obj,c)
            %  set all boundary nodes to specified constant c;
            %   defaults to one; all interior nodes are set to 0
            %    this really only makes sense for nodes
            switch obj.type
                case 'cell'
                    error('no such thing as boundary cells!')
                case 'node'
                    if nargin == 1
                        c = 1;
                    end
                    obj1 = obj.zero;
                    obj1.v([1 end],:) = c;
                    obj1.v(:,[1 end]) = c;
            end
        end

  
        %%
        function result = consistent(obj1,obj2)
            result = strcmp(obj1.type,obj2.type);
            result = result &&  ...
                obj1.grid.Ny == obj2.grid.Ny && obj1.grid.Nz == obj2.grid.Nz;
        end
        %%  this is a prototype arithmetic object: dot products will also
        %         ulimately be useful
        function objout = plus(obj1,obj2)
            if ~consistent(obj1,obj2)
                error('objects are not same size or type; cannot add')
            end
            objout = obj1;
            objout.v = objout.v+obj2.v;
        end
        %%  
        function objout = minus(obj1,obj2)
            if ~consistent(obj1,obj2)
                error('objects are not same size or type; cannot subtract')
            end
            objout = obj1;
            objout.v = objout.v-obj2.v;
        end
        %%
        function vecOut = vecAvg(obj,type)
            %    average scalar object (cell type only) onto TVector object; of edge or
            %     node type;   boundary edges are left set to 0???
            
            
            if ~strcmp(obj.type,'cell')
                error('vecAvg only makes sense for scalar fields of type cell')
            end
            if nargin == 1
                type = 'node';
            end
            [DY,DZ] = ndgrid(obj.grid.Dy,obj.grid.Dz);
            vol = DY.*DZ;
            temp = obj.v.*vol;
            switch type
                case 'node'
                    vecOut = TScalar2D(obj.grid,type);
                    vecOut.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                        temp(1:end-1,2:end)+temp(2:end,2:end);
                    vecOut.v(2:end-1,2:end-1) = vecOut.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                        vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
                case 'edge'
                    vecOut = TVector2D(obj.grid);
                    vecOut.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
                    vecOut.y(:,2:end-1) = vecOut.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                        vol(:,2:end));
                    vecOut.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
                    vecOut.z(2:end-1,:) = vecOut.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                        vol(2:end,:));
                    % Zeqiu notes: Entries of 1 in the Del^2 operator matrix will multiply Field vector to
                    % determine the BCs, so averaging of model parameters on corresponding
                    % boundary nodes or edges are zeros (for isotropy or diagonal part of anisotrpy),
                    % for anisotropy, averaging on dual boundary domain is not zero
                    % While for 1D case, entries of boundary elements are assigned with 1 in
                    % the averageing/parameter mapping, not in the operator matrix.
                    
            end
        end
        %%
        function obj = interpFunc(obj1,location,type)
            %  create a TVector object containiing weights needed for
            %  interpolation of yz component of obj1 to location
            obj = obj1.zero;
            switch type
                case 'node'
                    yC = cumsum([0;obj1.grid.Dy]);
                    zC = cumsum([0;obj1.grid.Dz]);
                case 'cell'
                    yC = cumsum([obj1.grid.DDy]);
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
            switch type
                case 'node'
                    obj.v(iy,iz) = wy*wz;
                    obj.v(iy+1,iz) = (1-wy)*wz;
                    obj.v(iy,iz+1) = wy*(1-wz);
                    obj.v(iy+1,iz+1) = (1-wy)*(1-wz);
                case 'cell'
                    obj.v(iy,iz) = wy*wz;
                    obj.v(iy+1,iz) = (1-wy)*wz;
                    obj.v(iy,iz+1) = wy*(1-wz);
                    obj.v(iy+1,iz+1) = (1-wy)*(1-wz);
            end            
        end
    end
end
