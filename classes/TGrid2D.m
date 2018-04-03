classdef TGrid2D < handle
    % Egbert 2014
    %  I only forsee one 3D grid type, so no abstract superclass here
   
    properties
        
        Ny, Nz   %number of grid cells in y,z directions
        Nza   %     number of air layers
        Dy, Dz    %     cell dimensions: x,y,z-direction
        NEdges, NNodes, NCells
        
        
        NzEarth      %   added for consistency with other variants of grid objects
        %    these certainly should not be grid properties!
        ind_i, ind_b, ind_n; %indexes of inner , boundary edges and inner nodes
        %   but these should be ...
        DDy, DDz
    end
    
    methods
        %*******************************************************************
        function obj = TGrid2D(Dy,Dz,Nza)
            if nargin ==3
                %   class constructor ... simple
                obj.Ny = length(Dy);
                obj.Nz = length(Dz);
                obj.Nza = Nza;
                obj.NzEarth = obj.Nz-obj.Nza;
                obj.Dy = Dy;
                obj.Dz = Dz;
                obj.dualLengths;
                obj.setIndices;
            end
        end
        %*******************************************************************
        function dualLengths(obj)
            %   default setting of dual lengths
            %            -- will reset DDZ at edges for mult-resolution case
            obj.DDy = ([obj.Dy; 0]+[0; obj.Dy])/2;
            obj.DDz = ([obj.Dz; 0]+[0; obj.Dz])/2;
        end
        %*******************************************************************
        function result = NumberEdges(obj)
            result=obj.Ny*(obj.Nz+1) + (obj.Ny+1)*obj.Nz;
        end
        %*******************************************************************
        function result = NumberCells(obj)
            result = obj.Ny*obj.Nz;
        end
        %*******************************************************************
        function result = NumberNodes(obj)
            result = (obj.Ny+1)*(obj.Nz+1);
        end
        %*******************************************************************
        function [ny,nz] = setLimits(obj,type)
            %   set arraay limits for components of types
            %                     cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
            %   should probably add range checks
            switch lower(type)
                case 'cell'
                    ny = obj.Ny;
                    nz = obj.Nz;
                case 'node'
                    ny = obj.Ny+1;
                    nz = obj.Nz+1;
                case 'yedge'
                    ny = obj.Ny;
                    nz = obj.Nz+1;
                case 'zedge'
                    ny = obj.Ny+1;
                    nz = obj.Nz;
            end
        end
        %*******************************************************************
        function [J,K] = gridIndex(obj,index,type)
            %  cells, nodes, edges, faces are enumerated as column vector
            %   elements in the standard way, consistent with
            %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
            %  dimension X(Nx,Ny,Nz).   Given array of integer indices of
            %   vector elements this function, computes the
            %  corresponding i,j,k for an array of given type
            %   allowabel types : cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
            
            [ny,~] = obj.setLimits(type);
            J  = mod(index,ny);
            J(J==0) = ny;
            K = ceil(index/ny);
        end
        %*******************************************************************
        function [index] = vectorIndex(obj,J,K,type)
            %  cells, nodes, edges, faces are enumerated as column vector
            %   elements in the standard way, consistent with
            %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
            %  dimension X(Nx,Ny,Nz).   Given array of integer indices of
            %   vector elements this function computes the indicies in the
            %   vector corresponding to a lis of indicies I,J,K for an array
            %   of given type
            %   allowabel types : cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
            
            [ny,~] = obj.setLimits(type);
            index = (K-1)*ny+J;
        end %vectorIndex
        %******************************************************************
        function index =  vectorIndexRanges(obj,type,rangeString)
            %   return vector indicies for field of given type (i.e., for
            %   the index in the array of all edges, for xedge, yedge or
            %   zedge) corresponding to indicies specified by string
            %   "rangeStrine" --- e.g. something like ':,1:2:end,3'
            %   for now just coding this for edges, but could easily be extended
            %   when someone has time ...
            %
            
            [ny,nz] = setLimits(obj,type);
            temp = zeros(ny,nz);
            %eval(['temp(' rangeString ') = 1']);
            eval(['temp(' rangeString ') = 1;']);
            index = find(temp);
            switch lower(type)
                case 'yedge'
                case 'zedge'
                    [ny,nz] = setLimits(obj,'yedge');
                    index = index+ny*nz;
                case 'xface'
                case 'yface'
                otherwise
                    %   for node or cell, just use index as is
            end
        end
        %%
        function setIndices(obj)
            % set here counters for Edges, Faces and Nodes of the grid
            % interior nodes, boundaries and nodes indeces are set
             
            obj.NCells = obj.Ny*obj.Nz;
            %number of nodes
            obj.NNodes = obj.NumberNodes;
            % number of edges
            nye = obj.Ny*(obj.Nz+1);
            nze =(obj.Ny+1)*obj.Nz;
            obj.NEdges = [nye;nze];
            
            % fill in edges with zeros
            y = zeros(obj.Ny,obj.Nz+1);
            z = zeros(obj.Ny+1,obj.Nz);
            % find ind of interior and boundary edges by setting
            % boundaries to 1
            y(1:obj.Ny,[1 obj.Nz+1]) = 1;
            z([1 obj.Ny+1],1:obj.Nz) = 1;
            
            %
            vec = [reshape(y,obj.NEdges(1),1);reshape(z,obj.NEdges(2),1)];
            
            % boundary edges
            obj.ind_b = find(vec);
            % inner edges
            obj.ind_i = find(vec == 0);
            % inner nodes indices
           
            j = reshape(repmat((2:obj.Ny)',1,obj.Nz-1),(obj.Ny-1)*(obj.Nz-1),1);           
%            j = reshape(repmat((2:obj.Ny),obj.Nz-1,1),(obj.Ny-1)*(obj.Nz-1),1);
            k = reshape(ones((obj.Ny-1),1)*(2:obj.Nz),(obj.Ny-1)*(obj.Nz-1),1);
            obj.ind_n = obj.vectorIndex(j,k,'node')';
        end;
 
        
        %%
        function gr = oldGrid(obj)
            %   for consistency with old grid structure, need to copy
            %   Tgrid3D object with slightly different property names
            gr = struct('Ny',obj.Ny,'Nz',obj.Nz,...
                'Nza',obj.Nza,'dx','dy',obj.Dy,'dz',obj.Dz,'NzEarth',obj.NzEarth);
        end
        %
        
        
    end     % methods
end    % classdef
