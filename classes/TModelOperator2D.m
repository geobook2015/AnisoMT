classdef TModelOperator2D < handle
    % Egbert
    %
    % this does full model operator, topology and metric elements for 2D
    % grids
    
    properties
        grid %   
        G    %   topology of gradient: sparse matrix mapping nodes to edges
        Gset
        ie               %   interior edge indices
        be               %   boundary edge indicies
        in               %   interior node indicies
        bn               %   boundary node indicies

        % metric element array
        EdgeLength             %   array defining edge length
        DualEdgeLength            %   array defining dual edge lengths
        NodeArea           %   array of node areas
        CellArea           %   array of cell areas
    end
    
    methods
        
        function obj = TModelOperator2D(grid)
            %  class constructor
            if nargin == 1
                obj.grid = grid;
            end
        end
        %******************************************************************
        function setIndicies(obj)
            x = TVector2D(obj.grid);
            [obj.ie,obj.be] = x.int_bdry_indices;
            x =TScalar2D(obj.grid,'node');
            [obj.in,obj.bn] = x.int_bdry_indices;
        end
        %**************************************************************************
        %   Better to not do set up automatically -- need more control about how
        %    these objects are used
        function setMetricElements(obj)
            obj.setCellArea;
            obj.setEdgeLength;
            obj.setNodeArea;
            obj.setDualEdgeLength;
        end
        %******************************************************************
        function setCellArea(obj)
            %   create 1-D array of cell volumes, ordered as indC
            [DY,DZ] = ndgrid(obj.grid.Dy,obj.grid.Dz);
            obj.CellArea = DY.*DZ;
        end
        %
        %******************************************************************
        function setNodeArea(obj)   
            %   create 1-D array of node volumens, ordered as indN
            %    Let's just do this with a short-cut here, so that it
            %    is sensible for the special case where all cells in each
            %    vertical layer have the same coarseness
            %    SET dual edge lengths first
            [DY,DZ] = ndgrid(obj.grid.DDy,obj.grid.DDz);
            obj.NodeArea = DY.*DZ;
        end
        %
        %******************************************************************
        function setEdgeLength(obj)
            %   creates 1-D array of edge lengths
            nye = obj.grid.Ny*(obj.grid.Nz+1);
            nze = (obj.grid.Ny+1)*obj.grid.Nz;
            obj.EdgeLength = zeros(nye+nze,1);
            %   y-edges
            i1 = 1; i2 = obj.grid.NEdges(1);
            obj.EdgeLength(i1:i2) = reshape(obj.grid.Dy*ones(1,obj.grid.Nz+1),nye,1);
            % z-edges
            i1 = i2+1; i2 = i2+nze;
            obj.EdgeLength(i1:i2) = reshape(ones(obj.grid.Ny+1,1)*obj.grid.Dz',nze,1);
        end
        %******************************************************************
        function setDualEdgeLength(obj)
            %   creates 1-D array of dual edge lengths, one for each face
            nye = obj.grid.Ny*(obj.grid.Nz+1);
            nze = (obj.grid.Ny+1)*obj.grid.Nz;
            obj.DualEdgeLength = zeros(nye+nze,1);
            % dual y-edges  (one for each z-edge
            i1 = 1; i2 = nye;
            % dual z-edges  (one for each y-edge
            obj.DualEdgeLength(i1:i2) = reshape(ones(obj.grid.Ny,1)*obj.grid.DDz',nye,1);
             % dual y-edges  (one for each z-edge
            i1 = i2+1; i2 = i2+nze;
            obj.DualEdgeLength(i1:i2) = reshape(obj.grid.DDy*ones(1,obj.grid.Nz),nze,1);
        end
        %******************************************************************
        function setGrad(obj)
            %   create sparse matrix for 2D gradient operator
            %  Y-EDGES
            gr = obj.grid;
            n = gr.NEdges(1);
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
            %   columns
            [J,K] = gr.gridIndex(1:n,'yedge');
            J2 = J+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'node'); gr.vectorIndex(J2,K,'node')], n*2,1);
            %   matrix entries
            entries = reshape([-ones(1,n); ones(1,n)],n*2,1);
            
            % Z-EDGES
            n = gr.NEdges(2);
            i1 = i2+1; i2 = i2+n;
            rowIndicies = [rowIndicies; reshape(ones(2,1)*(i1:i2),n*2,1)];
            %   columns
            [J,K] = gr.gridIndex(1:n,'zedge');
            K2 = K+1;
            colIndicies = [ colIndicies; ...
                reshape([ gr.vectorIndex(J,K,'node'); gr.vectorIndex(J,K2,'node')], n*2,1)];
            %   matrix entries
            entries = [entries; reshape([-ones(1,n); ones(1,n)],n*2,1)];
            %  define sparse matrix
            obj.G = sparse(rowIndicies,colIndicies,entries);
            obj.Gset = true; 
        end
        
    end % methods
end % classdef
