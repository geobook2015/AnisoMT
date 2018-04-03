classdef TModelOperator2D_Aniso < TModelOperator2D
% Zeqiu Guo, Gary Egbert 2015
% zeqiu_guo@hotmail.com; egbert@coas.oregonstate.edu
% CUGB, Beijing China; OSU, OR USA
    % revised to solve coupled system of three electrical components, not
    % Ex and Hx any more, following the staggered style of discritization
    % so here, Ex at corner (actually edge along the inviriant X
    % direction), Ez and Ez along edges in Y and Z directions respectively.
    % revised to define Ex at corner/node, Hx on face/cell center
    % after discussing with Gary;
    % For anisotropic case, Ex is defined on edges,Hx is defined at nodes
    % G[y,z] mapping from node to edge, G'[y,z] mapping from edge to node
    properties
        Gy, Gz, Gfn, Gef_yz, Gef_zy
        Gya, Gza, Gefa_yz, Gefa_zy
        Gnf
    end
    
    methods
        function obj = TModelOperator2D_Aniso(grid)
            obj = obj@TModelOperator2D(grid);
        end
        function setGradyz(obj)
            %   create sparse matrix for 2D gradient operator
            gr = obj.grid;
            
            % edge to face/cell in y-direction for z-edges for first
            % derivatives
            n = gr.Ny*gr.Nz; %Number of face/cell
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);   
            [J,K] = gr.gridIndex(1:n,'cell');
            J2 = J+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'zedge'); gr.vectorIndex(J2,K,'zedge')], n*2,1);
            entries = reshape([-ones(1,n); ones(1,n)],n*2,1);
            obj.Gef_yz = sparse(rowIndicies,colIndicies,entries);  
            
            % edge to face/cell in z-direction for y-edges for first
            % derivatives
            n = gr.Ny*gr.Nz; % Number of face/cell
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);   
            [J,K] = gr.gridIndex(1:n,'cell');
            K2 = K+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'yedge'); gr.vectorIndex(J,K2,'yedge')], n*2,1);
            entries = reshape([-ones(1,n); ones(1,n)],n*2,1);
            obj.Gef_zy = sparse(rowIndicies,colIndicies,entries);     
            
            
            % edge to face/cell in y-direction for z-edges for average
            n = gr.Ny*gr.Nz; %Number of face/cell
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);   
            [J,K] = gr.gridIndex(1:n,'cell');
            J2 = J+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'zedge'); gr.vectorIndex(J2,K,'zedge')], n*2,1);
            entries = reshape([1/2*ones(1,n); 1/2*ones(1,n)],n*2,1);
            obj.Gefa_yz = sparse(rowIndicies,colIndicies,entries);  
            
            % edge to face/cell in z-direction for y-edges for average
            n = gr.Ny*gr.Nz; % Number of face/cell
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);   
            [J,K] = gr.gridIndex(1:n,'cell');
            K2 = K+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'yedge'); gr.vectorIndex(J,K2,'yedge')], n*2,1);
            entries = reshape([1/2*ones(1,n); 1/2*ones(1,n)],n*2,1);
            obj.Gefa_zy = sparse(rowIndicies,colIndicies,entries);        
            
            % node to face/cell for average
            n = gr.NCells;
            i1 = 1; i2 = n;            
            rowIndicies = reshape(ones(4,1)*(i1:i2),n*4,1);   
            [ J, K ] = gr.gridIndex(1:n,'cell');  % index of interior node
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'node'); gr.vectorIndex(J+1,K,'node');...
                gr.vectorIndex(J,K+1,'node'); gr.vectorIndex(J+1,K+1,'node')], n*4,1);
            entries = reshape([1/4*ones(1,n); 1/4*ones(1,n);1/4*ones(1,n); 1/4*ones(1,n)],n*4,1);
            obj.Gnf = sparse(rowIndicies,colIndicies,entries);     
       
            
            %  Y-EDGES    node to edge for first derivative
            n = gr.NEdges(1);
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
            [J,K] = gr.gridIndex(1:n,'yedge');
            J2 = J+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'node'); gr.vectorIndex(J2,K,'node')], n*2,1);
            entriesy = reshape([-ones(1,n); ones(1,n)],n*2,1);
            obj.Gy = sparse(rowIndicies,colIndicies,entriesy);            
            
            % Z-EDGES    node to edge for first derivative
            n = gr.NEdges(2);
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
            [J,K] = gr.gridIndex(1:n,'zedge');
            K2 = K+1;
            colIndicies = reshape([ gr.vectorIndex(J,K,'node'); gr.vectorIndex(J,K2,'node')], n*2,1);
            entriesz = reshape([-ones(1,n); ones(1,n)],n*2,1);
            obj.Gz = sparse(rowIndicies,colIndicies,entriesz);
            
            %  Y-EDGES    node to edge  for average      
            n = gr.NEdges(1);
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
            [J,K] = gr.gridIndex(1:n,'yedge');
            J2 = J+1;
            colIndicies = reshape( ...
                [gr.vectorIndex(J,K,'node'); gr.vectorIndex(J2,K,'node')], n*2,1);
            entriesy = reshape([1/2*ones(1,n); 1/2*ones(1,n)],n*2,1);
            obj.Gya = sparse(rowIndicies,colIndicies,entriesy);            
            
            % Z-EDGES    node to edge for average
            n = gr.NEdges(2);
            i1 = 1; i2 = n;
            rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
            [J,K] = gr.gridIndex(1:n,'zedge');
            K2 = K+1;
            colIndicies = reshape([ gr.vectorIndex(J,K,'node'); gr.vectorIndex(J,K2,'node')], n*2,1);
            entriesz = reshape([1/2*ones(1,n); 1/2*ones(1,n)],n*2,1);
            obj.Gza = sparse(rowIndicies,colIndicies,entriesz);            
        end
        
    end % methods
end % classdef
