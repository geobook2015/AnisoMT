classdef TMT1DModel_ImpBCBottom < TMT1DModel 
    % subclass of 1D modeling class, modified for impedance BC at bottom
    properties
    end
    methods
        function obj = TMT1DModel_ImpBCBottom(grid,m)
            %   for now this sets up the MT1DModel object using grid and
            %   model parameter objects used for 3D forward model  (or 2D
            %   ...  does not really depend on grid being a class; could
            %   just be a structure       
            obj = obj@TMT1DModel(grid,m);
        end
  
        %%
        function setEquations(obj,omega,job)
            %  this just modifes the bottom BC
            if nargin == 2
                job = 'E'
            end
            obj.omega = omega;
            switch job
                case {'TE','E'}
                    %   set curl-curl operator for 1D problem; centered differences
                    sigmaAvg = obj.sigma.*obj.Dz;
                    sigmaAvg = (sigmaAvg(1:end-1)+sigmaAvg(2:end))./(obj.Dz(1:end-1)+obj.Dz(2:end));
                    I = reshape([1:obj.Nz;1:obj.Nz],2*obj.Nz,1);
                    J =  reshape([1:obj.Nz;2:obj.Nz+1],2*obj.Nz,1);
                    S =  reshape([-ones(1,obj.Nz);ones(1,obj.Nz)],2*obj.Nz,1);
                    G = sparse(I,J,S);
                    D = G(:,2:obj.Nz)';
                    %    set up the differential operator
                    %    (1)  the curl-curl part (independent of sigma, omega)
                    obj.A = spdiags(1./obj.DDz,0,obj.Nz-1,obj.Nz-1)*D* ....
                        spdiags(1./obj.Dz,0,obj.Nz,obj.Nz)*G;
                    [I,J,S] = find(obj.A);
                    %   now add rows to top and bottom (containing 0) to allow for boundary
                    %   nodes
                    I = [1; I+1; obj.Nz+1];
                    J = [1; J; obj.Nz+1];
                    S = [0; S ; 0];
                    obj.A = sparse(I,J,S);                    
                    
                    % revised to assume new bottom 1D BC
                    zn = sum(obj.Dz(end-1:end));
                    zn_1 = obj.Dz(end-1);
                    alpha = (2*zn-zn_1)/(zn*(zn-zn_1));
                    beta = -zn/(zn_1*(zn-zn_1));
                    gama = (zn-zn_1)/(zn*zn_1);
                    c = 1/(1i*sqrt(1i*obj.omega*obj.mu*obj.sigma(end))); %c=1/(i*k)
                    obj.A(end,end) = c*alpha-1;
                    obj.A(end,end-1) = c*beta;
                    obj.A(end,end-2) = c*gama;      
                    %   the diagonal part: sigma, omega, BCs: there will be one  equation for
                    %   each node, including top and bottom boundary nodes
%                     ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*sigmaAvg; 1];                    
                    ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*sigmaAvg; 0];
                    obj.A = obj.A + spdiags(ioms,0,obj.Nz+1,obj.Nz+1);
                case {'TM','B'}
                    I = reshape([1:obj.Nz;1:obj.Nz],2*obj.Nz,1);
                    J =  reshape([1:obj.Nz;2:obj.Nz+1],2*obj.Nz,1);
                    S =  reshape([-ones(1,obj.Nz);ones(1,obj.Nz)],2*obj.Nz,1);
                    G = sparse(I,J,S);
                    D = G(:,2:obj.Nz)';
                    %    set up the differential operator
                    %    (1)  the curl-curl part (independent of sigma, omega)
                    obj.A = spdiags(1./obj.DDz,0,obj.Nz-1,obj.Nz-1)*D* ....
                        spdiags(1./obj.sigma./obj.Dz,0,obj.Nz,obj.Nz)*G;
                    [I,J,S] = find(obj.A);
                    %   now add rows to top and bottom (containing 0) to allow for boundary
                    %   nodes
                    I = [1; I+1; obj.Nz+1];
                    J = [1; J; obj.Nz+1];
                    S = [0; S ; 0];
                    obj.A = sparse(I,J,S);                    
                    
                    % revised to assume new bottom 1D BC
                    zn = sum(obj.Dz(end-1:end));
                    zn_1 = obj.Dz(end-1);
                    alpha = (2*zn-zn_1)/(zn*(zn-zn_1));
                    beta = -zn/(zn_1*(zn-zn_1));
                    gama = (zn-zn_1)/(zn*zn_1);
                    c = -sqrt(1i/(obj.omega*obj.mu*obj.sigma(end))); %c= -k/(w*u*sigma)=1/(i*k), k=(i*w*u*sigma)^1/2
                    % There will no warnings if let the depth be great enough, not sure why?
                    obj.A(end,end) = c*alpha-1;
                    obj.A(end,end-1) = c*beta;
                    obj.A(end,end-2) = c*gama;
                                        
                    %   the diagonal part: sigma, omega, BCs: there will be one  equation for
                    %   each node, including top and bottom boundary nodes
                    ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*ones(obj.Nz-1,1); 0];
%                     ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*ones(obj.Nz-1,1); 1];                    
                    obj.A = obj.A + spdiags(ioms,0,obj.Nz+1,obj.Nz+1);   
            end
        end
    end
end
