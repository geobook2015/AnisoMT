classdef TMT1DModel < handle
    % 1D forward class -- basic version with specificied values at top and
    % bottom (normally 1 for top, 0 for bottom -- grid has to go deep
    % enough)
    properties
        Dz
        DDz
        Nz
        Nza
        sigma
        A
        e
        omega
        isign = -1;   % +1 for exp(+i omega t)
        mu		 = pi*4*1E-07; %12.566E-7 %  magnetic permeability so far a scalar values
    end
    methods
        function obj = TMT1DModel(grid,m)
            %   for now this sets up the MT1DModel object using grid and
            %   model parameter objects used for 3D forward model  (or 2D
            %   ...  does not really depend on grid being a class; could
            %   just be a structure
            if nargin >=1
                obj.Nz = length(grid.Dz);
                obj.Dz = grid.Dz;
                obj.Nza = grid.Nza;
                obj.DDz = (obj.Dz(1:end-1)+obj.Dz(2:end))/2;
            end
            if nargin == 2
                obj.setsigma(m);
            end
        end
        %%
        function setsigma(obj,m)
            %    this is old simple routine for setting up 1D
            %    conductivity; modified to work for 2D and 3D (and 1D)
            %    this just takes the first 1D verticle profile from the
            %    array ... for other options add code.
            switch ndims(m.v)
                case 2
                    [~,n2] = size(m.v);
                    if (n2>1)
                         temp = squeeze(m.v(1,:)).';
                    else
                        temp = m.v;
                    end
                case 3
                    temp = squeeze(m.v(1,1,:));
                otherwise
                    error('input model parameter must be of dimension 1, 2 or 3')
            end
            switch upper(m.paramType)
                case 'LOGE'
                    obj.sigma = exp([ones(obj.Nza,1)*m.AirCond; temp]);
                case 'LINEAR'
                    obj.sigma = [ones(obj.Nza,1)*exp(m.AirCond); temp];
            end
        end;
        %%
        function setEquations(obj,omega,job)
            %   I am changing "pol" to "job" ... really we are talking
            %   about solving 1D equations in terms of E or H, not TE or TM
            %   But will continue to support the old way of calling.
            if nargin == 2
                job = 'E';
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
                    %   the diagonal part: sigma, omega, BCs: there will be one  equation for
                    %   each node, including top and bottom boundary nodes
%                     ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*sigmaAvg; 1];                    
                    ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*sigmaAvg; 1];
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
                    %   the diagonal part: sigma, omega, BCs: there will be one  equation for
                    %   each node, including top and bottom boundary nodes
                    ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*ones(obj.Nz-1,1); 1];
%                     ioms = [ 1; obj.isign*1i*obj.omega*obj.mu*ones(obj.Nz-1,1); 1];                    
                    obj.A = obj.A + spdiags(ioms,0,obj.Nz+1,obj.Nz+1);   
            end
        end
        %%
        function Solve(obj)
           
            %  solve directly --simple preconditioner might prevent warinings
            %    about poor-conditioning ...
            d = diag(obj.A);
            d = spdiags(1./d,0,obj.Nz+1,obj.Nz+1);
            P = d*obj.A;
             %   Boundary data --- to start let's assume zero at the bottom
            b = d*[1; zeros(obj.Nz,1)];
            obj.e =P\b;
        end
        %%
        function Plot(obj,depthRange)   
            %  simple plot vs. depth
            figure('Position',[100,100,400,800])
            z = cumsum([0 ; obj.Dz]);
            z = z-z(obj.Nza+1);
            if nargin == 1
                i1 = 1;
                i2 = length(z);
            else
                i1 = find(min(abs(z-depthRange(1)))==abs((z-depthRange(1))));
                i2 = find(min(abs(z-depthRange(2)))==abs((z-depthRange(2))));
            end
            plot(real(obj.e(i1:i2)),z(i1:i2) ,'b-','linewidth',2);
            axis('ij')
            hold on
            plot(imag(obj.e(i1:i2)),z(i1:i2), 'r--','linewidth',2);
            set(gca,'FontSize',16,'FontWeight','demi')
            
        end
    end
end
