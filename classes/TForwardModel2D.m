classdef TForwardModel2D < handle
    
    %   Guess what?  Yet again I don't think this should be a subclass of
    %   anything.   This is supposed to be the "Solver Driver" object --
    %   a handle object that is created, contains matrices that define
    %   equations, etc.   This will have setup routines, routines for changing
    %   model parameters and frequency, routines for using sources and BCs to
    %   set up RHS terms for equations, and drivers to control solution.   It
    %   will manage efficient solutions for multiple RHS (without
    %   reinitiallizing things that don't change).   I envision this as a
    %   polymorphic class -- if we wanted to set up different equations and
    %   solve things in different ways some methods could be overloaded in
    %   subclasses.   This is only thought out at a general level, details
    %   need more work.   (BUT: the principals outlined above have been thought
    %   about and undrlie strucure of ModEM)
    
    properties
        grid   %   Yep, this underlies everything, but this is NOT going to be
        %     an extension of a grid object.
        modOp  %   this points to the modelOperator object, which is
        %    used to construct the equations (many ways to do this ...
        %     and there are variants on the modelOperator object (SG,
        %     MR, Cartesian, spherical)   Once grid is defined, these
        %     can be set.
        m      %   The model parameter object.   This DOES NOT necessarily
        %     depend on the grid  (in fact, in our MR implementation
        %     this depends on the underlying fine grid, but not the
        %     actual MR grid details!)
        omega  %   Frequency
        A      % coefficient matrix for solver -- only plan to support direct
               %    solution for 2D
        % Afactored = false
        %L      %  lower triangular factor
        %U      %   upper triangular factor
        mode = 'TE'  % options are TE/TM
    end
    properties %(Access = private)
        mu = pi*4E-7   % just single value for now
        epsilon = 1
        isign = -1   %   this should be used to define sign dependence of time
        %    dependence: exp(ISIGN * 1i* omega * t)
        %   Usually in ModEM we use ISIGN = -1, consitent with
        %   normal  data processing conventions (i.e., standard
        %   definition of most fft implementations)    
        PIm        %    Model parameter mapped to primary or secondary grid--e.g.,
        %      conductivity on nodes for TE, resistivity on edges for TM.  Not sure
        %      this is really useful to keep around or not.
    end
    methods
        function obj = TForwardModel2D(grid,mode,m,omega)
            if nargin >= 1
                obj.grid = grid;
                obj.modOp = TModelOperator2D(obj.grid);
            end
            
            if nargin >=2
                obj.mode = mode;
            end      %  lower triangular factor
        
            
            if nargin >= 3
               
                %     PIm is in this case conductivity mapped to edges ... but
                %     could be something different, such as resistivity mapped
                %     to faces.  PI(m) in the notation of Egbert and Kelbert
                %     (2012)
                obj.PIm = m.PDEmapping(mode).GetVector;
            end
            
            if nargin ==4
                obj.omega = omega;
            end  
        end
        %******************************************************************
        function setEquations(obj)
            %   this routine initializes matrices -- but only those parts
            %   that do not depend on omega or sigma;  the additional
            %   initialilzation involving these parameters is done
            %   separately, to avoid repeating these parts that remain the
            %   same for all conductivities, all periods.  To complete
            %   setup of operators need to call setFreqCond
            
            %   for 2D better to do all operator setup in setFreqCond
            obj.modOp.setGrad;
            obj.modOp.setMetricElements;
            obj.modOp.setIndicies;
           
        end
        %*******************************************************************

        function setFreqCond(obj,omega,m)
            %     This routine sets additional matrices which depend on
            %     conductivity, using matrices that are more generic
            %
            %     PIm is in this case conductivity mapped to edges ... but
            %     could be something different, such as resistivity mapped
            %     to faces.  PI(m) in the notation of Egbert and Kelbert
            %     (2012)
            
            %    for 2D version: factor matrix after forming
            if nargin == 3
                obj.PIm = m.PDEmapping(obj.mode).GetVector;
            end
            if nargin >= 2
                obj.omega  = omega;
            end
            switch obj.mode
                case 'TE'
                    %  for TE mode case set up del^2 part of operator
                    ne = length(obj.modOp.EdgeLength);
                    nn = size(obj.modOp.NodeArea);
                    nn = nn(1)*nn(2);
                    %   this is defined for all nodes ... need to reduce rows to
                    %   include only interior nodes
                    %   Here,not integral(Hdl)/s=j .or. integral(Hdl)=j*s
                    %   just discretization of delta^2(Ex)+iwu*sigma*Ex=0
                    %   reshape(obj.modOp.NodeArea,nn,1), obj.modOp.DualEdgeLength
                    %   are for division of distance from derivatives of
                    %   fields;
                    %   function of A: node/E-to-edge/H with matrix G,
                    %   edge/H-to-node/E with matrix G'
                    obj.A =  spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G'...
                        *spdiags(obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*obj.modOp.G;
                    Vioms = obj.isign * 1i*obj.omega*obj.mu*obj.PIm;
                    n = length(Vioms);
                    obj.A = obj.A + spdiags(Vioms,0,n,n);
                case 'TM'
                    ne = length(obj.modOp.EdgeLength);
                    nn = size(obj.modOp.NodeArea);
                    nn = nn(1)*nn(2);
                    % First edge-to-node with G', then node-to-edge with G,
                    % G'*G, nEdge*nEdge
%                     obj.A =  spdiags(obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*...
%                         obj.modOp.G*spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G';
%                     Vioms = obj.isign * 1i*obj.omega*obj.mu*obj.PIm; 
                    obj.A =  spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G'...
                        *spdiags(obj.PIm.*obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*obj.modOp.G;
                    % keep coefficients of boundary nodes as zeros,
                    % although makes no big difference, cause i*w*u is too
                    % small, but this term can't be igored for other nodes
                    cof = ones(obj.grid.Ny+1,obj.grid.Nz+1); %ones to set BCs,obj.PIm is zeros at boundary nodes
                    cof(2:end-1,2:end-1) = obj.isign * 1i*obj.omega*obj.mu;
                    Vioms = reshape(cof,numel(cof),1);
%                     Vioms = ones(size(obj.A,1),1)*obj.isign * 1i*obj.omega*obj.mu;
                    n = length(Vioms);
                    obj.A = obj.A + spdiags(Vioms,0,n,n);
            end 
            %    factor matrix
            %[obj.L,obj.U] = lu(obj.A(obj.modOp.in,obj.modOp.in));
            %obj.Afactored = true;
        end
        %*******************************************************************
        function [e,E] = Solve(obj,src)
            %  solve equations for specified source
            %    Source is an input source object defining all source terms
            %      including BC.   It is an instance of the abstract
            %      TSource class
            %   The output is a column vector, but of full grid dimension,
            %   including boundary data
            src.setRHS;   % this completes definition of source terms
            ii = obj.modOp.in;
            ib = obj.modOp.bn;
            e = src.B.GetVector;

            %   need to multiply source by volume elements
            %rhs = e(ii)-obj.A(ii,ib)*e(ib);
            rhs = -obj.A(ii,ib)*e(ib);
            %   in some applications should probably do LU decomp and save
            %   factors, to accelerate subsequent solution.
            %if ~obj.Afactored
            %    [obj.L,obj.U] = LU(obj.A(ii,ii));
            %    obj.Afactored = true;
            %end
            %opts.UT = true;
            %temp = linsolve(obj.U,rhs,opts);
            %opts.UT = false;opts.LT = true;
            %e(ii) = linsolve(obj.L,temp,opts);
            %save Ae.mat obj e rhs ii
            e(ii) = obj.A(ii,ii)\rhs;
            inorm = norm(obj.A(ii,ii)*e(ii)- rhs)/norm(rhs);
            if nargout ==2
                E = TScalar2D(obj.grid).SetVector(e);
            end
        end
        %******************************************************************
        function [Hy,Hz] = BfromE(obj,E)
            %  computes magnetics (on dual edes) given electric on nodes, interpolate
            %  result to interior nodes:   
            % NOTE realy this computes B, not   H!!  method name changed to
            % reflect this fact -- variable names not changed
            e = E.GetVector;
            ne = length(obj.modOp.EdgeLength);
            temp = 1./(-obj.isign*1i*obj.omega*obj.modOp.EdgeLength);
            h = spdiags(temp,0,ne,ne)*obj.modOp.G*e;
            Hyz = TVector2D(obj.grid).SetVector(h);
            mdz1 = repmat(obj.grid.Dz(1:end-1,:),1,obj.grid.Ny+1)';
            mdz2 = repmat(obj.grid.Dz(2:end,:),1,obj.grid.Ny+1)';
            mdy1 = repmat(obj.grid.Dy(1:end-1,:)',obj.grid.Nz+1,1)';
            mdy2 = repmat(obj.grid.Dy(2:end,:)',obj.grid.Nz+1,1)';
            Hy = TScalar2D(obj.grid);            
            Hy.v(:,2:end-1) = Hyz.z(:,1:end-1).*mdz1+Hyz.z(:,2:end).*mdz2; 
            Hy.v(:,2:end-1) = Hy.v(:,2:end-1)./(mdz1+mdz2);
            Hz = TScalar2D(obj.grid);
            Hz.v(2:end-1,:) = Hyz.y(1:end-1,:).*mdy1+Hyz.y(2:end,:).*mdy2;
            Hz.v(2:end-1,:) = Hz.v(2:end-1,:)./(mdy1+mdy2);
        end
        function [Ey,Ez] = EfromH(obj,H)
            %  computes electric (on dual edes) given magnetic on nodes, interpolate
            %  result to interior nodes
            h = H.GetVector;
            ne = length(obj.modOp.EdgeLength);
            temp = obj.PIm./obj.modOp.EdgeLength;
            e = spdiags(temp,0,ne,ne)*obj.modOp.G*h;
            Eyz = TVector2D(obj.grid).SetVector(e);
            mdz1 = repmat(obj.grid.Dz(1:end-1,:),1,obj.grid.Ny+1)';
            mdz2 = repmat(obj.grid.Dz(2:end,:),1,obj.grid.Ny+1)';
            mdy1 = repmat(obj.grid.Dy(1:end-1,:)',obj.grid.Nz+1,1)';
            mdy2 = repmat(obj.grid.Dy(2:end,:)',obj.grid.Nz+1,1)';
            Ey = TScalar2D(obj.grid);            
            Ey.v(:,2:end-1) = Eyz.z(:,1:end-1).*mdz1+Eyz.z(:,2:end).*mdz2; 
            Ey.v(:,2:end-1) = Ey.v(:,2:end-1)./(mdz1+mdz2);
            Ez = TScalar2D(obj.grid);
            Ez.v(2:end-1,:) = Eyz.y(1:end-1,:).*mdy1+Eyz.y(2:end,:).*mdy2;
            Ez.v(2:end-1,:) = Ez.v(2:end-1,:)./(mdy1+mdy2);
        end        
    end
end
