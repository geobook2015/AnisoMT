classdef TSource2DMT_1DBC < TSource
    %   specific instance of source class for 3D MT with BC from 1D model
    %     really this is for case SG only as coded ... not sure this can
    %     be easily generalized
    properties
        E1D   %   1D modeling object
        sigmaLeft   %  1D profile on left side
        sigmaRight  %  1D profile on right side
    end
    
    methods
        function obj = TSource2DMT_1DBC(grid,m)
            if nargin ==2
                obj.grid = grid;
                %   don't set conductivity in this case
                obj.E1D = TMT1DModel_ImpBCBottom(grid,m);
%                obj.E1D = TMT1DModel(grid,m); 
                switch upper(m.paramType)
                    case 'LOGE'
                        obj.sigmaLeft = exp([ones(obj.grid.Nza,1)*m.AirCond; m.v(1,:)']);
                        obj.sigmaRight = exp([ones(obj.grid.Nza,1)*m.AirCond; m.v(end,:)']);
                    case 'LINEAR'
                        obj.sigmaLeft = [ones(obj.grid.Nza,1)*exp(m.AirCond); m.v(1,:)'];
                        obj.sigmaRight = [ones(obj.grid.Nza,1)*exp(m.AirCond); m.v(end,:)'];
                end
            end
        end
        
        %*******************************************************************
        function SetSourceParams(obj,omega,Polarization)
            %  complete setup of 1D equations, solve
            obj.SourceParams = struct('omega',omega,'Polarization',Polarization);
            obj.nonZeroSource = false;
         end
        %*******************************************************************
        function setRHS(obj)
            %    this sets rhs in B
            obj.B = TScalar2D(obj.grid);
            %   first finish setting up 1D equations
            %   LHS
            obj.E1D.sigma = obj.sigmaLeft;
            obj.E1D.setEquations(obj.SourceParams.omega,obj.SourceParams.Polarization);
            obj.E1D.Solve;
            eLeft = obj.E1D.e;
            % RHS
            obj.E1D.sigma = obj.sigmaRight;
            obj.E1D.setEquations(obj.SourceParams.omega,obj.SourceParams.Polarization);
            obj.E1D.Solve;
            eRight = obj.E1D.e;
            
            %  top
            obj.B.v(:,1) = 1;
            % left
            obj.B.v(1,:) = eLeft;
            %  right
            obj.B.v(end,:) = eRight;
            %   bottom ... linearly interpolate
            xi  = [0; cumsum(obj.grid.Dy)];
            X = [0;xi(end)];
            Y = [eLeft(end);eRight(end)];
            obj.B.v(:,end) = interp1(X,Y,xi);
        end
        
        %*******************************************************************
        %    finally a routine to initialize solution vector
        function initSolnVector(obj)
        end
    end
end
