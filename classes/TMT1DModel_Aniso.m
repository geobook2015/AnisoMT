classdef TMT1DModel_Aniso < TMT1DModel
    properties
        Axx,Ayy,Axy,Ayx
        cxx,cyy,czz,cxy,cxz,cyz
        F
    end
    methods
        function obj = TMT1DModel_Aniso(grid,m)
            obj = obj@TMT1DModel(grid,m);
        end
        function setsigma(obj, m)
            switch ndims(m.cxx)
                case 2
                    [~,n2] = size(m.cxx);
                    if n2 > 1
                        obj.cxx = squeeze(m.cxx(1,:)).';
                        obj.cyy = squeeze(m.cyy(1,:)).';
                        obj.czz = squeeze(m.czz(1,:)).';
                        obj.cxy = squeeze(m.cxy(1,:)).';
                        obj.cxz = squeeze(m.cxz(1,:)).';
                        obj.cyz = squeeze(m.cyz(1,:)).';
                    else
                        obj.cxx = m.cxx;
                        obj.cyy = m.cyy;
                        obj.czz = m.czz;
                        obj.cxy = m.cxy;
                        obj.cxz = m.cxz;
                        obj.cyz = m.cyz;
                    end
                case 3
                    obj.cxx = squeeze(m.cxx(1,1,:));
                    obj.cyy = squeeze(m.cyy(1,1,:));
                    obj.czz = squeeze(m.czz(1,1,:));
                    obj.cxy = squeeze(m.cxy(1,1,:));
                    obj.cxz = squeeze(m.cxz(1,1,:));
                    obj.cyz = squeeze(m.cyz(1,1,:));
                otherwise
                    error('input model parameter must be of dimension 2 or 3')
            end
            switch upper(m.paramType)
                % currently just 'Linear' assumed for anisotropy
                case 'LOGE'
                    error('case LOGE not coded for anisotropic model parameters')
                case 'LINEAR'
                    obj.cxx=[ones(obj.Nza,1)*m.AirCond; obj.cxx];
                    obj.cyy=[ones(obj.Nza,1)*m.AirCond; obj.cyy];
                    obj.czz=[ones(obj.Nza,1)*m.AirCond; obj.czz];
                    obj.cxz=[zeros(obj.Nza,1); obj.cxz];
                    obj.cxy=[zeros(obj.Nza,1); obj.cxy];
                    obj.cyz=[zeros(obj.Nza,1); obj.cyz];
                    obj.Axx=obj.cxx-obj.cxz.^2./obj.czz;
                    obj.Ayy=obj.cyy-obj.cyz.^2./obj.czz;
                    obj.Axy=obj.cxy-obj.cxz.*obj.cyz./obj.czz;
                    obj.Ayx=obj.cxy-obj.cxz.*obj.cyz./obj.czz;
            end
        end
        function setEquations(obj,omega)
            if nargin == 2
                obj.omega = omega;
            end
%             obj.Axx(obj.Nza) = obj.Axx(obj.Nza+1);
%             obj.Ayy(obj.Nza) = obj.Ayy(obj.Nza+1);
%             obj.Axy(obj.Nza) = obj.Axy(obj.Nza+1);
%             obj.Ayx(obj.Nza) = obj.Ayx(obj.Nza+1);
            %   set curl-curl operator for 1D problem; centered differences
            sigAvg_xx=obj.Axx.*obj.Dz;
            sigAvg_xx=(sigAvg_xx(1:end-1)+sigAvg_xx(2:end))./(obj.Dz(1:end-1)+obj.Dz(2:end));
            sigAvg_yy=obj.Ayy.*obj.Dz;
            sigAvg_yy=(sigAvg_yy(1:end-1)+sigAvg_yy(2:end))./(obj.Dz(1:end-1)+obj.Dz(2:end));
            sigAvg_xy=obj.Axy.*obj.Dz;
            sigAvg_xy=(sigAvg_xy(1:end-1)+sigAvg_xy(2:end))./(obj.Dz(1:end-1)+obj.Dz(2:end));
            sigAvg_yx=obj.Ayx.*obj.Dz;
            sigAvg_yx=(sigAvg_yx(1:end-1)+sigAvg_yx(2:end))./(obj.Dz(1:end-1)+obj.Dz(2:end));
            
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
            
            % revised to use impedance bottom BCs
            da12=sqrt((obj.Axx(end)-obj.Ayy(end))^2+4*obj.Axy(end)*obj.Ayx(end));
            A1=(obj.Axx(end)+obj.Ayy(end)+da12)/2;
            A2=(obj.Axx(end)+obj.Ayy(end)-da12)/2;
            if da12 >0
                blt = -0.5*acos((obj.Axx(end)-obj.Ayy(end))/da12);
            else
                blt = 0;
            end
            ep1 = sqrt(-1i*obj.omega*obj.mu/A1);
            ep2 = sqrt(-1i*obj.omega*obj.mu/A2);
            Zb = -0.5*[ -(ep1-ep2)*sin(2*blt) (ep1+ep2)+(ep1-ep2)*cos(2*blt);
                -(ep1+ep2)+(ep1-ep2)*cos(2*blt) (ep1-ep2)*sin(2*blt) ]; % e(+iwt)
            I = [ 0 0 0 -1 0 0 0 0; 0 0 0 0 0 0 0 -1];
            d1 = obj.Dz(end-2);
            d2 = obj.Dz(end-1);
            d3 = obj.Dz(end);
            x0 = obj.Dz(end-2)/2;  x1 = obj.Dz(end-2)+obj.Dz(end-1)/2;
            x2 = sum(obj.Dz(end-2:end-1))+obj.Dz(end)/2; x3 = sum(obj.Dz(end-2:end));
            l0 = (x3-x1)*(x3-x2)/(x0-x1)/(x0-x2);
            l1 = (x3-x0)*(x3-x2)/(x1-x0)/(x1-x2);
            l2 = (x3-x0)*(x3-x1)/(x2-x0)/(x2-x1);

            D = [ 0 0 0 0 l0/d1/(1i*obj.omega*obj.mu) (-l0/d1+l1/d2)/(1i*obj.omega*obj.mu) (-l1/d2+l2/d3)/(1i*obj.omega*obj.mu) -l2/d3/(1i*obj.omega*obj.mu) ;
                 -l0/d1/(1i*obj.omega*obj.mu) (l0/d1-l1/d2)/(1i*obj.omega*obj.mu) (l1/d2-l2/d3)/(1i*obj.omega*obj.mu) l2/d3/(1i*obj.omega*obj.mu) 0 0 0 0 ];
            IZbD = -I + Zb*D;        %  [ Ex Hy ] = -Z*[ Hx Hy]                       
            
            %   the diagonal part: sigma, omega, BCs: there will be one  equation for
            %   each node, including top and bottom boundary nodes
%             ioms_xx = [ 1; obj.isign*1i*obj.omega*obj.mu*sigAvg_xx; 1];
%             ioms_yy = [ 1; obj.isign*1i*obj.omega*obj.mu*sigAvg_yy; 1];
            ioms_xx = [ 1; obj.isign*1i*obj.omega*obj.mu*sigAvg_xx; 0];
            ioms_yy = [ 1; obj.isign*1i*obj.omega*obj.mu*sigAvg_yy; 0];            
            ioms_xy = [ 0; obj.isign*1i*obj.omega*obj.mu*sigAvg_xy; 0];
            ioms_yx = [ 0; obj.isign*1i*obj.omega*obj.mu*sigAvg_yx; 0];
            
            obj.A = [obj.A, zeros(size(obj.A));zeros(size(obj.A)),obj.A];
            obj.A(obj.Nz+1,[ obj.Nz-2:obj.Nz+1 end-3:end ]) = IZbD(1,:);
            obj.A(end,[ obj.Nz-2:obj.Nz+1 end-3:end ]) = IZbD(2,:);            
            obj.A = obj.A + ...
                [spdiags(ioms_xx,0,obj.Nz+1,obj.Nz+1), spdiags(ioms_xy,0,obj.Nz+1,obj.Nz+1);
                spdiags(ioms_yx,0,obj.Nz+1,obj.Nz+1), spdiags(ioms_yy,0,obj.Nz+1,obj.Nz+1)];
        end
        
        function Solve(obj,pol)
            % Here, solution of 1D anisotropy is Ex and Ey, while
            % 2D anisotropy needs Ex and Hx as BCs, so need to
            % transform Ey to Hx with D(Ey)/D(z) = -iwu*Hx in
            % 1D condition.
            switch pol
                case 'X'
                    b = [1; zeros(obj.Nz*2+1,1)];
                case 'Y'
                    b = [zeros(obj.Nz+1,1); 1; zeros(obj.Nz,1)];
            end
            Ee = obj.A\b;
            Exd = (Ee(1:obj.Nz)+Ee(2:obj.Nz+1))/2;
            Eyd = (Ee(obj.Nz+2:obj.Nz*2+1)+Ee(obj.Nz+3:obj.Nz*2+2))/2;
            
            obj.F.Ex = Ee(1:obj.Nz+1);
            obj.F.Ey = Ee(obj.Nz+2:end);
            obj.F.Ez = -(obj.cxz.*Exd+obj.cyz.*Eyd)./obj.czz;
            obj.F.Hx = (obj.F.Ey(2:end)-obj.F.Ey(1:end-1))./(-obj.isign*obj.Dz*1i*obj.mu*obj.omega)*(-1);
            obj.F.Hy = (obj.F.Ex(2:end)-obj.F.Ex(1:end-1))./(-obj.isign*obj.Dz*1i*obj.mu*obj.omega);
            % Hz = 0 for 1D case
        end
        
        function Plot(obj,depthRange)
            %  simple plot vs. depth
            figure('Position',[100,100,800,800])
            z = cumsum([0 ; obj.Dz]);
            z = z-z(obj.Nza+1);
            if nargin == 1
                i1 = 1;
                i2 = length(z);
            else
                i1 = find(min(abs(z-depthRange(1)))==abs((z-depthRange(1))));
                i2 = find(min(abs(z-depthRange(2)))==abs((z-depthRange(2))));
            end
            subplot(121)
            plot(real(obj.F.Ex(i1:i2)),z(i1:i2) ,'b-','linewidth',2);
            axis('ij')
            hold on
            plot(imag(obj.F.Ex(i1:i2)),z(i1:i2), 'r--','linewidth',2);
            set(gca,'FontSize',16,'FontWeight','demi')
            legend('Real(Ex)','Image(Ex)')
            subplot(122)
            plot(real(obj.F.Ey(i1:i2)),z(i1:i2) ,'b-','linewidth',2);
            axis('ij')
            hold on
            plot(imag(obj.F.Ey(i1:i2)),z(i1:i2), 'r--','linewidth',2);
            set(gca,'FontSize',16,'FontWeight','demi')
            legend('Real(Ey)','Image(Ey)')
        end
    end
end
