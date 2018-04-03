classdef TModelParameterCell2D_Aniso < TModelParameterCell2D
% Zeqiu Guo, Gary Egbert 2015
% zeqiu_guo@hotmail.com; egbert@coas.oregonstate.edu
% CUGB, Beijing China; OSU, OR USA
    properties
        Mmatrix = false;
        cxx,cyy,czz,cxy,cyz,cxz
        %  | Cxx Cxy Cxz |
        %  | Cyx Cyy Cyz |
        %  | Czx Czy Czz |
        D,A,B,yyD,yzD,zzD,V        
    end
    methods
        function obj = TModelParameterCell2D_Aniso(GRID,paramType)
            if nargin >=1
                obj.grid = GRID;
                obj.ParamGrid = GRID;
                if nargin == 1
                    obj.paramType='LOGE';
                else
                    obj.paramType = paramType;
                end
            end
        end
        
        function [obj] = readVec(obj,cfile)
            [dy,dz,rho,nzAir,type] = read_2dmodel_Aniso(cfile);
            
            %   add air layers, following approach in ModEM (uggg)
            %nzAir = 10;
            dzAir = dz(nzAir:-1:1).*3.^(nzAir-1:-1:0)';
            dz = [dzAir; dz];            
            Grid = TGrid2D(dy,dz,nzAir);
            obj.AirCond = 1e-10;
            obj.paramType = type;
            obj.grid = Grid;            
            switch lower(obj.paramType)
                case 'loge'
                    %                     cellCond.Cxx =cellCond.Cxx + obj.AirCond;
                    %                     cellCond.Cxx(:,:,k1:k2) = exp(obj.Cxx);
                    %                     cellCond.Cyy =cellCond.Cyy + obj.AirCond;
                    %                     cellCond.Cyy(:,:,k1:k2) = exp(obj.Cyy);
                    %                     cellCond.Czz =cellCond.Czz + obj.AirCond;
                    %                     cellCond.Czz(:,:,k1:k2) = exp(obj.Czz);
                case 'linear'
                    obj.cxx = rho(:,:,1);
                    obj.cyy = rho(:,:,2);     
                    obj.czz = rho(:,:,3);
                    obj.cxy = rho(:,:,4);
                    obj.cyz = rho(:,:,5);
                    obj.cxz = rho(:,:,6);
                otherwise
                    error(['model parameter mapping not coded for paramtype' obj.paramType ])
            end
        end
        
        function sigMap = PDEmapping(obj,style)
            k1 = obj.grid.Nza;
            ny = obj.grid.Ny;            
            obj.cxx = [ ones(ny,k1)*obj.AirCond obj.cxx ];
            obj.cyy = [ ones(ny,k1)*obj.AirCond obj.cyy ];
            obj.czz = [ ones(ny,k1)*obj.AirCond obj.czz ];
            obj.cxy = [ zeros(ny,k1) obj.cxy ];
            obj.cxz = [ zeros(ny,k1) obj.cxz ];
            obj.cyz = [ zeros(ny,k1) obj.cyz ];
%             obj.cxx(:,k1) = obj.cxx(:,k1+1);
%             obj.cyy(:,k1) = obj.cyy(:,k1+1);
%             obj.czz(:,k1) = obj.czz(:,k1+1);
%             obj.cxy(:,k1) = obj.cxy(:,k1+1);
%             obj.cxz(:,k1) = obj.cxz(:,k1+1);
%             obj.cyz(:,k1) = obj.cyz(:,k1+1);            
            switch style
                % 'ThreeE'   'EHSta'   'EHFix'
                case 'ThreeE'
                    sigMap = sigAvgThreeE(obj);
                case 'EHSta'
                    sigMap = sigAvgEH(obj);
                case 'EHFix'
                    sigMap = sigAvgEH(obj);
                otherwise
                    error('The style you specified has not been implemented');
            end
        end
        function vecOut = sigAvgThreeE(obj)
            %    average scalar object (cell type only) onto TVector object; of edge or
            %     node type;   boundary edges are left set to 0???
            
            %  Cxx, Cxy, Cxz need to be mapped on nodes
            %         Cxy, Cxz, Cyy, Cyz, Czz need to be mapped on edges
            % First assigned to faces/cells
                     
            % m--node    e--edge
            [DY,DZ] = ndgrid(obj.grid.Dy,obj.grid.Dz);
            vol = DY.*DZ;
            
            temp = obj.cxx.*vol;
            vecOut.xxm = TScalar2D(obj.grid,'node');
            vecOut.xxm.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.xxm.v(2:end-1,2:end-1) = vecOut.xxm.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
            
            temp = obj.cxy.*vol;
            vecOut.xym = TScalar2D(obj.grid,'node');
            vecOut.xym.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.xym.v(2:end-1,2:end-1) = vecOut.xym.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
            
            vecOut.yzc = TScalar2D(obj.grid,'cell');
            vecOut.yzc.v = obj.cyz;
            
            temp = obj.cxy.*vol;
            vecOut.xye = TVector2D(obj.grid);
            vecOut.xye.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.xye.y(:,2:end-1) = vecOut.xye.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            
            temp = obj.cxz.*vol;
            vecOut.xzm = TScalar2D(obj.grid,'node');
            vecOut.xzm.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.xzm.v(2:end-1,2:end-1) = vecOut.xzm.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
            
            temp = obj.cxz.*vol;
            vecOut.xze = TVector2D(obj.grid);
            vecOut.xze.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.xze.z(2:end-1,:) = vecOut.xze.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.cyy.*vol;
            vecOut.yye = TVector2D(obj.grid);
            vecOut.yye.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.yye.y(:,2:end-1) = vecOut.yye.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            
            temp = obj.czz.*vol;
            vecOut.zze = TVector2D(obj.grid);
            vecOut.zze.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.zze.z(2:end-1,:) = vecOut.zze.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.cyz.*vol;
            vecOut.yze = TVector2D(obj.grid);
            vecOut.yze.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.yze.y(:,2:end-1) = vecOut.yze.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.yze.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.yze.z(2:end-1,:) = vecOut.yze.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
        end
        function vecOut = sigAvgEH(obj)
            %    average scalar object (cell type only) onto TVector object; of edge or
            %     node type;   boundary edges are left set to 0???
            
            %  V,A,B needs to be mapped on nodes
            % A,B,yyD,yzD,and zzD need to be mapped on edges
            % First assigned to faces/cells
            obj.D = obj.cyy.*obj.czz - obj.cyz.^2;
            obj.A = (obj.cxy.*obj.cyz - obj.cxz.*obj.cyy)./obj.D;
            obj.B = (obj.cxz.*obj.cyz - obj.cxy.*obj.czz)./obj.D;
            obj.yyD = obj.cyy./obj.D;
            obj.yzD = obj.cyz./obj.D;
            obj.zzD = obj.czz./obj.D;
            obj.V = obj.cxx + obj.A.*obj.cxz + obj.B.*obj.cxy;
            
            [DY,DZ] = ndgrid(obj.grid.Dy,obj.grid.Dz);
            vol = DY.*DZ;
            
            vecOut.yzDc = TScalar2D(obj.grid,'cell');
            vecOut.yzDc.v = obj.yzD; 
            vecOut.Ac = TScalar2D(obj.grid,'cell');
            vecOut.Ac.v = obj.A;              
            vecOut.Bc = TScalar2D(obj.grid,'cell');
            vecOut.Bc.v = obj.B;   
            
            vecOut.yyDc = TScalar2D(obj.grid,'cell');
            vecOut.yyDc.v = obj.yyD;              
            vecOut.zzDc = TScalar2D(obj.grid,'cell');
            vecOut.zzDc.v = obj.zzD;              
            
            temp = obj.V.*vol;
            vecOut.Vn = TScalar2D(obj.grid,'node');
            vecOut.Vn.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.Vn.v(2:end-1,2:end-1) = vecOut.Vn.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
            
            temp = obj.A.*vol;
            vecOut.An = TScalar2D(obj.grid,'node');
            vecOut.An.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.An.v(2:end-1,2:end-1) = vecOut.An.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));
            
            temp = obj.B.*vol;
            vecOut.Bn = TScalar2D(obj.grid,'node');
            vecOut.Bn.v(2:end-1,2:end-1) = temp(1:end-1,1:end-1)+temp(2:end,1:end-1)+...
                temp(1:end-1,2:end)+temp(2:end,2:end);
            vecOut.Bn.v(2:end-1,2:end-1) = vecOut.Bn.v(2:end-1,2:end-1)./(vol(1:end-1,1:end-1)+ ...
                vol(2:end,1:end-1)+vol(1:end-1,2:end)+vol(2:end,2:end));            
            
            temp = obj.A.*vol;
            vecOut.Am = TVector2D(obj.grid);
            vecOut.Am.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.Am.y(:,2:end-1) = vecOut.Am.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.Am.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.Am.z(2:end-1,:) = vecOut.Am.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.B.*vol;
            vecOut.Bm = TVector2D(obj.grid);
            vecOut.Bm.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.Bm.y(:,2:end-1) = vecOut.Bm.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.Bm.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.Bm.z(2:end-1,:) = vecOut.Bm.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.yyD.*vol;
            vecOut.yyDm = TVector2D(obj.grid);
            vecOut.yyDm.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.yyDm.y(:,2:end-1) = vecOut.yyDm.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.yyDm.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.yyDm.z(2:end-1,:) = vecOut.yyDm.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.yzD.*vol;
            vecOut.yzDm = TVector2D(obj.grid);
            vecOut.yzDm.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.yzDm.y(:,2:end-1) = vecOut.yzDm.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.yzDm.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.yzDm.z(2:end-1,:) = vecOut.yzDm.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
            
            temp = obj.zzD.*vol;
            vecOut.zzDm = TVector2D(obj.grid);
            vecOut.zzDm.y(:,2:end-1) = temp(:,1:end-1)+temp(:,2:end);
            vecOut.zzDm.y(:,2:end-1) = vecOut.zzDm.y(:,2:end-1)./(vol(:,1:end-1)+ ...
                vol(:,2:end));
            vecOut.zzDm.z(2:end-1,:) = temp(1:end-1,:)+temp(2:end,:);
            vecOut.zzDm.z(2:end-1,:) = vecOut.zzDm.z(2:end-1,:)./(vol(1:end-1,:)+ ...
                vol(2:end,:));
        end        
    end
end
