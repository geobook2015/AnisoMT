classdef TForwardModel2D_Aniso < TForwardModel2D
% Zeqiu Guo, Gary Egbert 2015
% zeqiu_guo@hotmail.com; egbert@coas.oregonstate.edu
% CUGB, Beijing China; OSU, OR USA
    properties
        xxm, xym, xye, xzm, xze, yye, zze, yze % vector of mapped model parameters
        Ax, Ay, Az
        Am,Bm,yyDm,yzDm,zzDm,Vn,An,Bn % vector of mapped model parameters
        Ae1,Ah1,Ae2,Ah2 % matrix for 2D Del^2 operator and first derivatives        
        yzc,yzDc,yyDc,zzDc,Ac,Bc
    end
    
    methods
        function obj = TForwardModel2D_Aniso(grid,m,omega)
            if nargin >= 1
                obj.grid = grid;
                obj.modOp = TModelOperator2D_Aniso(obj.grid);
            end
            
            if nargin >= 2                
                %     PIm is in this case conductivity mapped to edges ... but
                %     could be something different, such as resistivity mapped
                %     to faces.  PI(m) in the notation of Egbert and Kelbert
                %     (2012)
                obj.m = m;
            end
            
            if nargin ==3
                obj.omega = omega;
            end
        end
        
        function setEquations(obj)
            obj.modOp.setGrad;
            obj.modOp.setGradyz;
            obj.modOp.setMetricElements;
            obj.modOp.setIndicies;
        end
        
        function setFreqCond(obj,omega,m,style)
            %     This routine sets additional matrices which depend on
            %     conductivity, using matrices that are more generic
            %
            %     PIm is in this case conductivity mapped to edges ... but
            %     could be something different, such as resistivity mapped
            %     to faces.  PI(m) in the notation of Egbert and Kelbert
            %     (2012)
            
            %    for 2D version: factor matrix after forming
            
            if nargin >= 2
                obj.omega  = omega;
            end
            ne = length(obj.modOp.EdgeLength);
            ney = obj.grid.NEdges(1);
            nez = obj.grid.NEdges(2);
            nn = size(obj.modOp.NodeArea);
            nn = nn(1)*nn(2);
            nc = obj.grid.Ny*obj.grid.Nz;
            switch style
                case 'ThreeE'
                    % edge-to-node with G', node-to-edge with G,
                    % G'*G, nEdge*nEdge
                    if nargin == 4
                        vecOut = m.PDEmapping(style);
                        obj.xxm = vecOut.xxm.GetVector;
                        obj.xym = vecOut.xym.GetVector;
                        obj.xye = vecOut.xye.GetVector;
                        obj.xzm = vecOut.xzm.GetVector;
                        obj.xze = vecOut.xze.GetVector;
                        obj.yye = vecOut.yye.GetVector;
                        obj.zze = vecOut.zze.GetVector;
                        obj.yze = vecOut.yze.GetVector;
                        obj.yzc = vecOut.yzc.GetVector;
                    end                    
                
                    % Jx
                    dex =  spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G'...
                        *spdiags(obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*obj.modOp.G;
                    Vioms = obj.isign*1i*obj.omega*obj.mu*obj.xxm;
                    n = length(Vioms);
                    temp = dex + spdiags(Vioms,0,n,n);
                    
                    Vioms1 = obj.isign*1i*obj.omega*obj.mu*obj.xye(1:ney).*obj.modOp.EdgeLength(1:ney);
                    n = length(Vioms1);
                    temp1 = spdiags(reshape(repmat(1./obj.grid.DDy/2,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gya'*2*spdiags(Vioms1,0,n,n);
                    
                    Vioms2 = obj.isign*1i*obj.omega*obj.mu*obj.xze(ney+1:end).*obj.modOp.EdgeLength(ney+1:end);
                    n = length(Vioms2);
                    temp2 = spdiags(reshape(repmat(1./obj.grid.DDz'/2,obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gza'*2*spdiags(Vioms2,0,n,n);
                    
                    obj.Ax = [ temp temp1 temp2 ];
                    
                    % Jy
                    Vioms = obj.isign*1i*obj.omega*obj.mu*obj.xye(1:ney);
                    n = length(Vioms);
                    temp = spdiags(Vioms,0,n,n)*obj.modOp.Gya;
                    
                    dey = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1),0,ney,ney)*obj.modOp.Gef_zy'...
                        *spdiags(reshape(repmat(1./obj.grid.Dz',obj.grid.Ny,1),nc,1),0,nc,nc)*obj.modOp.Gef_zy;  
                    Vioms1 = obj.isign*1i*obj.omega*obj.mu*obj.yye(1:ney);
                    n = length(Vioms1);
                    temp1 = dey + spdiags(Vioms1,0,n,n);
                    

%                     dez =  spdiags(1./obj.modOp.EdgeLength(1:ney),0,ney,ney)*obj.modOp.Gy*spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)...
%                         *obj.modOp.Gz'*spdiags(obj.modOp.DualEdgeLength(ney+1:end),0,nez,nez);        
                    dez = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1),0,ney,ney)*obj.modOp.Gef_zy'...
                        *spdiags(reshape(repmat(1./obj.grid.Dy,1,obj.grid.Nz),nc,1),0,nc,nc)*obj.modOp.Gef_yz;
                    Vioms2 = obj.isign*1i*obj.omega*obj.mu*obj.yzc;
                    n = length(Vioms2);
                    temp2 = spdiags(reshape(repmat(1./obj.grid.DDz'/2,obj.grid.Ny,1),ney,1),0,ney,ney)*...
                        obj.modOp.Gefa_zy'*2*spdiags(Vioms2,0,n,n)*obj.modOp.Gefa_yz*spdiags(obj.modOp.EdgeLength(ney+1:end),0,nez,nez);
                    temp2 = temp2 - dez;
                    
                    obj.Ay = [ temp temp1 temp2 ];
                    
                    % Jz
                    Vioms = obj.isign*1i*obj.omega*obj.mu*obj.xze(ney+1:end);
                    n = length(Vioms);
                    temp =spdiags(Vioms,0,n,n)*obj.modOp.Gza;
                    
%                     dey =  spdiags(1./obj.modOp.EdgeLength(ney+1:end),0,nez,nez)*obj.modOp.Gz*spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)...
%                         *obj.modOp.Gy'*spdiags(obj.modOp.DualEdgeLength(1:ney),0,ney,ney);           
                    dey = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1),0,nez,nez)*obj.modOp.Gef_yz'...
                        *spdiags(reshape(repmat(1./obj.grid.Dz',obj.grid.Ny,1),nc,1),0,nc,nc)*obj.modOp.Gef_zy;  
                    Vioms1 = obj.isign*1i*obj.omega*obj.mu*obj.yzc;
                    n = length(Vioms1);
                    temp1 = spdiags(reshape(repmat(1./obj.grid.DDy/2,1,obj.grid.Nz),nez,1),0,nez,nez)*...
                        obj.modOp.Gefa_yz'*2*spdiags(Vioms1,0,n,n)*obj.modOp.Gefa_zy*spdiags(obj.modOp.EdgeLength(1:ney),0,ney,ney);
                    temp1 = temp1 - dey;
                    
                    dez = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1),0,nez,nez)*obj.modOp.Gef_yz'...
                        *spdiags(reshape(repmat(1./obj.grid.Dy,1,obj.grid.Nz),nc,1),0,nc,nc)*obj.modOp.Gef_yz;
                    Vioms2 = obj.isign*1i*obj.omega*obj.mu*obj.zze(ney+1:end);
                    n = length(Vioms2);
                    temp2 = dez + spdiags(Vioms2,0,n,n);
                    
                    obj.Az = [ temp temp1 temp2 ];
                    
                    obj.A = [ obj.Ax; obj.Ay; obj.Az ];                    
     
                case 'EHSta'
                    % edge-to-node with G', node-to-edge with G,
                    % G'*G, nEdge*nEdge
                    if nargin == 4
                        vecOut = m.PDEmapping(style);
                        obj.Vn = vecOut.Vn.GetVector;
                        obj.An = vecOut.An.GetVector;
                        obj.Bn = vecOut.Bn.GetVector;
                        obj.Am = vecOut.Am.GetVector;
                        obj.Bm = vecOut.Bm.GetVector;
                        obj.yyDm = vecOut.yyDm.GetVector;
                        obj.yzDm = vecOut.yzDm.GetVector;
                        obj.zzDm = vecOut.zzDm.GetVector;
                        obj.yzDc = vecOut.yzDc.GetVector;   
                        obj.yyDc = vecOut.yyDc.GetVector;
                        obj.zzDc = vecOut.zzDc.GetVector;                          
                    end
                    
                    % E-mode discritization
                    
                    obj.Ae1 =  spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G'...
                        *spdiags(obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*obj.modOp.G;
                    Vioms = obj.isign*1i*obj.omega*obj.mu*obj.Vn;
                    n = length(Vioms);
                    obj.Ae1 = obj.Ae1 + spdiags(Vioms,0,n,n);
                    
                    % number of DualEdgeLength in y direction equals to number of Z
                    % edges, and vice versa
                    
%                     temp1 = obj.isign*1i*obj.omega*obj.mu.*obj.Am(ney+1:end).*obj.modOp.EdgeLength(ney+1:end).*reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1);
%                     Ahy = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gza'...
%                         *spdiags(temp1,0,nez,nez)*obj.modOp.Gef_yz';
                    
                    temp1 =obj.modOp.EdgeLength(ney+1:end).*reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1);
                    Ahy = spdiags(obj.isign*1i*obj.omega*obj.mu.*obj.An.*reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gza'...
                        *spdiags(temp1,0,nez,nez)*obj.modOp.Gef_yz';


%                     Ahz = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gya'*spdiags(...
%                         obj.isign*1i*obj.omega*obj.mu.*obj.Bm(1:ney).*obj.modOp.EdgeLength(1:ney),0,ney,ney)...
%                         *spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1),0,ney,ney)*obj.modOp.Gef_zy';
                    Ahz = spdiags(obj.isign*1i*obj.omega*obj.mu.*obj.Bn.*reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gya'*...
                        spdiags(obj.modOp.EdgeLength(1:ney).*reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1),0,ney,ney)*obj.modOp.Gef_zy';
                    obj.Ah1 = -Ahy + Ahz; % cause G'*G corresponds to  - Delta^2
                    
                    % H-mode discritization
                    
                    temp1 = 1./obj.modOp.EdgeLength(ney+1:end);
                    Aez = spdiags(temp1,0,nez,nez)*obj.modOp.Gz*spdiags(obj.Bn,0,obj.grid.NNodes,obj.grid.NNodes);
                    Aez = obj.modOp.Gefa_yz*Aez;
                    
                    temp2 = 1./obj.modOp.EdgeLength(1:ney);
                    Aey = spdiags(temp2,0,ney,ney)*obj.modOp.Gy*spdiags(obj.An,0,obj.grid.NNodes,obj.grid.NNodes);
                    Aey = obj.modOp.Gefa_zy*Aey;
                    obj.Ae2 = - Aez + Aey; % cause G'*G corresponds to  - Delta^2
                    
                    temp1 = obj.yyDm(ney+1:end).*reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1);
                    Ah2t1 = spdiags(reshape(repmat(1./obj.grid.Dy,1,obj.grid.Nz),nc,1),0,nc,nc)*obj.modOp.Gef_yz*spdiags(temp1,0,nez,nez)*obj.modOp.Gef_yz';
                    
                    temp2 = obj.zzDm(1:ney).*reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1);
                    Ah2t2 = spdiags(reshape(repmat(1./obj.grid.Dz',obj.grid.Ny,1),nc,1),0,nc,nc)*obj.modOp.Gef_zy*spdiags(temp2,0,ney,ney)*obj.modOp.Gef_zy';
                    
                    temp3 = obj.yzDc.*reshape(repmat(1./obj.grid.Dz',obj.grid.Ny,1),nc,1);
                    Ah2t3 = obj.modOp.Gefa_yz*spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1),0,nez,nez)*...
                        obj.modOp.Gef_yz'*spdiags(temp3,0,nc,nc)*obj.modOp.Gef_zy*obj.modOp.Gefa_zy';
                    
                    temp4 = obj.yzDc.*reshape(repmat(1./obj.grid.Dy,1,obj.grid.Nz),nc,1);
                    Ah2t4 = obj.modOp.Gefa_zy*spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1),0,ney,ney)*...
                        obj.modOp.Gef_zy'*spdiags(temp4,0,nc,nc)*obj.modOp.Gef_yz*obj.modOp.Gefa_yz';
                    
                    % keep coefficients of boundary nodes as ones,
                    % because they are zeros from discretization of del^2
                    % for H mode;
                    cof = zeros(obj.grid.Ny,obj.grid.Nz); % ones, not zeros-->if Direclet BCs used
                    cof(:,[1 end]) = 1;
                    cof(2:end-1,2:end-1) = obj.isign * 1i*obj.omega*obj.mu;
                    Vioms = reshape(cof,numel(cof),1);
                    %             Vioms = obj.isign*1i*obj.omega*obj.mu*ones(obj.grid.NCells,1);
                    n = length(Vioms);
                    obj.Ah2 = Ah2t1 + Ah2t2 + Ah2t3 + Ah2t4 + spdiags(Vioms,0,n,n);
                    
%                     % newmann/second-type BCs for Hx at left- and right-boundaries
%                     gr = obj.grid;
%                     Ahb = zeros(nc);
%                     J=ones(1,gr.Nz);
%                     K=1:gr.Nz;
%                     indhL=gr.vectorIndex(J,K,'cell');
%                     J=ones(1,gr.Nz)*gr.Ny;
%                     K=1:gr.Nz;
%                     indhR=gr.vectorIndex(J,K,'cell');
%                     %             J=2:gr.Ny-1;
%                     %             K=ones(1,gr.Ny-2);
%                     %             indhT=gr.vectorIndex(J,K,'cell');
%                     for i=1:length(indhL)
%                         Ahb(indhL(i),indhL(i)) = -1;
%                         Ahb(indhL(i),indhL(i)+1) = 1;
%                     end
%                     for i=1:length(indhR)
%                         Ahb(indhR(i),indhR(i)) = 1;
%                         Ahb(indhR(i),indhR(i)-1) = -1;
%                     end
%                     %             for i=1:length(indhT)
%                     %                 Ahb(indhT(i),indhT(i)) = -1;
%                     %                 Ahb(indhT(i),indhT(i)+gr.Ny) = 1;
%                     %             end
%                     obj.Ah2 = obj.Ah2 + sparse(Ahb);
                    
                    obj.A = [ obj.Ae1 obj.Ah1; obj.Ae2 obj.Ah2 ];
                case 'EHFix'
                    if nargin == 4
                        vecOut = m.PDEmapping(style);
                        obj.Vn = vecOut.Vn.GetVector;
                        obj.An = vecOut.An.GetVector;
                        obj.Bn = vecOut.Bn.GetVector;
                        obj.Am = vecOut.Am.GetVector;
                        obj.Bm = vecOut.Bm.GetVector;
                        obj.yyDm = vecOut.yyDm.GetVector;
                        obj.yzDm = vecOut.yzDm.GetVector;
                        obj.zzDm = vecOut.zzDm.GetVector;
                        obj.yzDc = vecOut.yzDc.GetVector;   
                        obj.yyDc = vecOut.yyDc.GetVector;
                        obj.zzDc = vecOut.zzDc.GetVector;       
                        obj.Ac = vecOut.Ac.GetVector;                          
                        obj.Bc = vecOut.Bc.GetVector;                          
                    end
                    
                    % E-mode discretization                    
                    obj.Ae1 =  spdiags(1./reshape(obj.modOp.NodeArea,nn,1),0,nn,nn)*obj.modOp.G'...
                        *spdiags(obj.modOp.DualEdgeLength./obj.modOp.EdgeLength,0,ne,ne)*obj.modOp.G;
                    Vioms = obj.isign * 1i*obj.omega*obj.mu*obj.Vn;
                    n = length(Vioms);
                    obj.Ae1 = obj.Ae1 + spdiags(Vioms,0,n,n);
                    
                    temp1 = obj.isign*1i*obj.omega*obj.mu.*obj.Am(1:ney);
                    Ahy = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gya'*...
                        spdiags(temp1,0,ney,ney)*obj.modOp.Gy;
                    temp2 = obj.isign*1i*obj.omega*obj.mu.*obj.Bm(ney+1:end);
                    Ahz = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gza'*...
                        spdiags(temp2,0,nez,nez)*obj.modOp.Gz;

                    obj.Ah1 = Ahy-Ahz;
                    
                    % H-mode discretization
                    Aez = obj.modOp.Gz*spdiags(obj.Bn,0,obj.grid.NNodes,obj.grid.NNodes);
                    Aez = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*...
                        obj.modOp.Gza'*Aez;
                    
                    Aey = obj.modOp.Gy*spdiags(obj.An,0,obj.grid.NNodes,obj.grid.NNodes);
                    Aey = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*...
                        obj.modOp.Gya'*Aey;

%                     Aez = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*...
%                         obj.modOp.Gz'*spdiags(obj.Bm(ney+1:end),0,nez,nez)*obj.modOp.Gza;
% 
%                     Aey = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*...
%                         obj.modOp.Gy'*spdiags(obj.Am(1:ney),0,ney,ney)*obj.modOp.Gya;
                    
                    obj.Ae2 = Aey-Aez;
                    
                
                    Ah2t1 = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gy'*...
                        spdiags(obj.yyDm(1:ney)./obj.modOp.EdgeLength(1:ney),0,ney,ney)*obj.modOp.Gy;
                    
                    Ah2t2 = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gz'*...
                        spdiags(obj.zzDm(ney+1:end)./obj.modOp.EdgeLength(ney+1:end),0,nez,nez)*obj.modOp.Gz;
                    
                    %spdiags(reshape(repmat(1./obj.grid.Dy,1,obj.grid.Nz),nc,1),0,nc,nc)*...
                    
                    Ah2t3 = spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny+1,1),nn,1),0,nn,nn)*obj.modOp.Gza'*...
                        spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz),nez,1).*obj.modOp.EdgeLength(ney+1:end),0,nez,nez)*obj.modOp.Gef_yz'...
                        *spdiags(obj.yzDc,0,nc,nc)*obj.modOp.Gefa_yz*spdiags(1./obj.modOp.EdgeLength(ney+1:end),0,nez,nez)*obj.modOp.Gz;
                    
                    Ah2t4 = spdiags(reshape(repmat(1./obj.grid.DDy,1,obj.grid.Nz+1),nn,1),0,nn,nn)*obj.modOp.Gya'*...
                        spdiags(reshape(repmat(1./obj.grid.DDz',obj.grid.Ny,1),ney,1).*obj.modOp.EdgeLength(1:ney),0,ney,ney)*obj.modOp.Gef_zy'...
                        *spdiags(obj.yzDc,0,nc,nc)*obj.modOp.Gefa_zy*spdiags(1./obj.modOp.EdgeLength(1:ney),0,ney,ney)*obj.modOp.Gy;
                    
                    % keep coefficients of boundary nodes as ones,
                    % because they are their from discretization of del^2
                    % for H mode;
                    cof = zeros(obj.grid.Ny+1,obj.grid.Nz+1);
                    cof(:, [1 end])=1;
                    cof(2:end-1,2:end-1) = obj.isign * 1i*obj.omega*obj.mu;
                    Vioms = reshape(cof,numel(cof),1);
                    
                    %Vioms = obj.isign * 1i*obj.omega*obj.mu*ones(obj.grid.NNodes);
                    n = length(Vioms);
                    obj.Ah2 = Ah2t1 + Ah2t2 + Ah2t3 + Ah2t4 + spdiags(Vioms,0,n,n);
                    
%                     % newmann/second-type BCs for Hx at left- and right-boundaries
%                     gr = obj.grid;
%                     Ahb = zeros(gr.NNodes);
%                     J=ones(1,gr.Nz-1);
%                     K=2:gr.Nz;
%                     indhL=gr.vectorIndex(J,K,'node');
%                     J=ones(1,gr.Nz-1)*gr.Ny;
%                     K=2:gr.Nz;
%                     indhR=gr.vectorIndex(J,K,'node');
%                     %             J=2:gr.Ny-1;
%                     %             K=ones(1,gr.Ny-2);
%                     %             indhT=gr.vectorIndex(J,K,'cell');
%                     for i=1:length(indhL)
%                         Ahb(indhL(i),indhL(i)) = -1;
%                         Ahb(indhL(i),indhL(i)+1) = 1;
%                     end
%                     for i=1:length(indhR)
%                         Ahb(indhR(i),indhR(i)) = 1;
%                         Ahb(indhR(i),indhR(i)-1) = -1;
%                     end
%                     %             for i=1:length(indhT)
%                     %                 Ahb(indhT(i),indhT(i)) = -1;
%                     %                 Ahb(indhT(i),indhT(i)+gr.Ny) = 1;
%                     %             end
% %                     obj.Ah2 = obj.Ah2 + sparse(Ahb);
                    
                    obj.A = [ obj.Ae1 obj.Ah1; obj.Ae2 obj.Ah2 ];
                otherwise
                    error('The style you specified has not been implemented');
            end
        end
        function [ e,inorm ] = Solve(obj,src,style)
            %  solve equations for specified source
            %    Source is an input source object defining all source terms
            %      including BC.   It is an instance of the abstract
            %      TSource class
            %   The output is a column vector, but of full grid dimension,
            %   including boundary data
            src.setRHS(style);   % this completes definition of source terms
            ii = src.ii;
            ib = src.ib;
%             ii = obj.modOp.ie;
%             ib = obj.modOp.be;
            e = src.B;
            
            %   need to multiply source by volume elements
            rhs = e(ii)-obj.A(ii,ib)*e(ib);
%             rhs = -obj.A(ii,ib)*e(ib);
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
            display(inorm);
        end
    end
end
