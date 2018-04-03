%%   Driver script for running MT 2D anisotropic modeling
% load in model parameter, and grid
clc;clear;close all
% modelfile
modelFile = 'Modelfile_Aniso_new'
% method
style = 'ThreeE'  % 'ThreeE'  'EHSta'  'EHFix'
% period
T = logspace(-2,4,15);
% polarisation
XY={'X','Y'};
%   this block is specific to a particular model parameter implementation
m = TModelParameterCell2D_Aniso();
m = m.readVec(modelFile);
grid = TGrid2D(m.grid.Dy,m.grid.Dz,m.grid.Nza);
m.grid=grid;
m.ParamGrid = grid;
%%
%    create forward modeling operator: probably in normal usage just do grid
%    dependent setup; setup for model parameter and frequency separately

fwd  = TForwardModel2D_Aniso(grid,m);
fwd.setEquations;
Z_imp(2,2,length(T),grid.Ny)=0;
rxx(length(T),grid.Ny)=0;
ryy(length(T),grid.Ny)=0;
rxy(length(T),grid.Ny)=0;
ryx(length(T),grid.Ny)=0;
pxx(length(T),grid.Ny)=0;
pyy(length(T),grid.Ny)=0;
pxy(length(T),grid.Ny)=0;
pyx(length(T),grid.Ny)=0;
Exs{2}=0; Eys{2}=0;
Hxs{2}=0; Hys{2}=0;
ney = grid.NEdges(1);
nez = grid.NEdges(2);
nn = grid.NNodes;
nc = grid.Ny*grid.Nz;
%   set up source term
src = TSource2DMT_Aniso_1DBC(grid,m);
%  ultimately make a loop over periods and modes
for j=1:length(T)
    for i = 1:2
        pol = XY{i};
        omega = 2*pi/T(j);   %  angular frequency
        switch style
            case 'ThreeE'
                fwd.setFreqCond(omega,m,style); % 'ThreeE'  'EHSta'  'EHFix'
                src.SetSourceParams(omega,pol);
                %   solve
                e = fwd.Solve(src,style); % 'ThreeE'  'EHSta'  'EHFix'
                
%                 % Debug: deliver boundary values into interior region %
%                 e = [ reshape(repmat(src.E1D.F.Ex,1,grid.Ny+1)',(grid.Ny+1)*(grid.Nz+1),1);...
%                     reshape(repmat(src.E1D.F.Ey,1,grid.Ny)',(grid.Ny)*(grid.Nz+1),1);...
%                     reshape(repmat(src.E1D.F.Ez,1,grid.Ny+1)',(grid.Ny+1)*(grid.Nz),1) ];
%                 %-----------------------------------------------------%
                
                Ex = reshape(e(1:grid.NNodes),grid.Ny+1,grid.Nz+1);
                Ey = reshape(e(grid.NNodes+1:grid.NNodes+grid.NEdges(1)),grid.Ny,grid.Nz+1);
                Ez = reshape(e(grid.NNodes+grid.NEdges(1)+1:end),grid.Ny+1,grid.Nz);                
                Hy = 1/(-fwd.isign*1i*omega*fwd.mu)*spdiags(1./fwd.modOp.EdgeLength(ney+1:end),0,nez,nez)*fwd.modOp.Gz*e(1:grid.NNodes);
                Hz = -1/(fwd.isign*1i*omega*fwd.mu)*spdiags(1./fwd.modOp.EdgeLength(1:ney),0,ney,ney)*fwd.modOp.Gy*e(1:grid.NNodes);
                Hx = 1/(-fwd.isign*1i*omega*fwd.mu)*(spdiags(reshape(repmat(1./grid.Dy,1,grid.Nz),nc,1),0,nc,nc)*...
                    fwd.modOp.Gef_yz*e(grid.NNodes+grid.NEdges(1)+1:end) - spdiags(reshape(repmat(1./grid.Dz',grid.Ny,1),...
                    nc,1),0,nc,nc)*fwd.modOp.Gef_zy*e(grid.NNodes+1:grid.NNodes+grid.NEdges(1)));
                Hx = reshape(Hx, grid.Ny, grid.Nz);
                Hy = reshape(Hy, grid.Ny+1, grid.Nz);
                Hz = reshape(Hz, grid.Ny, grid.Nz+1);
                
                Exs{i} = (Ex(1:end-1,grid.Nza+1)+Ex(2:end,grid.Nza+1))/2;
                Eys{i} = Ey(:,grid.Nza+1);
                
%                                 Hxs{i} = (Hx(:,grid.Nza).*grid.Dz(grid.Nza)+Hx(:,grid.Nza+1).*grid.Dz(grid.Nza+1))/grid.DDz(grid.Nza+1)/2;
                                
                x0 = grid.Dz(grid.Nza+1)/2;  x1 = grid.Dz(grid.Nza+1)+grid.Dz(grid.Nza+2)/2;...
                x2 = sum(grid.Dz(grid.Nza+1:grid.Nza+2))+grid.Dz(grid.Nza+3)/2;
                Hxs{i} =  Hx(:,grid.Nza+1)*x1*x2/(x0-x1)/(x0-x2) + Hx(:,grid.Nza+2)*x0*x2/(x1-x0)/(x1-x2) + Hx(:,grid.Nza+3)*x0*x1/(x2-x0)/(x2-x1);
                
                temp1 = (Hy(1:end-1,grid.Nza)+Hy(2:end,grid.Nza))/2;
                temp2 = (Hy(1:end-1,grid.Nza+1)+Hy(2:end,grid.Nza+1))/2;                
                Hys{i} = (temp1.*grid.Dz(grid.Nza+1)+temp2.*grid.Dz(grid.Nza))/grid.DDz(grid.Nza+1)/2;
                
%                 Hys{i} = Hy(:,grid.Nza+1)+(Hy(:,grid.Nza+1)-Hy(:,grid.Nza+2))*grid.Dz(grid.Nza+1)/(grid.Dz(grid.Nza+2)+grid.Dz(grid.Nza+1));        

%                 temp1 = (Hy(1:end-1,grid.Nza+1)+Hy(2:end,grid.Nza+1))/2;
%                 temp2 = (Hy(1:end-1,grid.Nza+2)+Hy(2:end,grid.Nza+2))/2;         
%                 temp3 = (Hy(1:end-1,grid.Nza+3)+Hy(2:end,grid.Nza+3))/2;            
%                 Hys{i} =  temp1*x1*x2/(x0-x1)/(x0-x2) + temp2*x0*x2/(x1-x0)/(x1-x2) + temp3*x0*x1/(x2-x0)/(x2-x1);                
            case 'EHSta'
                fwd.setFreqCond(omega,m,style);
                src.SetSourceParams(omega,pol);
                %   solve
                e = fwd.Solve(src,style);
                Ex = reshape(e(1:grid.NNodes),grid.Ny+1,grid.Nz+1);
                Hx = reshape(e(grid.NNodes+1:end),grid.Ny,grid.Nz);
                Exs{i} = (Ex(1:end-1,grid.Nza+1)+Ex(2:end,grid.Nza+1))/2;
%                 Hxs{i} = (Hx(:,grid.Nza).*grid.Dz(grid.Nza)+Hx(:,grid.Nza+1).*grid.Dz(grid.Nza+1))/grid.DDz(grid.Nza+1)/2;
                x0 = grid.Dz(grid.Nza+1)/2;  x1 = grid.Dz(grid.Nza+1)+grid.Dz(grid.Nza+2)/2;...
                    x2 = sum(grid.Dz(grid.Nza+1:grid.Nza+2))+grid.Dz(grid.Nza+3)/2;
                Hxs{i} =  Hx(:,grid.Nza+1)*x1*x2/(x0-x1)/(x0-x2) + Hx(:,grid.Nza+2)*x0*x2/(x1-x0)/(x1-x2) + Hx(:,grid.Nza+3)*x0*x1/(x2-x0)/(x2-x1);                
                
                %%%-------------------------------------------------------------------
                % shouldn't let conductivity of boudary cells be zeros, different from
                % settings of BCs; and here is a simplified test, should do V*c/V average
                Bn = reshape(fwd.Bn, grid.Ny+1, grid.Nz+1);
                Bn(1,:) = Bn(2,:); Bn(end,:) = Bn(end-1,:);
                fwd.Bn = reshape(Bn,(grid.Ny+1)*(grid.Nz+1),1);
                
                yzDm = reshape(fwd.yzDm(1:ney), grid.Ny, grid.Nz+1);
                yzDm(1,:) = yzDm(2,:); yzDm(end,:) = yzDm(end-1,:);
                fwd.yzDm(1:ney) = reshape(yzDm,grid.Ny*(grid.Nz+1),1);
                
                zzDm = reshape(fwd.zzDm(ney+1:end), grid.Ny+1, grid.Nz);
                zzDm(1,:) = zzDm(2,:); zzDm(end,:) = zzDm(end-1,:);
                fwd.zzDm(ney+1:end) = reshape(zzDm,(grid.Ny+1)*grid.Nz,1);
                %%%-------------------------------------------------------------------
                
                % adding minus sign is because of transpose of derivative operator, not average operator
                temp1 = -fwd.modOp.Gya*spdiags(reshape(repmat(1./grid.DDz',grid.Ny+1,1),nn,1),0,nn,nn)*...
                    fwd.modOp.Gza'*spdiags(fwd.modOp.EdgeLength(ney+1:end).*fwd.yzDm(ney+1:end).*reshape(repmat(1./grid.DDy,1, ...
                    grid.Nz),nez,1),0,nez,nez)*fwd.modOp.Gef_yz'*e(grid.NNodes+1:end) + ...
                    -spdiags(fwd.zzDm(1:ney).*reshape(repmat(1./grid.DDz',grid.Ny,1),...
                    ney,1),0,ney,ney)*fwd.modOp.Gef_zy'*e(grid.NNodes+1:end) + ...
                    fwd.modOp.Gya*(e(1:grid.NNodes).*fwd.Bn);
                temp1 = reshape(temp1,grid.Ny,grid.Nz+1);
                
                temp2 = spdiags(1/(-fwd.isign*1i*omega*fwd.mu)./fwd.modOp.EdgeLength(ney+1:end),0,nez,nez)...
                    *fwd.modOp.Gz*e(1:grid.NNodes);
                temp2 = fwd.modOp.Gefa_yz*temp2;
                temp2 = reshape(temp2,grid.Ny,grid.Nz);
                
                Eys{i} = temp1(:,grid.Nza+1);
                %         Hys{i} = temp2(:,grid.Nza);
%                 Hys{i} = (temp2(:,grid.Nza).*grid.Dz(grid.Nza)+temp2(:,grid.Nza+1).*grid.Dz(grid.Nza+1))/grid.DDz(grid.Nza+1)/2;
                Hys{i} =  temp2(:,grid.Nza+1)*x1*x2/(x0-x1)/(x0-x2) + temp2(:,grid.Nza+2)*x0*x2/(x1-x0)/(x1-x2) + temp2(:,grid.Nza+3)*x0*x1/(x2-x0)/(x2-x1);                       
            case 'EHFix'
                fwd.setFreqCond(omega,m,style);
                src.SetSourceParams(omega,pol);
                %   solve
                e = fwd.Solve(src,style);
                Ex = reshape(e(1:grid.NNodes),grid.Ny+1,grid.Nz+1);
                Hx = reshape(e(grid.NNodes+1:end),grid.Ny+1,grid.Nz+1);
                Exs{i} = (Ex(1:end-1,grid.Nza+1)+Ex(2:end,grid.Nza+1))/2;
                Hxs{i} = (Hx(1:end-1,grid.Nza+1)+Hx(2:end,grid.Nza+1))/2;
                
                %%%-------------------------------------------------------------------
                % shouldn't let conductivity of boudary cells be zeros, different from
                % settings of BCs; and here is a simplified test, should do V*c/V average
                Bn = reshape(fwd.Bn, grid.Ny+1, grid.Nz+1);
                Bn(1,:) = Bn(2,:); Bn(end,:) = Bn(end-1,:);
                fwd.Bn = reshape(Bn,(grid.Ny+1)*(grid.Nz+1),1);
                
                yzDm = reshape(fwd.yzDm(1:ney), grid.Ny, grid.Nz+1);
                yzDm(1,:) = yzDm(2,:); yzDm(end,:) = yzDm(end-1,:);
                fwd.yzDm(1:ney) = reshape(yzDm,grid.Ny*(grid.Nz+1),1);
                
                zzDm = reshape(fwd.zzDm(ney+1:end), grid.Ny+1, grid.Nz);
                zzDm(1,:) = zzDm(2,:); zzDm(end,:) = zzDm(end-1,:);
                fwd.zzDm(ney+1:end) = reshape(zzDm,(grid.Ny+1)*grid.Nz,1);
                %%%-------------------------------------------------------------------
                
                temp1 = spdiags(fwd.yzDm(1:ney)./fwd.modOp.EdgeLength(1:ney),0,ney,ney)*fwd.modOp.Gy*e(grid.NNodes+1:end) + ...
                    fwd.modOp.Gya*spdiags(reshape(repmat(1./grid.DDz',grid.Ny+1,1),nn,1),0,nn,nn)*fwd.modOp.Gza'*...
                    spdiags(fwd.zzDm(ney+1:end),0,nez,nez)*fwd.modOp.Gz*e(grid.NNodes+1:end)+ ... 
                    fwd.modOp.Gefa_zy'*fwd.modOp.Gefa_yz*fwd.modOp.Gza*(e(1:grid.NNodes).*fwd.Bn);
                temp1 = reshape(temp1,grid.Ny,grid.Nz+1);
                
%                 temp2 = spdiags(1/(-fwd.isign*1i*omega*fwd.mu)./fwd.modOp.EdgeLength(ney+1:end),0,nez,nez)...
%                     *fwd.modOp.Gz*e(1:grid.NNodes);
%                 temp2 = spdiags(reshape(repmat(1./grid.DDz',grid.Ny+1,1),nn,1),0,nn,nn)...
%                     *fwd.modOp.Gza'*spdiags(fwd.modOp.EdgeLength(ney+1:end),0,nez,nez)*temp2;
%                 temp2 = reshape(temp2,grid.Ny+1,grid.Nz+1);

                temp2 = fwd.modOp.Gefa_yz*spdiags(1/(-fwd.isign*1i*omega*fwd.mu)./fwd.modOp.EdgeLength(ney+1:end),0,nez,nez)...
                    *fwd.modOp.Gz*e(1:grid.NNodes);
                temp2 = reshape(temp2,grid.Ny,grid.Nz);
                
                Eys{i} = temp1(:,grid.Nza+1);
                
%                 Hys{i} = (temp2(1:end-1,grid.Nza+1)+temp2(2:end,grid.Nza+1))/2;
                
                Hys{i} = (temp2(:,grid.Nza).*grid.Dz(grid.Nza+1)+temp2(:,grid.Nza+1).*grid.Dz(grid.Nza))/grid.DDz(grid.Nza+1)/2;
            otherwise
                error('The style you specified has not been implemented')
        end
    end
    for k = 1:grid.Ny
        Z_imp(:,:,j,k) = [ Exs{1}(k) Exs{2}(k); Eys{1}(k) Eys{2}(k) ]/[ Hxs{1}(k) Hxs{2}(k); Hys{1}(k) Hys{2}(k) ];
    end
    rxx(j,:) = abs(Z_imp(1,1,j,:)).^2/(omega*fwd.mu);
    ryy(j,:) = abs(Z_imp(2,2,j,:)).^2/(omega*fwd.mu);
    rxy(j,:) = abs(Z_imp(1,2,j,:)).^2/(omega*fwd.mu);
    ryx(j,:) = abs(Z_imp(2,1,j,:)).^2/(omega*fwd.mu);
    
    pxx(j,:) = atan2(imag(Z_imp(1,1,j,:)),real(Z_imp(1,1,j,:)))*180/pi;
    pyy(j,:) = atan2(imag(Z_imp(2,2,j,:)),real(Z_imp(2,2,j,:)))*180/pi;
    pxy(j,:) = atan2(imag(Z_imp(1,2,j,:)),real(Z_imp(1,2,j,:)))*180/pi;
    pyx(j,:) = atan2(imag(Z_imp(2,1,j,:)),real(Z_imp(2,1,j,:)))*180/pi;
end