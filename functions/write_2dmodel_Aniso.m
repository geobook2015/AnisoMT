function write_2_model_Aniso
%revised to write 2D data body for 1D layer model
% Zeqiu, 8/6/2015

%function status = write_WS3d_model(fname,x,y,z,rho,type,origin,rotation)
%revised to write generally anisotropic model files 
%for 1D layer model temporarily, take considerations from Wilson and Thiel (2002) 
%about rotation of resistivity tensor and conductivity tensor
% Zeqiu Guo, 28/4/2015

% writes a 3D resistivity model in Weerachai Siripunvaraporn's format;
% allows for natural log resistivity by setting type = 'LOGE'
%  (c) Anna Kelbert, 2009
%  open file for output
clear;clc;
close all;
NL=4;
res=[10000.  10000.  10000.;
    300.  10000.  300;
    30.   1000.   30.;
    30.  30.    30.];
ustr=[0 30 60 0];
udip=[0 0 0 0];
usla=[0 0 0 0];
sg(NL,3,3)=0;
con(NL,3)=0;
for layer=1:NL    
    for j=1:3
        con(layer,j)=1./res(layer,j);
    end    
    sig_tem=spdiags(con(layer,:)',0,3,3);
    rstr=pi*ustr(layer)/180;
    rdip=pi*udip(layer)/180;
    rsla=pi*usla(layer)/180;
    Rsz_r=[cos(rstr) sin(rstr) 0; -sin(rstr) cos(rstr) 0; 0 0 1];
    Rsz_p=[cos(-rstr) sin(-rstr) 0; -sin(-rstr) cos(-rstr) 0; 0 0 1];
    Rdx_r=[1 0 0; 0 cos(rdip) sin(rdip); 0 -sin(rdip) cos(rdip)];
    Rdx_p=[1 0 0; 0 cos(-rdip) sin(-rdip); 0 -sin(-rdip) cos(-rdip)];
    Rlz_r=[cos(rsla) sin(rsla) 0; -sin(rsla) cos(rsla) 0; 0 0 1];
    Rlz_p=[cos(-rsla) sin(-rsla) 0; -sin(-rsla) cos(-rsla) 0; 0 0 1];    
    sg(layer,:,:)=Rsz_p*Rdx_p*Rlz_p*sig_tem*Rlz_r*Rdx_r*Rsz_r;
end
Thick=[10 30 160 250]'*1000; %Thickness/(m)
Thk=[250*ones(1,40) 500*ones(1,60) 1000*ones(1,160) 12500*ones(1,20)];
dep=cumsum(Thick);
Axx=zeros(length(Thk),1);
Ayy=zeros(length(Thk),1);
Azz=zeros(length(Thk),1);
Axy=zeros(length(Thk),1);
Ayz=zeros(length(Thk),1);
Axz=zeros(length(Thk),1);
for i=1:NL
    if i~=1
        Axx(dep(i-1)<cumsum(Thk))=sg(i,1,1);
        Ayy(dep(i-1)<cumsum(Thk))=sg(i,2,2);
        Azz(dep(i-1)<cumsum(Thk))=sg(i,3,3);
        Axy(dep(i-1)<cumsum(Thk))=sg(i,1,2);
        Ayz(dep(i-1)<cumsum(Thk))=sg(i,2,3);
        Axz(dep(i-1)<cumsum(Thk))=sg(i,1,3);        
    elseif i==1
        Axx(cumsum(Thk)<=dep(1))=sg(1,1,1);
        Ayy(cumsum(Thk)<=dep(1))=sg(1,2,2);
        Azz(cumsum(Thk)<=dep(1))=sg(1,3,3);
        Axy(cumsum(Thk)<=dep(1))=sg(1,1,2);
        Ayz(cumsum(Thk)<=dep(1))=sg(1,2,3);
        Axz(cumsum(Thk)<=dep(1))=sg(1,1,3);
    end
end
rho=[Axx Ayy Azz Axy Ayz Axz];

y=[20000 10000 5000 2500 1250 1250 1250 1250 2500 5000 5000 2500 1250 1250 1250 ...
    1250 2500 5000 5000 2500 1250 1250 1250 1250 2500 5000 10000 20000]/5;
y=reshape(repmat(y,5,1),5*length(y),1);

z=Thk;
fid = fopen('../Test/Modelfile_Aniso','w');
ny=length(y);
nz=length(z);
type = ' ';
comment = 'Written by Matlab write_2dmodel_Aniso script     --Anisotropy--'; 
fprintf(fid, '# %s\n', comment);
fprintf(fid, '%d %d %d %s\n', ny, nz, 0, type);
for j = 1:ny
    status = fprintf(fid,'%G ',y(j));
end
fprintf(fid, '\n');
for j = 1:nz
    status = fprintf(fid,'%G ',z(j));
end
fprintf(fid, '\n');
for m=1:6
    for k = 1:nz
        fprintf(fid,'\n');
        for j = 1:ny
                fprintf(fid,'%17.15E ',rho(k,m));
         end
    end
end
fclose(fid);
return
