function [y,z,rho,nzAir,type,par] = read_2dmodel_Aniso(fname)
% reads an anistropic 2D conductivity model

fid = fopen(fname);
Aniso = false;
line = fgetl(fid); % comment
if strfind(line,'Anisotropy')
    Aniso = true;
end
line = fgetl(fid);
[n] = sscanf(line,'%d',[3 1]);
if strfind(line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end
if strfind(line,'new')
    ftype = 'new';
end
ny  =   n(1);   nz  =   n(2);

nzAir = 6; % even number larger than one, resulted by interpolation

y   =   fscanf(fid,'%f',ny);
z   =   fscanf(fid,'%f',nz);

if ~Aniso
    k = 1;
    rho(1:ny,1:nz) = 0;
    while k <= nz
        for j = 1:ny
                rho(j,k) = fscanf(fid,'%f',1);
        end
        k = k+1;
    end
else
    rho(1:ny,1:nz,6) = 0;
    for m=1:6
        k = 1;
        while k <= nz
            for j = 1:ny
                    rho(j,k,m) = fscanf(fid,'%f',1);
            end
            k = k+1;
        end
    end
end
if ftype == 'new'
    k=1;
    while k <= nz
        for j = 1:ny
                par.sig1(j,k) = fscanf(fid,'%f',1);
        end
        k = k+1;
    end    
    k=1;
    while k <= nz
        for j = 1:ny
                par.sig2(j,k) = fscanf(fid,'%f',1);
        end
        k = k+1;
    end 
    k=1;
    while k <= nz
        for j = 1:ny
                par.sig3(j,k) = fscanf(fid,'%f',1);
        end
        k = k+1;
    end    
    k=1;
    while k <= nz
        for j = 1:ny
                par.strike(j,k) = fscanf(fid,'%i',1);
        end
        k = k+1;
    end  
    k=1;
    while k <= nz
        for j = 1:ny
                par.dip(j,k) = fscanf(fid,'%i',1);
        end
        k = k+1;
    end    
    k=1;
    while k <= nz
        for j = 1:ny
                par.slant(j,k) = fscanf(fid,'%i',1);
        end
        k = k+1;
    end        
end
fclose(fid);
