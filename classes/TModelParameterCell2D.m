classdef TModelParameterCell2D < TModelParameter
    % class for storing 2D MT model parameters
    %   tuned specifically for WiFi3D (in the sense that the grid
    %   object is explicitly of type TGrid3D
    %   slighlty changed from Object Oriented inversion
    %    Need to clean up treatment of grid object ... not abstract
    %    enough, because originally this was not treated as an actual
    %    object, just a structure with known field names
    %   NOTE that grid.Nz and grid.Dz include air layers!  These have
    %  to be added in read routine
    
    %  NOTE: for 2D I am simplifying things, since  I don't see supporting 
    %    different grids or equation formulations for 2D.  Thus we assume
    %     the standard 2D FD approach and implement the PDEmapping routine
    %     in this class already.  Also, we allow a mode flag to control TE vs.
    %     TM
    
    properties %(SetAccess = private)
        v           %     array of conductive cells
        paramType   %   linear or log
        AirCond=log(1e-6)     %    air conductivity
        ParamGrid        %  this model parameter does depend on a TGrid2D object,
        %   which in the simplest case will be also the
        %   grid used for solving the system of equations
        %   but this might not be the same grid used for
        %   the model mapping in other cases, so lets give
        %   the property a unique name!
        grid
    end
    
    methods
        %*******************************************************************
        function obj = TModelParameterCell2D(GRID,paramType,v,AirCond)
            % class constructor
            if nargin >=1
                obj.grid = GRID;
                obj.ParamGrid = GRID;
                if nargin == 1
                    obj.paramType='LOGE';
                else
                    obj.paramType = paramType;
                end
                if nargin >=3
                    obj.v = v;
                else
                    obj.v = zeros(obj.grid.Ny,obj.grid.NzEarth);
                end
                if nargin ==4
                    obj.AirCond = AirCond;
                end
            end
        end
        
        %*******************************************************************
        function result =length(obj)
            %  finds length of MT2DmodelParam object
            [ny,nz] = size(obj.v);
            result = ny*nz;
        end
        %*******************************************************************
        function result = maxAbs(obj)
            % maximum model parameter magnitude (absolute value)
            result = max(max(abs(obj.v)));
        end
        %*******************************************************************
        function obj = SetFrom3D(obj,Cond,XY,slice)
            %  fill in the contents of a 2D MT model parameter object from a
            %       3D conductivity model  input structure
            gr = MT2Dgrid;
            MT2DgridFrom3Dgrid(gr,Cond.grid,XY)
            obj.grid = gr;
            obj.paramType = 'LOGE';
            if strcmp(Cond.paramType,'LOG10')
                obj.AirCond = log(10)*Cond.AirCond;
                if strcmpi(XY,'x')
                    obj.v = log(10)*squeeze(Cond.v(slice,:,:));
                else
                    obj.v = log(10)*squeeze(Cond.v(:,slice,:));
                end
            else
                obj.AirCond = log(Cond.AirCond);
                if strcmpi(XY,'x')
                    obj.v = log(squeeze(Cond.v(slice,:,:)));
                else
                    obj.v = log(squeeze(Cond.v(:,slice,:)));
                end
            end
        end
        %******************************************************************
        function result = PDEmapping(obj,mode)
            %   map model parameter to 2D TScalar or TVector object, giving
            %    conductivity on nodes (mode = 'TE') or resistivity on
            %    edges (mode = 'TM')
            
            if nargin == 1
                mode = 'TE';
            end
            %   compute conductivity on all 
            k1  = obj.grid.Nza+1;
            k2 = obj.grid.Nz;
            cellCond = TScalar2D(obj.grid,'cell');
            switch lower(obj.paramType)
                case 'loge'
                    cellCond.v = cellCond.v + exp(obj.AirCond);
                    cellCond.v(:,k1:k2) = exp(obj.v);
                case 'linear'
                    cellCond.v = cellCond.v + exp(obj.AirCond);
                    cellCond.v(:,k1:k2) = obj.v;
                otherwise
                    error(['model parameter mapping not coded for paramtype' obj.paramType ])
            end
            switch mode
                case 'TE'
                    result = cellCond.vecAvg('node');
                case 'TM'
                    cellCond.v = 1./cellCond.v;
                    %   I am keeping resistivity at surface equalt to Earth
                    %   resistivity, not averaging with air
                    cellCond.v(:,k1-1) = cellCond.v(:,k1);
                    result = cellCond.vecAvg('edge');
                otherwise
                    error('mode must be TE or TM')
            end
        end
        %*******************************************************************
        function [obj] = readVec(obj,cfile)
            %   reads from ascii 2D model conductivity file cfile, puts
            %   contents into already created obj
            % Usage:  [m] = readCond(cfile)
            
            [dy,dz,rho,type] = read_mackie2d_model(cfile);
            
            nzAir = 10;     
            dzAir(1:nzAir) = 0;
            dzAir(1) = max(dz(1),10.0);
            for k = 2:nzAir
                dzAir(k) = dzAir(k-1)*3;
            end
            obj.ParamGrid = TGrid2D(dy,[dzAir(end:-1:1)' ; dz],nzAir);
            
            
            if strfind(type,'LOGE')
                obj.v = - rho;
            else
                obj.v = 1./rho;               
            end
            obj.paramType = type;
            
        end
        %*******************************************************************
        function [status] = writeVec(obj,cfile)
            %
            % Usage:  [status] = writeCond_2D(cfile,cond)
            %
            %  writes ModelParam object, provided as structure
            %  status is total number of bytes written
            
            nza = obj.grid.Nza;
            nz  = obj.grid.Nz;
            
            dy = obj.grid.Dy;
            dz = obj.grid.Dz(nza+1:nz);
            
            type = obj.paramType;
            
            if strcmp(type,'LOGE')
                rho = - obj.v;
            else
                rho = 1./(obj.v);
            end
            
            status = write_mackie2d_model(cfile,dy,dz,rho,type);
        end
        %*******************************************************************
        function obj = plus(obj1,obj2)
            obj = obj1;
            obj.v = obj1.v+obj2.v;
        end
        %*******************************************************************
        function obj = minus(obj1,obj2)
            obj = obj1;
            obj.v = obj1.v-obj2.v;
        end
        %*******************************************************************
        function [m] = linComb(c1,m1,c2,m2)
            %
            %   computes linear combination of model parameters
            %                 m = c1*m1+c2*m2
            %   This version is for a single model parameter
            %
            %  Usage: [m] = linCombMod(c1,m1,c2,m2);
            
            m = m1;
            m.v = c1*m1.v + c2*m2.v;
        end
        %*******************************************************************
        function [mOut] = mtimes(c,mIn)
            %  Initial implementation of scalar multiplication of a model
            %   space object by a real scalar; no error type checking;
            %   scalar has to be first argument
            %
            %  Usage: [mOut] = mtimes(c,mIn)
            if ~isa(c,'double')
                error('Error: can only multiply data vectors by doubles')
            end
            
            mOut = mIn;
            mOut.v = c*mOut.v;
        end
        %*******************************************************************
        function [ip] = dot(m1,m2)
            %  dot product for (real) model parameters m1, m2
            %
            %  Usage:  [ip] = dot(m1,m2);
            
            ip = sum(sum(m1.v .* m2.v));
        end
        %*******************************************************************
        function [mOut] = InitHalfSpace(mIn,sig)
            %   Usage: [mOut] = InitHalfSpace(mIn,sig);
            %   Initialize a model parameter mOut as half space
            %     mIn = template for model parameter
            %     sig = conductivity (NOT log conductivity,
            %		or resistivity) for half space
            %
            %  set prior, and starting conductivity (uniform half space)
            mOut = mIn;
            if strcmp(mIn.paramType,'LOGE')
                mOut.v = log(sig)*ones(size(mIn.v));
            else
                mOut.v = sig*ones(size(mIn.v));
            end
        end
        %*******************************************************************
        function [obj] = ZeroVec(obj)
            %   Usage: [mOut] = ZeroVec(obj);
            %    zero existing model parameter
            obj.v = zeros(size(obj.v));
        end
        %*******************************************************************
        function [mOut] = oneD_average(mIn)
            %   Usage: [mOut] = oneD_averae(mIn);
            %    average mIn to a 1-D model
            mOut = mIn;
            mOut.v = ones(mIn.grid.Ny,1)*mean(mIn.v,1);
        end
        %*******************************************************************
        function obj = setModelParam(obj,v,paramType,AirCond,Grid)
            %  sets MT3DmodelParam properties--assume calling routine makes
            %  sure that AirCond and paramType are consistent
            obj.v = v;
            if nargin > 2
                obj.paramType = paramType;
                if nargin > 3
                    obj.AirCond = AirCond;
                    if nargin > 4
                        obj.ParamGrid = Grid;
                    end
                end
            end
        end
        %*******************************************************************
        function mVec = GetVector(obj)
            %   extracts model parameter values and returns as a standard vector
            mVec = obj.v;
            [ny,nz] = size(mVec);
            mVec = reshape(mVec,ny*nz,1);
        end
        %*******************************************************************
        function obj = SetVector(obj,v)
            %   inserts  standard vector v into model parameter obj (as an
            %   array)
            vobj.v = ...
                reshape(v,[obj.grid.Ny,obj.grid.NzEarth]);
        end
        
        %*******************************************************************
        function [obj] = zero(obj)
            %   Usage: [mOut] = InitHalfSpace(mIn,sig);
            %    zero existing model parameter
            obj.v = zeros(size(obj.v));
        end
        
        %***********************************************
        function diff  = relDiff(obj1,obj2)
            temp = obj1-obj2;
            diff = dot(temp,temp)/dot(obj1,obj1);
            diff = sqrt(diff);
        end
        %*******************************************************************
        function [h,hCB] = plotCond(m,OPTIONS)
            % Plots conductivity model
            % Usage [h,hCB] = plotCond(m,grid,OPTIONS)
            %   m is the 2D conductivity model to plot
            %   OPTIONS is a structure of plotting OPTIONS
            %          .nySkip = number of cells to omit from
            %                    each end of the profile
            %          .nZplot = number of vertical layers to plot
            %                    (starting from the earth surface)
            %          .ncax = color axis limits (vector with 2 elements)
            %          .title = plot title
            
            if nargin ==  1
                nZplot = m.grid.Nz-m.grid.Nza;
                nYskip = 5;
                cax = [-3.,0];
                OPTIONS = struct('Contour',0,'nYskip',nYskip,'nZplot',nZplot,...
                    'cax',cax,'title','');
            end
            if ~isfield(OPTIONS,'Contour')
                OPTIONS.Contour = 0;
            end
            if strcmp(m.paramType,'LOGE')
                m.v = log10(exp(m.v));
            else
                m.v = log10(m.v);
            end
            
            figure('Position',[100,100,600,400],...
                'PaperPosition',[1,1,6,4])
            y = cumsum([0; m.grid.Dy]);
            z = cumsum([0; m.grid.Dz]);
            zCenter = z-z(m.grid.Nza+1);
            zCenter = zCenter(m.grid.Nza+1:end)/1000;
            yCenter = y(OPTIONS.nYskip+1:end-OPTIONS.nYskip)/1000;
            yCenter = yCenter-mean(yCenter);
            Cond = m.v;
            Cond = [Cond Cond(:,end)];
            Cond = [Cond ;Cond(end,:)];
            Cond = Cond(OPTIONS.nYskip+1:end-OPTIONS.nYskip,1:OPTIONS.nZplot);
            if OPTIONS.Contour
                %   positive contours: solid line
                cL  = OPTIONS.cLev;
                cL = cL(cL>0);
                [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                    cL,'Linewidth',1,'color','b');
                set(gca,'FontWeight','bold','FontSize',16);
                clabel(h,hCB,cL(2:2:end),'FontSize',14,'fontweight','bold');
                hold on
                cL  = OPTIONS.cLev;
                cL = cL(cL<0);
                [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                    cL,'Linewidth',1,'color','r','linestyle','--');
                clabel(h,hCB,cL(2:2:end),'FontSize',14,'fontweight','bold');
                [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                    [0,0],'Linewidth',2,'Color',[0,0,0]);
                if isempty(h)
                    clabel(h,hCB,'manual','FontSize',14,'fontweight','bold')
                end
                axis('ij');
            else
                h = pcolor(yCenter,zCenter(1:OPTIONS.nZplot),Cond'); ...
                    shading flat;
                axis('ij');
                set(gca,'FontWeight','demi','FontSize',13);
                caxis(OPTIONS.cax);
                c = colmap;
                c = c(end:-1:1,:);
                X = 1:17;
                XI = 1:.25:17;
                c = interp1(X,c,XI);
                colormap(c)
                hCB = colorbar;
                yt = floor(OPTIONS.cax(1)):1:ceil(OPTIONS.cax(2));
                ytLabel = 10.^(-yt);
                
                set(hCB,'FontWeight','demi','FontSize',12,...
                    'Ytick',yt,'YtickLabel',ytLabel);
            end
            ylabel('Depth (km)');
            xlabel('km');
            title(OPTIONS.title);
        end
    end    % methods
end   %   classdef
