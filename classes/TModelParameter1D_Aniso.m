classdef TModelParameter1D_Aniso < TModelParameter1D 
    properties
        cxx,cyy,czz,cxy,cyz,cxz
    end
    methods
        function obj = TModelParameter1D_Aniso(GRID,paramType,v,AirCond)
            % class constructor
            if nargin >=1
                obj.grid = GRID;
                obj.ParamGrid = GRID;
                if nargin == 1
                    obj.paramType='LINEAR';
                else
                    obj.paramType = paramType;
                end
                if nargin >=3
                    obj.cxx = v.cxx;
                    obj.cyy = v.cyy;
                    obj.czz = v.czz;
                    obj.cxy = v.cxy;
                    obj.cxz = v.cxz;
                    obj.cyz = v.cyz;                    
                end
                if nargin ==4
                    obj.AirCond = AirCond;
                end
            end
        end
        
        function obj = setFrom2Dor3D(obj,m)
            obj.ParamGrid  = TGrid1D;
            obj.ParamGrid.setFrom2Dor3D(m.ParamGrid);
            obj.AirCond = m.AirCond;
            obj.paramType = m.paramType;
            switch ndims(m.cxx)
                case 2
                    [~,n2] = size(m.cxx);
                    if (n2>1)
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
        end
    end
end

