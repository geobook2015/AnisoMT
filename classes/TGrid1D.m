classdef TGrid1D < handle
    %  just makes it cleaner to set the grid up as a class -- so simple not
    %  much more than a structure
    
    properties
        
        Nz   %number of cells in z direction ... including air layers
        Nza   %     number of air layers
        Dz    %     cell dimensions
        NzEarth
    end
    
    methods
        %*******************************************************************
        function obj = TGrid1D(Dz,Nza)
            if nargin ==2
                %   class constructor ... simple
                obj.Nz = length(Dz);
                obj.Nza = Nza;
                obj.NzEarth = obj.Nz-obj.Nza;
                obj.Dz = Dz;
            end
        end
        %*******************************************************************
        function setFrom2Dor3D(obj,grid2Dor3D)
            obj.Nz = grid2Dor3D.Nz;
            obj.Nza = grid2Dor3D.Nza;
            obj.Dz = grid2Dor3D.Dz;
            obj.NzEarth = grid2Dor3D.NzEarth;
        end
    end
end