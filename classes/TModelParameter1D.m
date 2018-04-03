classdef TModelParameter1D 
    % class for storing 1D MT model parameters: not a subclass of the
    % abstract TModelParameter class ... at this point!
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
    
    properties
        v           %   1D  array of conductive cells
        paramType   %   linear or log
        AirCond=log(1e-10)     %    air conductivity
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
        function obj = TModelParameter1D(GRID,paramType,v,AirCond)
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
                    obj.v = zeros(obj.grid.NzEarth,1);
                end
                if nargin ==4
                    obj.AirCond = AirCond;
                end
            end
        end
        %******************************************************************
            function obj = setFrom2Dor3D(obj,m)
            obj.ParamGrid  = TGrid1D;
            obj.ParamGrid.setFrom2Dor3D(m.ParamGrid);
            obj.AirCond = m.AirCond;
            obj.paramType = m.paramType;
            switch ndims(m.v)
                case 2
                    [~,n2] = size(m.v);
                    if (n2>1)
                        temp = squeeze(m.v(1,:)).';
                    else
                        temp = m.v;
                    end
                    obj.v = temp;
                case 3
                    obj.v = squeeze(m.v(1,1,:));
                otherwise
                    error('input model parameter must be of dimension 2 or 3')
            end
            end
    end
end
