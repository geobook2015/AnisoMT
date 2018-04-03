classdef TSource < handle
    %   abstract base class to define forcing (sources and boundary
    %   conditions) for induction  equations  (should this really be
    %   handle class?)
    %   it needs to know a full A
    properties
        b        % RHS--sources and BCs (1D array)
        B        % RHS--sources and BCs (TVector)
        e0       % initial solution (TVector)
        grid
        SourceParams    %   structure/object containing defining source
                        %  parameters (e.g., for MT  omega,
                        %  polarization)
        nonZeroSource
    end
  
methods

end
   methods (Abstract)
       %    after creating the object, a routine to set BC/source for 
       %    a specific source
       result = SetSourceParams(obj,info)
       %    then a routine to create rhs vector
       result = setRHS(obj)
       %    and perhaps a routine to set initial solution for iterative
       %    solver
       %result = initSolnVector(obj)
   end       
end