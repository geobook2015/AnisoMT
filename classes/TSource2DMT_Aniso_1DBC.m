classdef TSource2DMT_Aniso_1DBC < TSource2DMT_1DBC
% Zeqiu Guo, Gary Egbert 2015
% zeqiu_guo@hotmail.com; egbert@coas.oregonstate.edu
% CUGB, Beijing China; OSU, OR USA
    properties
        cxxL,cxxR,cyyL,cyyR,czzL,czzR
        cxyL,cxyR,cyzL,cyzR,cxzL,cxzR
        ii,ib
    end
    methods
        function obj = TSource2DMT_Aniso_1DBC(grid,m)
            % Subclass constructor functions must explicitly call superclass constructors if the superclass
            % constructors require input arguments. The subclass constructor must specify these arguments
            % in the call to the superclass constructor using the constructor output argument and
            % the returned object must be assigned to the constructor output argument.
            if nargin ==2
                obj.grid = grid;
                obj.E1D = TMT1DModel_Aniso(grid,m);
                
                switch upper(m.paramType)
                    case 'LOGE'
                        obj.sigmaLeft = exp([ones(obj.grid.Nza,1)*m.AirCond; m.v(1,:)']);
                        obj.sigmaRight = exp([ones(obj.grid.Nza,1)*m.AirCond; m.v(end,:)']);
                    case 'LINEAR'
%                         obj.cxxL = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.cxx(1,1);  m.cxx(1,:)'];
%                         obj.cxxR = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.cxx(end,1); m.cxx(end,:)'];
%                         obj.cyyL = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.cyy(1,1); m.cyy(1,:)'];
%                         obj.cyyR = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.cyy(end,1); m.cyy(end,:)'];
%                         obj.czzL = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.czz(1,1); m.czz(1,:)'];
%                         obj.czzR = [ones(obj.grid.Nza-1,1)*exp(m.AirCond); m.czz(end,1); m.czz(end,:)'];
%                         obj.cxyL = [zeros(obj.grid.Nza-1,1); m.cxy(1,1); m.cxy(1,:)'];
%                         obj.cxyR = [zeros(obj.grid.Nza-1,1); m.cxy(end,1); m.cxy(end,:)'];
%                         obj.cxzL = [zeros(obj.grid.Nza-1,1); m.cxz(1,1); m.cxz(1,:)'];
%                         obj.cxzR = [zeros(obj.grid.Nza-1,1); m.cxz(end,1); m.cxz(end,:)'];
%                         obj.cyzL = [zeros(obj.grid.Nza-1,1); m.cyz(1,1); m.cyz(1,:)'];
%                         obj.cyzR = [zeros(obj.grid.Nza-1,1); m.cyz(end,1); m.cyz(end,:)'];                        
                        
%                         obj.cxxL = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.cxx(1,1);  m.cxx(1,:)' ];
%                         obj.cxxR = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.cxx(end,1);  m.cxx(end,:)' ];
%                         obj.cyyL = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.cyy(1,1);  m.cyy(1,:)' ];
%                         obj.cyyR = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.cyy(end,1);  m.cyy(end,:)' ];
%                         obj.czzL = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.czz(1,1);  m.czz(1,:)' ];
%                         obj.czzR = [ ones(obj.grid.Nza-1,1)*m.AirCond; m.czz(end,1);  m.czz(end,:)' ];
%                         obj.cxyL = [ zeros(obj.grid.Nza-1,1); m.cxy(1,1); m.cxy(1,:)' ];
%                         obj.cxyR = [ zeros(obj.grid.Nza-1,1); m.cxy(end,1); m.cxy(end,:)' ];
%                         obj.cxzL = [ zeros(obj.grid.Nza-1,1); m.cxz(1,1); m.cxz(1,:)' ];
%                         obj.cxzR = [ zeros(obj.grid.Nza-1,1); m.cxz(end,1); m.cxz(end,:)' ];
%                         obj.cyzL = [ zeros(obj.grid.Nza-1,1); m.cyz(1,1); m.cyz(1,:)' ];
%                         obj.cyzR = [ zeros(obj.grid.Nza-1,1); m.cyz(end,1); m.cyz(end,:)' ];

                        obj.cxxL = [ ones(obj.grid.Nza,1)*m.AirCond;  m.cxx(1,:)' ];
                        obj.cxxR = [ ones(obj.grid.Nza,1)*m.AirCond;  m.cxx(end,:)' ];
                        obj.cyyL = [ ones(obj.grid.Nza,1)*m.AirCond;  m.cyy(1,:)' ];
                        obj.cyyR = [ ones(obj.grid.Nza,1)*m.AirCond;  m.cyy(end,:)' ];
                        obj.czzL = [ ones(obj.grid.Nza,1)*m.AirCond;  m.czz(1,:)' ];
                        obj.czzR = [ ones(obj.grid.Nza,1)*m.AirCond;  m.czz(end,:)' ];
                        obj.cxyL = [ zeros(obj.grid.Nza,1); m.cxy(1,:)' ];
                        obj.cxyR = [ zeros(obj.grid.Nza,1); m.cxy(end,:)' ];
                        obj.cxzL = [ zeros(obj.grid.Nza,1); m.cxz(1,:)' ];
                        obj.cxzR = [ zeros(obj.grid.Nza,1); m.cxz(end,:)' ];
                        obj.cyzL = [ zeros(obj.grid.Nza,1); m.cyz(1,:)' ];
                        obj.cyzR = [ zeros(obj.grid.Nza,1); m.cyz(end,:)' ];
                end
            end
        end
        function setRHS(obj,style) 
            gr = obj.grid;
            switch style
                case 'ThreeE'
                    obj.B = zeros(gr.NNodes+sum(gr.NEdges),1);                    
                    %    LHS
                    obj.E1D.cxx = obj.cxxL;
                    obj.E1D.cyy = obj.cyyL;
                    obj.E1D.czz = obj.czzL;
                    obj.E1D.cxy = obj.cxyL;
                    obj.E1D.cyz = obj.cyzL;
                    obj.E1D.cxz = obj.cxzL;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    
                    ExL = obj.E1D.F.Ex;
                    indexL = gr.vectorIndex(ones(1,gr.Nz+1),1:gr.Nz+1,'node');
                    obj.B(indexL) = ExL;
                    
                    EyL = obj.E1D.F.Ey;
                    indeyL = gr.NNodes+gr.vectorIndex(ones(1,gr.Nz+1),1:gr.Nz+1,'yedge');
                    obj.B(indeyL) = EyL;
                    
                    EzL = obj.E1D.F.Ez;
                    indezL = gr.NNodes+gr.NEdges(1)+gr.vectorIndex(ones(1,gr.Nz),1:gr.Nz,'zedge');
                    obj.B(indezL) = EzL;
                    
                    %    RHS
                    obj.E1D.cxx = obj.cxxR;
                    obj.E1D.cyy = obj.cyyR;
                    obj.E1D.czz = obj.czzR;
                    obj.E1D.cxy = obj.cxyR;
                    obj.E1D.cyz = obj.cyzR;
                    obj.E1D.cxz = obj.cxzR;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    
                    ExR = obj.E1D.F.Ex;
                    indexR = gr.vectorIndex(ones(1,gr.Nz+1)*(gr.Ny+1),1:gr.Nz+1,'node');
                    obj.B(indexR) = ExR;
                    
                    EyR = obj.E1D.F.Ey;
                    indeyR = gr.NNodes+gr.vectorIndex(ones(1,gr.Nz+1)*gr.Ny,1:gr.Nz+1,'yedge');
                    obj.B(indeyR) = EyR;
                    
                    EzR = obj.E1D.F.Ez;
                    indezR = gr.NNodes+gr.NEdges(1)+gr.vectorIndex(ones(1,gr.Nz)*(gr.Ny+1),1:gr.Nz,'zedge');
                    obj.B(indezR) = EzR;
                    
                    %    BOTTOM, revised to use impedance Bottom BCs
                    
                    % Ex Bottom
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [ExL(end);ExR(end)];
                    temp = interp1(X,Y,xi);
                    indexB = (gr.Ny+1)*gr.Nz+2:gr.NNodes-1;
                    obj.B(indexB) = temp(2:end-1);
                    % Ey Bottom
                    xi  = [0; cumsum(obj.grid.DDy)];
                    X = [0;xi(end)];
                    Y = [ EyL(end); EyR(end)];
                    temp = interp1(X,Y,xi);
                    indeyB = gr.NNodes+gr.Ny*gr.Nz+2:gr.NNodes+gr.NEdges(1)-1;
                    obj.B(indeyB) = temp(3:end-2);
                    % Ez Bottom
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [EzL(end);EzR(end)];
                    temp = interp1(X,Y,xi);
                    indezB = gr.NNodes+gr.NEdges(1)+(gr.Ny+1)*(gr.Nz-1)+2 : ...
                        gr.NNodes+sum(gr.NEdges)-1;
                    obj.B(indezB) = temp(2:end-1);                    

                    %    TOPE
                    
                    % Ex Top
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [ExL(1);ExR(1)];
                    temp = interp1(X,Y,xi);
                    indexT = 2:gr.Ny;
                    obj.B(indexT) = temp(2:end-1);
                    % Ey Top
                    xi  = [0; cumsum(obj.grid.DDy)];
                    X = [0;xi(end)];
                    Y = [ EyL(1); EyR(1)];
                    temp = interp1(X,Y,xi);
                    indeyT = gr.NNodes+2:gr.NNodes+gr.Ny-1;
                    obj.B(indeyT) = temp(3:end-2);
                    % Ez Top
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [EzL(1);EzR(1)];
                    temp = interp1(X,Y,xi);
                    indezT = gr.NNodes+gr.NEdges(1)+2:gr.NNodes+gr.NEdges(1)+gr.Ny;
                    obj.B(indezT) = temp(2:end-1);
                    
                    obj.ib = sort([ indexL indexR indexT indexB indeyL indeyR ...
                        indeyT indeyB indezL indezR indezT indezB ]);
                    AllEle = 1:sum(gr.NEdges)+gr.NNodes;
                    AllEle(obj.ib) = [];
                    obj.ii = AllEle;
                case 'EHSta'
                    obj.B = zeros(gr.NNodes+gr.NCells,1);                    
                    %    LHS
                    obj.E1D.cxx = obj.cxxL;
                    obj.E1D.cyy = obj.cyyL;
                    obj.E1D.czz = obj.czzL;
                    obj.E1D.cxy = obj.cxyL;
                    obj.E1D.cyz = obj.cyzL;
                    obj.E1D.cxz = obj.cxzL;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    ExL = obj.E1D.F.Ex;
                    HxL = obj.E1D.F.Hx;
                    indeL = gr.vectorIndex(ones(1,gr.Nz-1),2:gr.Nz,'node');
                    obj.B(indeL) = ExL(2:end-1);
                    
                    indhL = gr.NNodes+gr.vectorIndex(ones(1,gr.Nz),1:gr.Nz,'cell');
                    obj.B(indhL) = HxL;
                    
                    %    RHS
                    obj.E1D.cxx = obj.cxxR;
                    obj.E1D.cyy = obj.cyyR;
                    obj.E1D.czz = obj.czzR;
                    obj.E1D.cxy = obj.cxyR;
                    obj.E1D.cyz = obj.cyzR;
                    obj.E1D.cxz = obj.cxzR;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    ExR = obj.E1D.F.Ex;
                    HxR = obj.E1D.F.Hx;                    
                    indeR = gr.vectorIndex(ones(1,gr.Nz-1)*(gr.Ny+1),2:gr.Nz,'node');
                    obj.B(indeR) = ExR(2:end-1);
                    
                    indhR = gr.NNodes+gr.vectorIndex(ones(1,gr.Nz)*(gr.Ny),1:gr.Nz,'cell');
                    obj.B(indhR) = HxR;
                    
                    %    BOTTOM
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [ExL(end);ExR(end)];
                    temp = interp1(X,Y,xi);
                    indeB = (gr.Ny+1)*gr.Nz+1:gr.NNodes;
                    obj.B(indeB) = temp;                    
                    
                    xi  = [0; cumsum(obj.grid.DDy)];
                    X = [0;xi(end)];
                    Y = [ HxL(end); HxR(end)];
                    temp = interp1(X,Y,xi);
                    indhB = gr.NNodes+gr.Ny*(gr.Nz-1)+2:gr.NNodes+gr.NCells-1;
                    obj.B(indhB) = temp(3:end-2);                    

                    %    TOPE
                    
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [ExL(1);ExR(1)];
                    temp = interp1(X,Y,xi);
                    indeT = 1:gr.Ny+1;
                    obj.B(indeT) = temp;                    
                    
                    xi  = [0; cumsum(obj.grid.DDy)];
                    X = [0;xi(end)];
                    Y = [ HxL(1); HxR(1)];
                    temp = interp1(X,Y,xi);
                    indhT = gr.NNodes+2:gr.NNodes+gr.Ny-1;
                    obj.B(indhT) = temp(3:end-2);                        
                    
%                     ExT = ones(gr.Ny+1,1)*obj.E1D.F.Ex(1);
%                     indeT = 1:gr.Ny+1;
%                     obj.B(indeT) = ExT;
%                     indhT = gr.NNodes+2:gr.NNodes+gr.Ny-1;
%                     obj.B(indhT) = obj.E1D.F.Hx(1);
                    
                    obj.ib = sort([ indeL indeR indeT indeB indhL indhR indhT indhB ]);
%                     obj.ib = sort([ indeL indeR indeT indeB indhT indhB ]);                    
                    AllEle = 1:gr.NCells+gr.NNodes;
                    AllEle(obj.ib) = [];
                    obj.ii = AllEle;
                case 'EHFix'
                    obj.B = zeros(2*gr.NNodes,1);                   
                    %    LHS
                    obj.E1D.cxx = obj.cxxL;
                    obj.E1D.cyy = obj.cyyL;
                    obj.E1D.czz = obj.czzL;
                    obj.E1D.cxy = obj.cxyL;
                    obj.E1D.cyz = obj.cyzL;
                    obj.E1D.cxz = obj.cxzL;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    ExL = obj.E1D.F.Ex;
                    % another solution is to solve 1D anisotropy 
                    % problem of coupled Ex and Hx fields, then
                    % Ex and Hx would be on nodes 1D profiles
                    HxL = obj.E1D.F.Hx;
                    Hxl = (HxL(1:end-1).*obj.E1D.Dz(1:end-1)+HxL(2:end).*...
                        obj.E1D.Dz(2:end))./obj.E1D.DDz/2;                         
                    Hsl=HxL(1)+(HxL(1)-HxL(2))*obj.E1D.Dz(1)/obj.E1D.Dz(2);
                    HxL = [Hsl; Hxl; HxL(end)];
                    indeL = gr.vectorIndex(ones(1,gr.Nz-1),2:gr.Nz,'node');
                    obj.B(indeL) = ExL(2:end-1);
                    
                    indhL = gr.NNodes+gr.vectorIndex(ones(1,gr.Nz+1),1:gr.Nz+1,'node');
                    obj.B(indhL) = HxL;
                    
                    %    RHS
                    obj.E1D.cxx = obj.cxxR;
                    obj.E1D.cyy = obj.cyyR;
                    obj.E1D.czz = obj.czzR;
                    obj.E1D.cxy = obj.cxyR;
                    obj.E1D.cyz = obj.cyzR;
                    obj.E1D.cxz = obj.cxzR;
                    obj.E1D.setEquations(obj.SourceParams.omega);
                    obj.E1D.Solve(obj.SourceParams.Polarization);
                    ExR = obj.E1D.F.Ex;
                    HxR = obj.E1D.F.Hx;
                    Hxr = (HxR(1:end-1).*obj.E1D.Dz(1:end-1)+HxR(2:end).*...
                        obj.E1D.Dz(2:end))./obj.E1D.DDz/2;                         
                    Hsr=HxR(1)+(HxR(1)-HxR(2))*obj.E1D.Dz(1)/obj.E1D.Dz(2);    
                    HxR = [Hsr; Hxr; HxR(end)];                    
                    indeR = gr.vectorIndex(ones(1,gr.Nz-1)*(gr.Ny+1),2:gr.Nz,'node');
                    obj.B(indeR) = ExR(2:end-1);
                    
                    indhR =  gr.NNodes+gr.vectorIndex(ones(1,gr.Nz+1)*(gr.Ny+1),1:gr.Nz+1,'node');
                    obj.B(indhR) = HxR;
                    
                    %    BOTTOM
                    xi  = [0; cumsum(obj.grid.Dy)];
                    X = [0;xi(end)];
                    Y = [ExL(end);ExR(end)];
                    temp = interp1(X,Y,xi);
                    indeB = (gr.Ny+1)*gr.Nz+1:gr.NNodes;
                    obj.B(indeB) = temp;                 
                    
                    Y = [HxL(end);HxR(end)];
                    temp = interp1(X,Y,xi);
                    indhB = gr.NNodes+(gr.Ny+1)*gr.Nz+2:2*gr.NNodes-1;
                    obj.B(indhB) = temp(2:end-1);                            
                    
                    %    TOP
                    ExT = ones(gr.Ny+1,1)*obj.E1D.F.Ex(1);
                    indeT = 1:gr.Ny+1;
                    obj.B(indeT) = ExT;
                                   
                    indhT = gr.NNodes+2:gr.NNodes+gr.Ny;
%                     xi  = [0; cumsum(obj.grid.Dy)];
%                     X = [0;xi(end)];                    
%                     Y = [HxL(1);HxR(1)];
%                     temp = interp1(X,Y,xi);
                    obj.B(indhT) = (Hsl+Hsr)/2;
                    
                    obj.ib = sort([ indeL indeR indeT indeB indhL indhR indhT indhB ]);
                    AllEle = 1:2*gr.NNodes;
                    AllEle(obj.ib) = [];
                    obj.ii = AllEle;
                otherwise
                    error('The style you specified has not been implemented')
            end
        end
    end
end
