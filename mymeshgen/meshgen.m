function [nelem,nelemq,Mesht,Meshq,bsido,nboun,bcvals,bctype]=meshgen(Mesht,hglob,problem,probdata)

% An interface to the 2D unstructured meshing program written by Darren Engwirda
% outputs: 
% nelem: Number of triangles
% nelemq: Number of quadrilaterials (0)
% Mesht: Structure containing coordinates and connectivities for triangular mesh
% Meshq: as above for quadrilaterials
% bsido: Boundary condition data
% nboun: Number of boundary conditions
% bcvals: Boundary condition values

node=probdata.node;
cnect=probdata.cnect;
bcflags=probdata.bcflags;
bctype=probdata.bctype;
face=probdata.face;

bcvals =[ 0; 0; 1];

% decide which mesh generator to use
mesh2dmeshfaces=probdata.meshfaces;
% hdata = [];
% hdata.fun= probdata.meshfun;
% hdata.args =probdata.meshfunarg;
% %hdata.fun = @h1;
% %hdata.args = {hglob};
% options.dhmax = hglob;

hdata.fun = probdata.meshfun;
hdata.args = probdata.meshfunarg;
opts=[];
opts.dhmax = hglob;
% if mesh2dmeshfaces==0
% [p,t]=mesh2d(node,cnect,hdata,options);
% [nelem dum]=size(t);
% fnum=ones(nelem,1);
% else
% % do not modfiy this code
% [p,t,fnum,stats] = meshfaces(node,cnect,face,hdata,options);
% end

if mesh2dmeshfaces==0
%    [p,t]=mesh2d(node,cnect,hdata,options);
    opts.kind = 'delaunay' ;
    opts.disp = inf;
    [p,etri,t,fnum] = refine2(node,cnect,[],opts,hglob) ;
    [p,etri,t,fnum] = smooth2(p,etri,t,fnum, opts) ;  
    [p,t] = removeOuterbox(p,t);
    [nelem dum]=size(t);
    fnum=ones(nelem,1);
else
% do not modfiy this code
%[p,t,fnum,stats] = meshfaces(node,cnect,face,hdata,options);
    opts.kind = 'delaunay' ;
    opts.disp = inf;
    [p,etri,t,fnum] = refine2(node,cnect,face,opts,hglob) ;
    [p,etri,t,fnum] = smooth2(p,etri,t,fnum, opts) ;
    [p,t] = removeOuterbox(p,t);
end

%create the list of boundaries
% do not modifiy this code
[bsido,p]=my_data(p,t,cnect,bcflags,node);

[nboun dum]=size(bsido);

Mesht.Coordinates=p;
Mesht.Elements=t;
Mesht.Fnum=fnum;
Meshq.Coordinates=p;
Meshq.Elements=[];
Meshq.Fnum=[];
nelemq=0;
[nelem dum]=size(t);
plot_Mesh(Mesht);
figure

function h = h1(x,y,hglob)

% User defined size function for uniform spacing
r=sqrt(x.^2+y.^2);
h=hglob*sqrt(10*r);
[m n]=size(r);
minh=1e-4*hglob*ones(m,n);
h=max(h,minh);

h=hglob*ones(size(x));

% L shape
rg=sqrt(x.^2+y.^2);
[m n]=size(rg);
rmin=0.5;
rgmax=rmin*ones(m,n);
%%Grading Factor
GF=64;
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
h=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin);
h=min(h,hgmax)


%h=min(hglob*sqrt(dx),hglob*sqrt(dy));
%[m n]=size(dx);
%minh=0.1*hglob*ones(m,n);
%h=max(h,minh);

h=hglob*ones(size(x));
