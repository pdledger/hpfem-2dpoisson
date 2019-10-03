function probdata=problem0()


% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=0;
probdata.order=order;

% set mesh spacing
h=0.1;
probdata.h=h;

% Define problem geometry for mesh generator

% set the coordinates of boundary nodes
node=[0 0; 1 0; 1 1; 0 1];

% set the connectivity of the boundary nodes
cnect=[ 1 2; 2 3; 3 4; 4 1];

probdata.node=node;
probdata.cnect=cnect;
probdata.face=[];

% use mesh2d and not meshfaces
probdata.meshfaces=0;


% persibe function for mesh refinement
arg={h};
probdata.meshfun=@meshref
probdata.meshfunarg=arg;


% define boundary dat

% bctype
% bctype 1 Neumann type
% bctype 2 Dirichlet

% bcflag
% The flag to denote the individual boundary segments
% each boundary segment must have a flag and a type!

bcflags=[ 2; 2; 3; 2];
bctype= [0;  2; 2];
% ie each boundary is Dirichlet type and there are different Dirichlet functions

probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
probdata.dirfun=@problem0dir
probdata.dirfunarg=[];

% define the function handle for Neumann BC's
probdata.neufun=[]
probdata.neufunarg=[];

% define the function handle for the source term
probdata.srcfun=@problem0src
probdata.srcfunarg=[];

% define the materials
probdata.mat=1;

% exact solution
probdata.exact=1;
probdata.presscoeff=0;
marg.a=1;
marg.b=1;
probdata.exactfun=@problem0exact;
probdata.exactfunarg=marg;

probdata.dexactfun=@dproblem0exact;
probdata.dexactfunarg=marg;

% straight geometry
lingeom=1;
rin=0; % not used
bcscatter=0; % not used
probdata.lingeom=lingeom;
probdata.bcscatter=bcscatter;



function bcval=problem0dir(x,y,index,arg)
if index==2
bcval=0;
elseif index==3
bcval=1;
else
disp('Not defined for this index')
end

function src=problem0src(x,y,arg)
src=0;

function mat=problem0exact(X,Y,arg)
v=1;
a=arg.a;
b=arg.b;
mat=zeros(size(X));
for k=1:2:20
mat=mat+(4*v/pi)*(sin(k*pi*X/b).*sinh(k*pi*Y/b))./(k*sinh(k*pi*a/b));
end

function dmat=dproblem0exact(X,Y,arg)
a=arg.a;
b=arg.b;
v=1;
dmat=zeros(2,1);
for k=1:2:20
dmat(1)=dmat(1)+(4*v/pi)*(k*pi/b*cos(k*pi*X/b).*sinh(k*pi*Y/b))./(k*sinh(k*pi*a/b));
dmat(2)=dmat(2)+(4*v/pi)*(sin(k*pi*X/b).*k*pi/b*cosh(k*pi*Y/b))./(k*sinh(k*pi*a/b));
end


function h = meshref(x,y,arg)
hglob=arg;

% User defined size function for uniform spacing
r=sqrt(x.^2+y.^2);
h=hglob*sqrt(10*r);
[m n]=size(r);
minh=1e-4*hglob*ones(m,n);
h=max(h,minh);

h=hglob*ones(size(x));

