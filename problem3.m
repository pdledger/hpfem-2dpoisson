function probdata=problem3()

% Smooth problem

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=2;
probdata.order=order;

% set mesh spacing
h=0.3;
probdata.h=h;

% Define problem geometry for mesh generator

% set the coordinates of boundary nodes
node=[0 0; pi 0; pi pi; 0 pi];

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

% for Neumann
%bcflags=[ 1; 1; 1; 1];
%bctype= [1];
% for Dirichlet
bcflags=[ 1; 1; 1; 1];
bctype= [2];


% ie each boundary is Neumann type and there is only one index

probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
arg=[];
probdata.dirfun=@problem3dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
probdata.neufun=@problem3neu
probdata.neufunarg=[];

% define the function handle for the source term
probdata.srcfun=@problem3src
probdata.srcfunarg=[];

% define the materials
probdata.mat=1;

% exact solution
probdata.exact=1;
probdata.presscoeff=0;
probdata.exactfun=@problem3exact;
probdata.exactfunarg=[];

probdata.dexactfun=@dproblem3exact;
probdata.dexactfunarg=[];

% tag nodes
tag=[];
probdata.tag=tag;

% tag edges
tagedge=[];
probdata.tagedge=tagedge;	   

%straight geometry
lingeom=1;
rin=0; % not used
bcscatter=0; % not used
probdata.lingeom=lingeom;
probdata.bcscatter=bcscatter;


function bcval=problem3dir(x,y,index,arg)
bcval=0;

function neuval=problem3neu(xp,yp,index,nm,arg)
gradu=[cos(xp)*sin(yp); sin(xp)*cos(yp)];
neuval=nm(1)*gradu(1)+nm(2)*gradu(2); 
 
function src=problem3src(x,y,arg)
src=2*sin(x)*sin(y);

function vex=problem3exact(xp,yp,arg)
vex=sin(xp)*sin(yp);

function dvex=dproblem3exact(xp,yp,arg)
dvex=[cos(xp)*sin(yp); sin(xp)*cos(yp)];

function h = meshref(x,y,arg)
hglob=arg;

% User defined size function for uniform spacing
h=hglob*ones(size(x));

