function probdata=problem5()

% Flow past a symmetric airfoil using the velocity potential

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=3;
probdata.order=order;

% set mesh spacing
h=0.5;
probdata.h=h;

%set grading factor GF=1,2,4 etc
GF=5;

% Define problem geometry for mesh generator

% size of outer box
ymax=3;

%create airfoil geometry
t=0.44
np=20;
dx=1/np;
node=[];
for i=0:np
x=i*dx;
y= (t / 0.20) * (0.29690 * sqrt(x) - 0.12600 *x - 0.35160 * x^2 + 0.28430 * x^3 - 0.10150 * x^4);
node=[node; x-0.5 y];
end
%np+1 points
for i=np-1:-1:1
x=i*dx;
y= (t / 0.20) * (0.29690 * sqrt(x) - 0.12600 *x - 0.35160 * x^2 + 0.28430 * x^3 - 0.10150 * x^4);
node=[node; x-0.5 -y];
end
npp=20;


%include outer (square) surface
size(node)
node   = [node;  %2*np
         ymax  0;  %2*np+1
         ymax  ymax;  %2*np+2
         -ymax  ymax;  %2*np+3
         -ymax  -ymax;  %2*np+4
          ymax  -ymax];  %2*np+5
cnect  = [ 1:(2*np); 2:(2*np) 1]';
[nbodybc dum]=size(cnect);
cnect  = [cnect; 2*np+1 2*np+2];
cnect  = [cnect; 2*np+2 2*np+3; 2*np+3 2*np+4; 2*np+4 2*np+5; 2*np+5 2*np+1];  
[ntotalbc dum]=size(cnect);

probdata.node=node;
probdata.cnect=cnect;
probdata.face=[];


% use mesh2d and not meshfaces
probdata.meshfaces=0;

% persibe function for mesh refinement
arg={h,GF};
probdata.meshfun=@meshref5
probdata.meshfunarg=arg;


% straight geometry
lingeom=1;
probdata.lingeom=lingeom;

% define boundary dat

% bctype
% bctype 1 Neumann type
% bctype 2 Dirichlet

% bcflag
% The flag to denote the individual boundary segments
% each boundary segment must have a flag and a type!

% for Neumann on body
for i=1:nbodybc
bcflags(i)=1;
end
bctype(1)=1; % Neumann

% for Dirichlet farfield 
for i=nbodybc+1:ntotalbc
bcflags(i)=2;
end
bctype(2)=2; % Dirichlet

probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
attachang=0;
vinf=[cosd(attachang) sind(attachang)];
probdata.vinf=vinf;
arg=[];
arg.vinf=vinf;
probdata.dirfun=@problem5dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
%specify the direction of vinf
probdata.neufun=@problem5neu
probdata.neufunarg=arg;

% define the function handle for the source term
probdata.srcfun=@problem5src
probdata.srcfunarg=[];

% define the materials
probdata.mat=1;

% exact solution
arg.vinf=vinf;
probdata.exactfun=@problem5exact;
probdata.exactfunarg=arg;

probdata.dexactfun=@dproblem5exact;
probdata.dexactfunarg=arg;

% tag nodes
tag=[];
probdata.tag=tag;

% tag edges
tagedge=[];
probdata.tagedge=tagedge;	   

% control outputs
probdata.exact=0;
probdata.presscoeff=1;
probdata.presscoefbc=1

% pressure coeefficent output 
probdata.presscoefout=2;

function bcval=problem5dir(x,y,index,arg)
vinf=arg.vinf;
bcval=(x*vinf(1)+y*vinf(2));

function neuval=problem5neu(xp,yp,index,nm,arg)
gradu=[0; 0];
neuval=nm(1)*gradu(1)+nm(2)*gradu(2); 
 
function src=problem5src(x,y,arg)
src=0;

function val=problem5exact(xp,yp,arg)
val=0;


function dvex=dproblem5exact(xp,yp,arg)
dvex=[0; 0];

function h = meshref5(x,y,hglob,GF)

% User defined size function for uniform spacing
h=hglob*ones(size(x));

% refinement towards leading edge
lepx=-0.5*ones(size(x));
lepy=0*ones(size(x));
rg=sqrt((x-lepx).^2+(y-lepy).^2)
[m n]=size(rg);
rmin=2;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
hle=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin)

% refinement towards trailing edge
trpx=0.5*ones(size(x));
trpy=0*ones(size(x));
rg=sqrt((x-trpx).^2+(y-trpy).^2)
[m n]=size(rg);
rmin=2;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
htr=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin)

%pause
h=min(hle,min(htr,hgmax));
