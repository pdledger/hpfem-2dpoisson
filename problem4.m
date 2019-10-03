function probdata=problem4()

% Flow past a circular cylinder using the velocity potential

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=2;
probdata.order=order;

% set mesh spacing
h=0.5;
probdata.h=h;

% choose either exact or approx farfield bc
%bcchoice=1; % approx bc
 bcchoice=2; % exact bc

% lingeom = 0 Curved
% lingeom = 1 Linear
lingeom=0; 
rin=1;
probdata.rin=rin;
probdata.lingeom=lingeom;
 
% Define problem geometry for mesh generator

% size of outer box
ymax=4;

%create inner circle of points
np=5;
a=rin; 
dphi=pi/np
theta  = (-dphi:-dphi:-2*pi)';
node   = [a*cos(theta) a*sin(theta)];

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

% for Neumann on body
for i=1:nbodybc
bcflags(i)=1;
end
bctype(1)=1; % Neumann
probdata.bcscatter=1;

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
arg.bcchoice=bcchoice;
arg.a=a;

probdata.dirfun=@problem4dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
%specify the direction of vinf
probdata.neufun=@problem4neu
probdata.neufunarg=arg;

% define the function handle for the source term
probdata.srcfun=@problem4src
probdata.srcfunarg=[];

% define the materials
probdata.mat=1;

% exact solution
arg.vinf=vinf;
arg.a=a;
probdata.exactfun=@problem4exact;
probdata.exactfunarg=arg;

probdata.dexactfun=@dproblem4exact;
probdata.dexactfunarg=arg;


% exact solution "on"
probdata.exact=1;
probdata.presscoeff=1;

% tag nodes
tag=[];
probdata.tag=tag;

% tag edges
tagedge=[];
probdata.tagedge=tagedge;	   

% control outputs
probdata.exact=1;
probdata.presscoef=1;
probdata.presscoefbc=1

% pressure coeefficent output 
probdata.presscoefout=1;


function bcval=problem4dir(x,y,index,arg)
vinf=arg.vinf;
bcchoice=arg.bcchoice;
a=arg.a;
if bcchoice==1
% for approx bc
bcval=x*vinf(1)+y*vinf(2);
else
% for exact bc
mvinf=norm(vinf);
r=sqrt(x^2+y^2);
bcval=mvinf*x*(1+(a^2/r^2));
end


function neuval=problem4neu(xp,yp,index,nm,arg)
gradu=[0; 0];
neuval=nm(1)*gradu(1)+nm(2)*gradu(2); 
 
function src=problem4src(x,y,arg)
src=0;

function val=problem4exact(xp,yp,arg)
a=arg.a;
vinf=arg.vinf;
mvinf=norm(vinf);
r=sqrt(xp^2+yp^2);
val=mvinf*xp*(1+(a^2/r^2));


function dvex=dproblem4exact(xp,yp,arg)
a=arg.a;
vinf=arg.vinf;
mvinf=norm(vinf);
r=sqrt(xp^2+yp^2);
vx=((1+(a^2/(r^2)))+xp*((-2*a^2*xp)/(r^4)) )*mvinf;
vy=mvinf*xp*a^2*(-2*yp)/r^4;
dvex=[vx; vy];

function h = meshref(x,y,arg)
hglob=arg;

% User defined size function for uniform spacing
h=hglob*ones(size(x));

