function probdata=problem1()

% L Shape domain

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=3;
probdata.order=order;

% set mesh spacing
h=0.5;
probdata.h=h;

% set grading factor for mesh (eg 1,2,4,8,16,64)
GF=16

% Define problem geometry for mesh generator

% set the coordinates of boundary nodes
node=[0 0; 0 -1; 1 -1; 1 1; -1 1; -1 0];

% set the connectivity of the boundary nodes
cnect=[ 1 2; 2 3; 3 4; 4 5; 5 6; 6 1];

probdata.node=node;
probdata.cnect=cnect;
probdata.face=[];

% use mesh2d and not meshfaces
probdata.meshfaces=0;

% persibe function for mesh refinement
arg={h,GF};
probdata.meshfun=@meshref;
probdata.meshfunarg=arg;

% define boundary dat

% bctype
% bctype 1 Neumann type
% bctype 2 Dirichlet

% bcflag
% The flag to denote the individual boundary segments
% each boundary segment must have a flag and a type!

bcflags=[ 2; 2; 2; 2; 2; 2];
bctype= [0;  2];
% ie each boundary is Dirichlet type and there are different Dirichlet functions

probdata.bcflags=bcflags;
probdata.bctype=bctype;

% define the function handle for Dirichlet BC's
arg=[];
probdata.dirfun=@problem1dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
probdata.neufun=[]
probdata.neufunarg=[];

% define the function handle for the source term
probdata.srcfun=@problem1src
probdata.srcfunarg=[];

% define the materials
probdata.mat=1;

% exact solution
probdata.exact=1;
probdata.presscoeff=0;
probdata.exactfun=@problem1exact;
probdata.exactfunarg=[];

probdata.dexactfun=@dproblem1exact;
probdata.dexactfunarg=[];


probdata.tag=[0 0; -1 0; 0 -1; 1 -1; 1 1; -1 1];
probdata.tag=[0 0];

lingeom=1;
rin=0; % not used
bcscatter=0; % not used
probdata.lingeom=lingeom;
probdata.bcscatter=bcscatter;

% tag edges
tagedge=[];
probdata.tagedge=tagedge;	   


function bcval=problem1dir(x,y,index,arg)
r=sqrt(x^2+y^2);
ang=atan2(y,x);
if abs(atan2(y,x)-pi) < 1e-10 
 ang=pi;
 elseif abs(atan2(y,x)+pi) < 1e-10
 ang=pi;
end
thvertaxis=pi/2-ang;
bcval=r^(2/3)*sin(2*thvertaxis/3+pi/3);


function src=problem1src(x,y,arg)
src=0;

function vex=problem1exact(xp,yp,arg)
rp=sqrt(xp^2+yp^2);
ang=atan2(yp,xp);
if abs(atan2(yp,xp)-pi) < 1e-10 
 ang=pi;
 elseif abs(atan2(yp,xp)+pi) < 1e-10
 ang=pi;
 end
thvertaxis=pi/2-ang;
vex=rp^(2/3)*sin(2*thvertaxis/3+pi/3);



function dvex=dproblem1exact(xp,yp,arg)
rp=sqrt(xp^2+yp^2);
ang=atan2(yp,xp);
if abs(atan2(yp,xp)-pi) < 1e-10 
 ang=pi;
 elseif abs(atan2(yp,xp)+pi) < 1e-10
 ang=pi;
 end
thvertaxis=pi/2-ang;
vex=rp^(2/3)*sin(2*thvertaxis/3+pi/3);
dvdr=(2/3)*rp^(2/3-1)*sin(2*thvertaxis/3+pi/3);
dvdtheta=rp^(2/3)*cos(2*thvertaxis/3+pi/3)*(-2/3);
      
drdx=0.5/rp*2*xp;
drdy=0.5/rp*2*yp;
      
dthetadx=-xp^2/(xp^2+yp^2)*yp/xp^2;
dthetady=xp^2/(xp^2+yp^2)/xp;
      
dvex=[dvdr*drdx+dvdtheta*dthetadx; dvdr*drdy+dvdtheta*dthetady];
      
function h = meshref(x,y,hglob,GF)

% User defined size function for uniform spacing
r=sqrt(x.^2+y.^2);
h=hglob*sqrt(10*r);
[m n]=size(r);
minh=1e-4*hglob*ones(m,n);
%h=max(h,minh);

%h=hglob*ones(size(x));

% refinement towards centre
rg=sqrt(x.^2+y.^2);
[m n]=size(rg);
rmin=0.5;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
hc=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin);

h=min(hc,hgmax);
