function probdata=problem2()

% Kellog problem

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order
order=0;
probdata.order=order;

% set mesh spacing
h=0.35;
probdata.h=h;

% set grading factor for mesh (eg 1,2,4,8,16,64)
GF=2

% Define problem geometry for mesh generator

% set the coordinates of boundary nodes
node=[0 0 ; 1 0;  1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1];

% set the connectivity of the boundary nodes
cnect=[ 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 2; 1 2; 1 4; 1 6; 1 8];

% Set the subdomain information
face{1}=[ 9 1 2 10];
face{2}=[ 10 3 4 11];
face{3}=[ 12 11 5 6];
face{4}=[  9 8 7 12];

probdata.node=node;
probdata.cnect=cnect;
probdata.face=face;

% use meshfaces
probdata.meshfaces=1;


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

bcflags=[ 2; 2; 2; 2; 2;2 ;2 ; 2; 0; 0; 0 ; 0];
bctype= [2 2];

% ie each boundary is Neumann type and there is only one index
probdata.bcflags=bcflags;
probdata.bctype=bctype;


% define material
% 0 < gamma  < 2
gamma=0.5%0.5%6%5%0.5 %limit case of 0.1 
[R,sigma,rho]=problem2param(gamma)      

mat=[R; 1 ; R; 1];
%mat=[161.4476387975881e-0; 1 ; 161.4476387975881e0; 1];
probdata.mat=mat;


% define the function handle for Dirichlet BC's
arg=[];
arg.R=R;
arg.sigma=sigma;
arg.rho=rho;
arg.gamma=gamma;
probdata.dirfun=@problem2dir
probdata.dirfunarg=arg;

% define the function handle for Neumann BC's
probdata.neufun=[]
probdata.neufunarg=[];

% define the function handle for the source term
probdata.srcfun=@problem2src
probdata.srcfunarg=[];


% exact solution
probdata.exact=1;
probdata.presscoeff=0;
probdata.exactfun=@problem2exact;
probdata.exactfunarg=arg;

probdata.dexactfun=@dproblem2exact;
probdata.dexactfunarg=arg;

% tag nodes
tag=[0 0];
probdata.tag=tag;

% tag edges
tagedge=[ -1 0 0 0;
           0 -1 0 0;
	   0 0 1 0;
	   0 0 0 1];
%tagedge=[];
probdata.tagedge=tagedge;

% linear geometry
lingeom=1;
rin=0; % not used
bcscatter=0; % not used
probdata.lingeom=lingeom;
probdata.bcscatter=bcscatter;


function bcval=problem2dir(x,y,index,arg)

R=arg.R;
sigma=arg.sigma;
rho=arg.rho;
tau=arg.gamma;

r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
 if abs(atan2(y,x)-pi) < 1e-10 
 theta=pi;
 elseif abs(atan2(y,x)+pi) < 1e-10
 theta=pi;
 end

if theta < 0
theta=theta+2*pi;
end
if theta <= pi/2
mu=cos((pi/2-sigma)*tau)*cos((theta-pi/2+rho)*tau);
elseif theta <= pi
		mu=cos(rho*tau)*cos((theta-pi+sigma)*tau);
elseif theta <= 3*pi/2
			mu=cos(sigma*tau)*cos((theta-pi-rho)*tau);
else
			mu=cos((pi/2-rho)*tau)*cos((theta-3*pi/2-sigma)*tau);
end
bcval=(r^tau)*mu;

 
function src=problem2src(x,y,arg)
src=0;

function vex=problem2exact(x,y,arg)

R=arg.R;
sigma=arg.sigma;
rho=arg.rho;
tau=arg.gamma;


r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
 if abs(atan2(y,x)-pi) < 1e-10 
 theta=pi;
 elseif abs(atan2(y,x)+pi) < 1e-10
 theta=pi;
 end

if theta < 0
theta=theta+2*pi;
end
if theta <= pi/2
mu=cos((pi/2-sigma)*tau)*cos((theta-pi/2+rho)*tau);
elseif theta <= pi
		mu=cos(rho*tau)*cos((theta-pi+sigma)*tau);
elseif theta <= 3*pi/2
			mu=cos(sigma*tau)*cos((theta-pi-rho)*tau);
else
			mu=cos((pi/2-rho)*tau)*cos((theta-3*pi/2-sigma)*tau);
end
vex=(r^tau)*mu;

function dvex=dproblem2exact(x,y,arg)

R=arg.R;
sigma=arg.sigma;
rho=arg.rho;
tau=arg.gamma;


r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
 if abs(atan2(y,x)-pi) < 1e-10 
 theta=pi;
 elseif abs(atan2(y,x)+pi) < 1e-10
 theta=pi;
 end

if theta < 0
theta=theta+2*pi;
end
if theta <= pi/2
mu=cos((pi/2-sigma)*tau)*cos((theta-pi/2+rho)*tau);
dmudtheta=-cos((pi/2-sigma)*tau)*sin((theta-pi/2+rho)*tau)*tau;
elseif theta <= pi
		mu=cos(rho*tau)*cos((theta-pi+sigma)*tau);
		dmudtheta=-cos(rho*tau)*sin((theta-pi+sigma)*tau)*tau;
elseif theta <= 3*pi/2
			mu=cos(sigma*tau)*cos((theta-pi-rho)*tau);
			dmudtheta=-cos(sigma*tau)*sin((theta-pi-rho)*tau)*tau;
else
			mu=cos((pi/2-rho)*tau)*cos((theta-3*pi/2-sigma)*tau);
		    dmudtheta=-cos((pi/2-rho)*tau)*sin((theta-3*pi/2-sigma)*tau)*tau;
end
vex=(r^tau)*mu;


drdx=0.5/r*2*x;
drdy=0.5/r*2*y;
      
dthetadx=-1/(x^2+y^2)*y;
dthetady=x/(x^2+y^2);

dvex=[tau*r^(tau-1)*mu*drdx+r^tau*dmudtheta*dthetadx; 
      tau*r^(tau-1)*mu*drdy+r^tau*dmudtheta*dthetady];
 
function [R,sigma,rho]=problem2param(gamma)      
% range for rho
minrho=max(0,pi*gamma-pi)/2/gamma
maxrho=min(pi*gamma,pi)/2/gamma

rho=(minrho+maxrho)/2

% range for sigma
minrho=-min(pi,2*pi-pi*gamma)/2/gamma
maxrho=-max(0,pi-pi*gamma)/2/gamma

sigma=(minrho+maxrho)/2

R=-tan((pi/2-sigma)*gamma)*cot(rho*gamma)

xnew=[R,rho,sigma]';
x=xnew*0.1
iter=0;
% iterate using Newton's Method
while norm(xnew-x)> 1e-10 & iter < 100
x=xnew;
iter=iter+1;
%compute
f=[R+tan((pi/2-x(3))*gamma)*cot(x(2)*gamma);
    1/R+tan(x(2)*gamma)*cot(x(3)*gamma);
    R+tan(x(3)*gamma)*cot((pi/2-x(2))*gamma)];
    
J=[1 -tan((pi/2-x(3))*gamma)*gamma*csc(gamma*x(2))^2 sec((pi/2-x(3))*gamma)^2*(-gamma)*cot(x(2)*gamma);
   -1/R^2 sec(x(2)*gamma)^2*gamma*cot(x(3)*gamma)   -csc(x(3)*gamma)*gamma*tan(x(2)*gamma);
   1 gamma*csc((pi/2-x(2))*gamma)^2*tan(x(3)*gamma) sec(x(3)*gamma)^2*gamma*cot((pi/2-x(2))*gamma)];
   
% solve for update
delta=J\(-f);

xnew=x+delta;
end
   
R=xnew(1);
rho=xnew(2);
sigma=xnew(3);


function h = meshref(x,y,hglob,GF)

% User defined size function for uniform spacing
r=sqrt(x.^2+y.^2);
h=hglob*sqrt(10*r);
[m n]=size(r);
minh=1e-4*hglob*ones(m,n);
h=max(h,minh);

h=hglob*ones(size(x));

% Kellog domain with refinement towards centre
rg=sqrt(x.^2+y.^2);
[m n]=size(rg);
rmin=0.5;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
hc=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin);

% Kellog domain with refinement towards x axis
rg=sqrt(x.^2);
[m n]=size(rg);
rmin=0.5;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
hx=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin);

% Kellog domain with refinement towards y axis
rg=sqrt(y.^2);
[m n]=size(rg);
rmin=0.5;
rgmax=rmin*ones(m,n);
%%Grading Factor
hmin=hglob/GF;
hgmax=hglob*ones(m,n);
hy=hmin*ones(m,n)+(hglob-hmin).*sqrt(abs(rg)./rmin);

%h=min(min(min(h,hx),hy),hgmax);
h=min(hc,hgmax);

