function main

close all
format long e

% this is the main part of a hpfem Poisson Solver without
% error estimation or adaptivity
% from this function other functions are called to perform
% the finite element analysis
% the program offers the possibility for either meshes of
% triangles or meshes of quadrilaterials (or hybrid combinations)

% add this directory and all subfolders
addpath(genpath('./'))

% problem - 0 Trough problem
% problem - 1 L Shape domain
% problem - 2 Kellog problem
% problem - 3 smooth problem
% problem - 4 flow past cylinder, curved geometry
% problem - 5 flow past a 00** airfoil

% select problem type
problem=4;

% Select Mesh type
% Meshtype 1 - trinagles
% Meshtype 2 - quadrilaterials (problems 0 & 3 only)
Meshtype=1;
if Meshtype == 2 & problem~=0 & problem~=3
    Meshtype=1;
end

% The problem files contain the settings (eg source functions, boundary
% conditions, mesh details and order for the problem)

if problem==0

% define all data for this problem
probdata=problem0();

elseif problem==1

% define all data for this problem
probdata=problem1();

elseif problem==2

probdata=problem2();


elseif problem==3

% define all data for this problem
probdata=problem3();

elseif problem==4

% define all data for this problem
probdata=problem4();

elseif problem==5

% define all data for this problem
probdata=problem5();

end

% Set order of elements;
order=probdata.order;
maxorder=order;

% set mesh spacing
h=probdata.h;

% generate an unstructured mesh of trinagles using the mesh2d program
% the domain definition and boundary condition can be modfied in meshgen.m
% comment the next lines to switch off triangle mesh generator.
if Meshtype==1
Mesht=[];

[nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals,bctype]=meshgen(Mesht,h,problem,probdata);

elseif Meshtype==2

% At present 
% Simple elliptic quadrilaterial mesh generator
% This generator can only mesh rectangular domains. The coordinates and
% boundary conditions can be modifed in ellip.m
% un comment the next line to switch on the quadrilaterial mesh generator.
if problem==0 | problem ==3
h=0.2
[Meshq,bsido,nboun,Mesht,nelemt,nelemq,bcvals]=ellip(h,probdata);
elseif problem==1
% apply also 90 degree (pi/2) rotation to generated mesh
    [nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals]=read_Mesh('edge.plt',-pi/2);
else
    disp('problem not defined');
end
end



% Assign constant order to all elements
for i=1:nelemt+nelemq
orderel(i)=order;
end


% precompute all required basis functions for a sufficently high p

% generate integration points on the interval [0,1];
disp(['Computing basis functions at all integration points for order= ',num2str(maxorder)])
disp('...this may take a little time!')

[x,w]=gaussquad(0.,1,2*(maxorder+2)); 

% evaluate the basis functions at the integration points across a reference element
if nelemq~=0
[basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1]=myevalh1(maxorder,x);
[basisxqhdiv,basisyqhdiv,sdivbasisqhdiv,rbasisxqhdiv,rbasisyqhdiv,l2basisq,rbasisql2]=myevalqhdiv(maxorder+1,x,w);


else
basisxh1=[];
basisyh1=[];
sbasish1=[];
rbasish1=[];
divbasisxxh1=[];divbasisxyh1=[];divbasisyxh1=[];divbasisyyh1=[];
basisxqhdiv=[];basisyqhdiv=[];sdivbasisqhdiv=[];rbasisxqhdiv=[];rbasisyqhdiv=[];l2basisq=[];rbasisql2=[];
lbasisxh1=[];lbasisyh1=[];lsbasish1=[];lrbasish1=[];
ldivbasisxxh1=[];ldivbasisxyh1=[];ldivbasisyxh1=[];ldivbasisyyh1=[];
lbasisxqhdiv=[];lbasisyqhdiv=[];lsdivbasisqhdiv=[];
lrbasisxqhdiv=[];lrbasisyqhdiv=[];ll2basisq=[];lrbasisql2=[];
end

% compute the integration weights and locations 
[intxi,inteta,intw]=gautri(2*(maxorder+2));

% NB for edges use different set to those used with quadrilaterials! 
[xt,wt]=gaussquad(-1.,1,2*(maxorder+2)); 

% work out basis at each integration point
if nelemt~= 0
[basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1]=myevalth1(maxorder,intxi,inteta,xt);

[basisxthdiv,basisythdiv,sdivbasisthdiv,rbasisxthdiv,rbasisythdiv,l2basist,rbasistl2,gl2basistx,gl2basisty]=myevalthdiv(maxorder+1,intxi,inteta,xt);
else
basisxth1=[];
basisyth1=[];
sbasisth1=[];
rbasisth1=[];
divbasisxxth1=[];
divbasisxyth1=[];
divbasisyxth1=[];
divbasisyyth1=[];
basisxthdiv=[]; basisythdiv=[];
sdivbasisthdiv=[];
rbasisxthdiv=[];rbasisythdiv=[];
l2basist=[];rbasistl2=[];gl2basistx=[];gl2basisty=[];
lbasisxth1=[];lbasisyth1=[];lsbasisth1=[];
lrbasisth1=[];
ldivbasisxxth1=[];ldivbasisxyth1=[];ldivbasisyxth1=[];ldivbasisyyth1=[];
lbasisxthdiv=[];lbasisythdiv=[];lsdivbasisthdiv=[];
lrbasisxthdiv=[];lrbasisythdiv=[];ll2basist=[];lrbasistl2=[];
end
disp('finished computing all basis functions')

[nelemt dum]=size(Mesht.Elements);

% extract the correct basis functions
if nelemt~=0
[lbasisxth1,lbasisyth1,lsbasisth1,lrbasisth1,ldivbasisxxth1,ldivbasisxyth1,ldivbasisyxth1,ldivbasisyyth1,lbasisxthdiv,lbasisythdiv,lsdivbasisthdiv,lrbasisxthdiv,lrbasisythdiv,ll2basist,lrbasistl2]=extractbasis(basisxth1,...
basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1,order,maxorder,intw,intxi,inteta,xt,wt...
,basisxthdiv,basisythdiv,sdivbasisthdiv,rbasisxthdiv,rbasisythdiv,l2basist,rbasistl2);
else
[lbasisxh1,lbasisyh1,lsbasish1,lrbasish1,ldivbasisxxh1,ldivbasisxyh1,ldivbasisyxh1,ldivbasisyyh1,lbasisxqhdiv,lbasisyqhdiv,lsdivbasisqhdiv,lrbasisxqhdiv,lrbasisyqhdiv,ll2basisq,lrbasisql2]=extractbasisq(x,w,...
basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1,order,maxorder,...
basisxqhdiv,basisyqhdiv,sdivbasisqhdiv,rbasisxqhdiv,rbasisyqhdiv,l2basisq,rbasisql2);
end

% Determine the numbering of the edges in the mesh
[ne,edges,glob,dir]=edgeno(Mesht,Meshq,bsido,nelemt,nelemq,nboun);



% default values -1 no curved edges
lingeom=probdata.lingeom;
coefft=-1*ones(nelemt,3,2);
coeff=-1*ones(nelemq,4,2);
if lingeom==0
rin=probdata.rin;
bcscatter=probdata.bcscatter;
[coefft,coeff,Mesht,Meshq]=getcoeff(nelemq,nelemt,Mesht,Meshq,edges,glob,ne,coefft,coeff,bcscatter,rin);
end

% Call hp-fem solver and compute solution
plotson=1;
mainhpfem(nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals,ne,edges,glob,dir,...
order,probdata,orderel,maxorder,lbasisxh1,lbasisyh1,lsbasish1,lrbasish1,ldivbasisxxh1,ldivbasisxyh1,...
ldivbasisyxh1,ldivbasisyyh1,x,w,lbasisxth1,lbasisyth1,lsbasisth1,lrbasisth1,ldivbasisxxth1,ldivbasisxyth1,ldivbasisyxth1,ldivbasisyyth1,...
intxi,inteta,intw,xt,wt,plotson,coeff,coefft);

