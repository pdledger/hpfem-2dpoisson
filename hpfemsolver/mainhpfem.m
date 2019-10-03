function mainhpfem(nelemt,nelemq,Mesht,Meshq,bsido,nboun,bcvals,ne,edges,glob,dir,...
order,probdata,orderel,maxorder,basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,...
divbasisyxh1,divbasisyyh1,x,w,basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1,...
intxi,inteta,intw,xt,wt,plotson,coeff,coefft)

% set the volume charge density (assumed constant throughout the domain)
rho=0;

%order=maxorder;

% Determine the number of unknowns
bctype=probdata.bctype;
[nunk,nunki,unkv,unkvh,unkvi,unkqih,ndir,helpbc,orderel]=nounk(ne,edges,nelemq,nelemt,Mesht,nboun,bsido,glob,Meshq,order,bctype,orderel);

% Assembly process (part 1)
disp('assemlbling system')

map=[];
irn=[];
icn=[];
nz=0;
% quadrilaterials
[rhs,I,J,X,dirval,nz]=assembleqh1(nelemq,Meshq,...        
order,glob,nunk,dir,nunki,edges,unkv,unkvh,unkvi,ndir,unkqih,x,w,helpbc,bcvals,rho,probdata,basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,...
divbasisyxh1,divbasisyyh1,coeff);

%triangles
[rhs,stiff,dirval]=assembleth1(nelemt,Mesht,...        
order,nz,glob,nunk,dir,nunki,edges,nelemq,unkv,unkvh,unkvi,ndir,helpbc,I,J,X,bcvals,rhs,dirval,rho,probdata,basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1,...
intxi,inteta,intw,xt,wt,coefft);

disp('assembly complete, solving system')

%sol=stiff\rhs;
sol=pcg(stiff,rhs,1e-10,nunki);

disp('Determing interior degrees of freedom')
if order > 0
% Assemble process (part 2)
% triangles
[sol]=assembletint(nelemt,Mesht,order,glob,nunk,dir,nunki,edges,nelemq,unkv,...
unkvh,unkvi,intxi,inteta,intw,xt,wt,sbasisth1,rbasisth1,sol,basisxth1,...
basisyth1,bcvals,helpbc,dirval,rho,probdata,coefft);    

%quadrilaterials
[sol]=assembleqh1int(nelemq,Meshq,order,map,irn,icn,nz,glob,nunk,dir,nunki,...
edges,unkv,unkvh,unkvi,ndir,x,w,sbasish1,rbasish1,basisxh1,basisyh1,sol,...
unkqih,helpbc,bcvals,rho,dirval,probdata,coeff);
end


disp('interiors complete, post processing begins')

% include Dirichlet values in the solution vector
[sol,unkv,nunki,unkvh]=streat(sol,dirval,unkv,unkvh,nunki,order);

% split the mesh ready for a contour plot
if plotson==1
[Meshsq,Meshst,splitq,orders,splitt]=split(ne,order,Meshq,nelemq,Mesht,nelemt,...
dir,glob,edges,coeff,coefft);

% make contour plot of solution on the split mesh
contourneub(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,dir,glob,Meshq,Mesht,sol,splitt,unkv,...
unkvh,unkvi,unkqih,edges,coeff,coefft);
end

% measure the L2 norm of the error
if probdata.exact==1
errengy=myerror(sol,sbasish1,nelemq,Meshq,order,unkv,unkvh,...
unkqih,unkvi,glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,...
basisxth1,basisyth1,sbasisth1,nelemt,basisxh1,basisyh1,probdata,coeff,coefft);
end	 
	 
if probdata.presscoeff==1
pressurecoef(nelemt,nelemq,order,splitq,dir,glob,Meshq,Mesht,sol,splitt,unkv,xt,wt,...
intxi,inteta,intw,unkvh,unkvi,sbasisth1,edges,probdata,coeff,coefft)
end
