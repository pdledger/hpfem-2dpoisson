function [rhs,stiff,dirval]=assembleth1(nelemt,Mesht,...        
order,nz,glob,nunk,dir,nunki,edges,nelemq,unkv,unkvh,unkvi,ndir,helpbc,I,J,X,bcvals,rhs,dirval,rho,probdata,basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1,...
intxi,inteta,intw,xt,wt,coefft)

% assemble the triangular elements performing static condensation on the fly
% outputs:
% rhs: global rhs vector
% stiff: assembled stiffness matrix
% intxi: integration points for triangular elements, xi component
% inteta: integration points for triangular elements, eta component
% intw: integration points for triangular elements, weights
% xt: integration points for the interval [-1,1]
% wt: integration weights for the interval [-1,1]
% sbasisth1: basis functions evaluated at integrations inside (triangular) reference element
% rbasisth1: basis functions evaluated at integration points on edges of (triangular) reference element
% basisxth1: x component of the gradient of basis functions evaluated at integration
% points inside triangular element
% basisyth1: y component of the gradient of basis functions evaluated at integration
% points inside triangular element
% dirval: Dirichlet boundary condition vector

% extract coordinate and connectivity information
coord=Mesht.Coordinates;
intmat=Mesht.Elements;
fnum=Mesht.Fnum;

% Material data
mat=probdata.mat;

[npoin dum]=size(coord);

% compute the integration weights and locations 
%[intxi,inteta,intw]=gautri(2*(order+2));

% NB for edges use different set to those used with quadrilaterials! 
%[xt,wt]=gaussquad(-1.,1,2*(order+2)); 

% work out basis at each integration point
%if nelemt ~= 0
%[basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1]=myevalth1(order,intxi,inteta,xt);
%else
%basisxth1=[];
%basisyth1=[];
%sbasisth1=[];
%rbasisth1=[];
%divbasisxxth1=[];
%divbasisxyth1=[];
%divbasisyxth1=[];
%divbasisyyth1=[];
%end

% get dirichlet values on vertices and edges
dirval=dirichlett(bcvals,npoin,ndir,unkv,helpbc,xt,wt,rbasisth1,unkvh,Mesht,...
edges,glob,nelemt,order,nelemq,dirval,probdata,coefft);

% set size of local arrays
esizeh1=3+3*order;
esizeh1full=3+3*order+(order-1)*order/2;



for i=1:nelemt

% extract elemental coordinates and elemental bc vector
   for j=1:3
      for k=1:2
         xyt(j,k)=coord(intmat(i,j),k);
      end
   end
   
   geomflag=0;
   for j=1:3
      for k=1:2
      localblendt(j,k)=coefft(i,j,k);
      if localblendt(j,k) > 0
      geomflag=1;
      end
      end
   end

   
   for j=1:3
      lbc(j)=edges(glob(i+nelemq,j),3);
   end
   
   lmat=mat(fnum(i));

% compute the elemental mass matrix 
  massel=massth1(xyt,order,intxi,inteta,intw,sbasisth1,localblendt,geomflag);

% compute the elemental stiffness matrix  
  stiffel=stiffth1(xyt,order,intxi,inteta,intw,basisxth1,basisyth1,...
  lmat,localblendt,geomflag);

% compute the elemental rhs vector   
  rhsel=rhsth1(xyt,order,intxi,inteta,intw,basisxth1,basisyth1,xt,...
  wt,sbasisth1,rbasisth1,lbc,rho,bcvals,probdata,localblendt,geomflag);

% perform static condenstation
  [kcc,rc]=static1h1(massel,rhsel,order,esizeh1full,esizeh1,stiffel);

% begin assembly proceedure 
% create a list of unknown numbers in the order of the elemental matrix
% lowest order block
% Vertex based dof's  
 for j=1:3
     lunkv(j)=unkv(intmat(i,j));
     ldrv(j)=1;
  end

  % Higher order edges
if order >0
  for j=1:3
  for p=0:order-1
  lunkv(3+j+3*(p))=unkvh(glob(i+nelemq,j),p+1);
  ldrv(3+j+3*(p))=dir(i+nelemq,j)^(p+2);
  end
  end
  % Interior block
  for j=1:order*(order-1)/2
  lunkv(3+3*order+j)=0;
  ldrv(3+3*(order)+j)=1;    
  end
end  

% Now assemble matrices in sparse format
  
  for j=1:esizeh1 % nb changed from cont
     row=lunkv(j);  
     dr=ldrv(j);
     if row > 0
        for k=1:esizeh1 % nb changed from cont
           col=lunkv(k);
           dc=ldrv(k);
           if col > 0
             nz=nz+1;
             len=length(X);
             if nz > len
               I(2*len)=0;
               J(2*len)=0;
               X(2*len)=0;
               Xabs(2*len)=0;
             end
             I(nz)=row;
             J(nz)=col;
             X(nz)=kcc(j,k)*dr*dc;
	      elseif col < 0
	     rhs(row)=rhs(row)-(kcc(j,k))*dr*dc*dirval(abs(col));
	   end
        end
	rhs(row)=rhs(row)+(dr*rc(j));
     end
  end 
end


% build stiffness matrix from I,J and X
stiff=sparse(I(1:nz),J(1:nz),X(1:nz),nunk,nunk);

clear I
clear J
clear X
