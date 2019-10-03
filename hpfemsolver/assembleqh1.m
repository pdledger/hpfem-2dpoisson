function [rhs,I,J,X,dirval,nz]=assembleqh1(nelemq,Meshq,...        
order,glob,nunk,dir,nunki,edges,unkv,unkvh,unkvi,ndir,unkqih,x,w,helpbc,bcvals,rho,probdata,basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,...
divbasisyxh1,divbasisyyh1,coeff)

% assemble the quadrilaterial elements performing static condensation on the fly
% outputs:
% rhs: global rhs vector
% I: stiffness row vector information
% J: stiffness column vector information
% X: stiffness non-zero information(I,J,X will be combined to form stiffness matrix)
% sbasish1: basis functions evaluated at integrations inside reference element
% rbasish1: basis functions evaluated at integration points on edges of reference element
% basisxh1: x component of the gradient of basis functions evaluated at integration
% points inside element
% basisyh1: y component of the gradient of basis functions evaluated at integration
% points inside element
% dirval: Dirichlet boundary condition vector
% nz: Number of non-zeros


% extract coordinate and connectivity informaion
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
[npoin dum]=size(coord);
fnum=Meshq.Fnum;

% Material data
mat=probdata.mat;

% evaluate the basis functions at the integration points across a reference element
%[basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1]=myevalh1(order,x);

% Number arrays ready for assembly
rhs=zeros(nunk,1);
I=zeros(nunk,1);
J=zeros(nunk,1);
X=zeros(nunk,1);
nz=0;

% get dirichlet values on vertices and edges
dirval=dirichletq(bcvals,npoin,ndir,unkv,helpbc,x,w,rbasish1,unkvh,Meshq,edges,glob,nelemq,order,probdata);

% set size of elemental arrays
esizeh1=4+4*order;
esizeh1full=(order+1+1)^2;

for i=1:nelemq

% extract elemental coordinates and elemental bc vector
   for j=1:4
      for k=1:2
         xy(j,k)=coord(intmaq(i,j),k);
      end
   end
   
  for j=1:4
      for k=1:2
      localblend(j,k)=coeff(i,j,k);
      end
   end

   
   for j=1:4
      lbc(j)=edges(glob(i,j),3);
   end

%  material constant in this element
   lmat=mat(fnum(i));

% compute the elemental mass matrix (update for curved geometry) 
  massel=massqh1(xy,order,x,w,sbasish1);

% compute the elemental stiffness matrix  (update for cuvred geometry_
  stiffel=stiffqh1(xy,order,x,w,basisxh1,basisyh1,lmat);
  
% compute the elemental rhs vector  (update for cuvred geometry_
  rhsel=rhsqh1(xy,order,x,w,rbasish1,sbasish1,rho,lbc,probdata);

% perform static condenstation
  [kcc,rc]=static1h1(massel,rhsel,order,esizeh1full,esizeh1,stiffel);

% begin assembly proceedure 
% create a list of unknown numbers in the order they appear in the elemental matrix
% Vertex based dof's  
 for j=1:4
     lunkv(j)=unkv(intmaq(i,j));
     ldrv(j)=1;
  end
% Higher order functions
if order >0
  % edge block
  for j=1:4
  for p=0:order-1
  lunkv(4+j+4*(p))=unkvh(glob(i,j),p+1);
  ldrv(4+j+4*(p))=dir(i,j)^(p+2);
  end
  end
  % Interior block
  for j=1:order^2
  lunkv(4+4*order+j)=0;
  ldrv(4+4*(order)+j)=1;    
  end
end  

% Now assemble matrices in sparse format
  

  for j=1:esizeh1 % nb changed from cont
     row=lunkv(j);  
     dr=ldrv(j);
     if row > 0
% Edge elements
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

