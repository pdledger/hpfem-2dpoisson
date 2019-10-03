function [sol]=assembleqh1int(nelemq,Meshq,order,map,irn,icn,nz,glob,nunk,dir,...
nunki,edges,unkv,unkvh,unkvi,ndir,x,w,sbasish1,rbasish1,basisxh1,basisyh1,sol,...
unkqih,helpbc,bcvals,rho,dirval,probdata,coeff)

% determine the interior dofs for triangular elements 
% outputs:
% sol: global sol vector including interior dofs for triangles

% extract coordinate and connectivity information
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
fnum=Meshq.Fnum;

% Material data
mat=probdata.mat;

% set size of local arrays
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
      lbc(j)=edges(glob(i,j),3);
   end
   
%  material constant in this element
   lmat=mat(fnum(i));

% compute the elemental mass matrix 
  massel=massqh1(xy,order,x,w,sbasish1);

% compute the elemental stiffness matrix  
  stiffel=stiffqh1(xy,order,x,w,basisxh1,basisyh1,lmat);
  
% compute the elemental rhs vectors  
  rhsel=rhsqh1(xy,order,x,w,rbasish1,sbasish1,rho,lbc,probdata);

% vertex functions
 for j=1:4
     lunkv(j)=unkv(intmaq(i,j));
     ldrv(j)=1;
  end
% Higher order functions
if order >0
% Higher order edges
  for j=1:4
  for p=0:order-1
  lunkv(4+j+4*(p))=unkvh(glob(i,j),p+1);
  ldrv(4+j+4*(p))=dir(i,j)^(p+2);
  end
  end
% Interior functions
  for j=1:order^2
  lunkv(4+4*order+j)=unkqih(i,j);
  ldrv(4+4*(order)+j)=1;    
  end
end  

lval=zeros(esizeh1,1);
for j=1:esizeh1
if lunkv(j) > 0
lval(j)=sol(lunkv(j))*ldrv(j);
elseif lunkv(j)< 0
lval(j)=dirval(abs(lunkv(j)));
end
end

% perform static condenstation
[lvali]=static2h1(massel,rhsel,order,lval,esizeh1full,esizeh1,stiffel);

% determine interior dofs  
  for j=1:order^2
  row=lunkv(j+esizeh1);
  dr=ldrv(j+esizeh1);
  if row > 0
  sol(row)=lvali(j)*dr;
  end
  end
  
end


