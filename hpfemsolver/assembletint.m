function [sol]=assembletint(nelemt,Mesht,...
        order,glob,nunk,dir,nunki,edges,nelemq,unkv,unkvh,unkvi,...
       intxi,inteta,intw,xt,wt,sbasisth1,...
       rbasisth1,sol,basisxth1,basisyth1,bcvals,helpbc,dirval,rho,probdata,coefft)

% determine the interior dofs for triangular elements 
% outputs:
% sol: global sol vector including interior dofs for triangles

% extract coordinate and connectivity information
coord=Mesht.Coordinates;
intmat=Mesht.Elements;

% set size of local arrays
esizeh1=3+3*order;
esizeh1full=3+3*order+(order-1)*order/2;

fnum=Mesht.Fnum;

% Material data
mat=probdata.mat;

for i=1:nelemt

% extract elemental coordinates and elemental bc vector
   for j=1:3
      for k=1:2
         xyt(j,k)=coord(intmat(i,j),k);
      end
   end
   
   for j=1:3
      lbc(j)=edges(glob(i+nelemq,j),3);
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

   lmat=mat(fnum(i));
   
% compute the elemental mass matrix 
  massel=massth1(xyt,order,intxi,inteta,intw,sbasisth1,localblendt,geomflag);

% compute the elemental stiffness matrix  
  stiffel=stiffth1(xyt,order,intxi,inteta,intw,basisxth1,basisyth1,lmat,localblendt,geomflag);

% compute the elemental rhs vector  
 
rhsel=rhsth1(xyt,order,intxi,inteta,intw,basisxth1,basisyth1,xt,wt,sbasisth1,rbasisth1,lbc,...
rho,bcvals,probdata,localblendt,geomflag);
   
% lowest order block
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
  lunkv(3+3*order+j)=unkvi(i,j);
  ldrv(3+3*(order)+j)=1;    
  end
end  

lval=zeros(esizeh1,1);
for j=1:esizeh1
if lunkv(j)> 0
lval(j)=sol(lunkv(j))*ldrv(j);
elseif lunkv(j) < 0
lval(j)=dirval(abs(lunkv(j)));
end
end

% perform static condenstation
[lvali]=static2h1(massel,rhsel,order,lval,esizeh1full,esizeh1,stiffel);

  
% determine interior dofs  
  for j=1:(order-1)*order/2 % nb changed from cont
     row=lunkv(j+esizeh1);  
     dr=ldrv(j+esizeh1);
     if row > 0
	sol(row)=lvali(j)*dr;
     end
  end 


end

