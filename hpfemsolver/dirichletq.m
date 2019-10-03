function dirval=dirichletq(bcvals,npoin,ndir,unkv,helpbc,x,w,rbasish1,unkvh,Meshq,edges,glob,nelemq,order,probdata)

% Determine the Dirichlet bounadary condition values at vertices and on edges
% assuming contant values.
% output:
% dirval: vector of Dirichlet conditions

% Intialize arrays
dirval=zeros(ndir,1);
nip=length(x);
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;

% Determine vertex based Dirichlet condition values
for i=1:npoin
if unkv(i) < 0
fun=probdata.dirfun;
arg=probdata.dirfunarg;
xp=coord(i,1);
yp=coord(i,2);
index=helpbc(i);
value=fun(xp,yp,index,arg);
dirval(abs(unkv(i)))=value;

end
end

% Determine edge based Dirichlet values by doing an L2 minimisation
if order > 0
for i=1:nelemq

   for j=1:4
      for k=1:2
         xyq(j,k)=coord(intmaq(i,j),k);
      end
   end

   for j=1:4
   if edges(glob(i,j),3) > 1
   % this is a Dirichlet edge
   a=zeros(order);
   r=zeros(order,1);
   
   for p=1:nip
   % obtain basis functions
     if j==1
     xi=x(p);
     et=0;
   elseif j==2
     xi=x(p);
     et=1;
   elseif j==3
     xi=1;
     et=x(p);
   else
     xi=0;
     et=x(p);
   end
   ph=rbasish1(:,p+(nip*(j-1)));

   J=mapquad(xyq,xi,et);
   detJ=det(J);
   Jinv=(inv(J))';
   
   if j==1 || j==2
   dete=sqrt(J(1,1)^2+J(2,1)^2);
   else
   dete=sqrt(J(1,2)^2+J(2,2)^2);
   end

%  compute first order approximation
   vlow=0;
%  vertex basis functions not associated with this edge vanish
%  when evaluated on Dirichlet edge   
   for pp=1:4
   if unkv(intmaq(i,pp)) < 0
   vlow=vlow+(ph(pp)*dirval(abs(unkv(intmaq(i,pp)))));
   end
   end


   [xp,yp]=getcoord(xyq,xi,et);

% evaluate Dirichlet function
fun=probdata.dirfun;
arg=probdata.dirfunarg;
index=edges(glob(i,j),3);
vex=fun(xp,yp,index,arg);

   
%  compute L-2 projection for obtianing Dirichlet values
   for pp=1:order
   for ppp=1:order
   a(pp,ppp)=a(pp,ppp)+(w(p)*dete*ph(j+4*(pp))*ph(j+4*(ppp)));
   end
   r(pp)=r(pp)-(w(p)*dete*ph(j+4*(pp))*(vlow-vex));
   end
   
   end
   s=a\r;
   for pp=1:order
   if abs(unkvh(glob(i,j),pp))~=0
   dirval(abs(unkvh(glob(i,j),pp)))=s(pp);
   end
   end
   end
   
   end
end   
end
