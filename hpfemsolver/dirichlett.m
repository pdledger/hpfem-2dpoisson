function dirval=dirichlett(bcvals,npoin,ndir,unkv,helpbc,xt,wt,rbasisth1,unkvh,...
Mesht,edges,glob,nelemt,order,nelemq,dirval,probdata,coefft)

% Determine the Dirichlet bounadary condition values at vertices and on edges
% assuming contant values.
% output:
% dirval: vector of Dirichlet conditions

% Intialize arrays
nip=length(xt);
coord=Mesht.Coordinates;
intmat=Mesht.Elements;

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
for i=1:nelemt

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
   if edges(glob(i+nelemq,j),3) > 1
   % this is a Dirichlet edge
   a=zeros(order);
   r=zeros(order,1);
   
   for p=1:nip   
   %sam as fortran code...
   if j==1
     eta=sqrt(3.)*0.5*(xt(p)+1);
     xi=(-1+xt(p))*(-0.5);
   elseif j==2
     eta=sqrt(3.)*(-0.5)*(xt(p)-1);
     xi=-1*(1+xt(p))*0.5;
   else
    eta=0.;
    xi=xt(p);
   end
   ph=rbasisth1(:,p+(nip*(j-1)));

   J=maptri(xyt,xi,eta,localblendt);
   detJ=det(J);
   Jinv=(inv(J))';

   %NB diffferent from fortran code!!!
   if j==1
     ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
     ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
     dete=sqrt((ddgx^2)+(ddgy^2));
   elseif j==2
      ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
      ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
      dete=sqrt((ddgx^2)+(ddgy^2));
   else
      ddgx=J(1,1);
      ddgy=J(2,1);
      dete=sqrt((ddgx^2)+(ddgy^2));
   end   
   [xp,yp]=getcoordt(xyt,xi,eta,localblendt);

%  compute first order approximation
   vlow=0;
%  vertex basis functions not associated with this edge vanish
%  when evaluated on Dirichlet edge   
   for pp=1:3
   if unkv(intmat(i,pp)) < 0
   vlow=vlow+(ph(pp)*dirval(abs(unkv(intmat(i,pp)))));
   end
   end

% evaluate Dirichlet function
fun=probdata.dirfun;
arg=probdata.dirfunarg;
index=edges(glob(i,j),3);
vex=fun(xp,yp,index,arg);
 
   
%  compute L-2 projection for obtianing Dirichlet values
   for pp=1:order
   for ppp=1:order
   a(pp,ppp)=a(pp,ppp)+(wt(p)*dete*ph(j+3*(pp))*ph(j+3*(ppp)));
   end
   r(pp)=r(pp)-(wt(p)*dete*ph(j+3*(pp))*(vlow-vex));
   end
   
   end
   s=a\r;
   
   for pp=1:order
       if abs(unkvh(glob(i+nelemq,j),pp))~=0
   dirval(abs(unkvh(glob(i+nelemq,j),pp)))=s(pp);
       end
   end
   end
   
   end
end   
end
