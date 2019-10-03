function [rhsel]=rhsth1(xyt,order,intxi,inteta,intw,basisxth1,basisyth1,xt,wt,...
sbasisth1,rbasisth1,lbc,rho,bcvals,probdata,localblendt,geomflag)

% Determine the elemental rhs vector for triangular elements
% output:
% rhsel: elemental rhs vector for triangular elements

% intialisations
nip=length(xt);
esizeh1=3+3*order+(order-1)*order/2;
rhsel=zeros(esizeh1,1);
norml=[0.5*sqrt(3.) 0.5; -0.5*sqrt(3.) 0.5; 0 -1.]';

bctype=probdata.bctype;

if geomflag==0

% apply Neumann condition if appropriate
for j=1:3
  if lbc(j)~=0 & bctype(lbc(j))==1

% This is a boundary edge and a boundary integral is to be applied   

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
   
   J=maptrilin(xyt,xi,eta);
   detJ=det(J);
   Jinv=(inv(J))';

   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
   %NB diffferent from fortran code!!!
   if j==1
     ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
     ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
     v=sqrt((ddgx^2)+(ddgy^2));
   elseif j==2
      ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
      ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
      v=sqrt((ddgx^2)+(ddgy^2));
   else
      ddgx=J(1,1);
      ddgy=J(2,1);
      v=sqrt((ddgx^2)+(ddgy^2));
   end   
   [xp,yp]=getcoordt(xyt,xi,eta,localblendt);
   
    fun=probdata.neufun;
    arg=probdata.neufunarg;
    index=lbc(j);
    val=fun(xp,yp,index,nm,arg);
    for k=1:esizeh1
     rhsel(k)=rhsel(k)+(val*ph(k)*wt(p)*v);
     end

   
   
   end
 end
end       

% apply source term
nip=length(intxi);
for p=1:nip
ph=sbasisth1(:,p);
J=maptrilin(xyt,intxi(p),inteta(p));
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end

[xp,yp]=getcoordtlin(xyt,intxi(p),inteta(p));

    fun=probdata.srcfun;
    arg=probdata.srcfunarg;
    val=fun(xp,yp,arg);

srho=val;

Jinv=(inv(J))';
for i=1:esizeh1
rhsel(i)=rhsel(i)+(intw(p)*detJ*ph(i)*srho);
end
end


else

% apply Neumann condition if appropriate
for j=1:3
  if lbc(j)~=0 & bctype(lbc(j))==1

% This is a boundary edge and a boundary integral is to be applied   

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

   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
   %NB diffferent from fortran code!!!
   if j==1
     ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
     ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
     v=sqrt((ddgx^2)+(ddgy^2));
   elseif j==2
      ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
      ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
      v=sqrt((ddgx^2)+(ddgy^2));
   else
      ddgx=J(1,1);
      ddgy=J(2,1);
      v=sqrt((ddgx^2)+(ddgy^2));
   end   
   [xp,yp]=getcoordt(xyt,xi,eta,localblendt);
   
    fun=probdata.neufun;
    arg=probdata.neufunarg;
    index=lbc(j);
    val=fun(xp,yp,index,nm,arg);
    for k=1:esizeh1
     rhsel(k)=rhsel(k)+(val*ph(k)*wt(p)*v);
     end

   
   
   end
 end
end       

% apply source term
nip=length(intxi);
for p=1:nip
ph=sbasisth1(:,p);
J=maptri(xyt,intxi(p),inteta(p),localblendt);
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end

[xp,yp]=getcoordt(xyt,intxi(p),inteta(p),localblendt);

    fun=probdata.srcfun;
    arg=probdata.srcfunarg;
    val=fun(xp,yp,arg);

srho=val;

Jinv=(inv(J))';
for i=1:esizeh1
rhsel(i)=rhsel(i)+(intw(p)*detJ*ph(i)*srho);
end
end

end
