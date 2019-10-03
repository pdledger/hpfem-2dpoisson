function rhsel=rhsqh1(xy,order,x,w,rbasish1,sbasish1,rho,lbc,probdata)

% Determine the elemental rhs vector for quadrilaterial elements
% output:
% rhsel: elemental rhs vector for quadrilaterial elements

% intialisations
nip=length(x);
esize=(order+1+1)^2;
rhsel=zeros(esize,1);
norml=[0 -1; 0 1; 1 0; -1 0]';
bctype=probdata.bctype;

% apply Neumann condition if appropriate
for j=1:4
   if lbc(j)~=0 & bctype(lbc(j))==1

% This is a boundary edge and a boundary integral is to be applied   

   for p=1:nip
   
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
   
   J=mapquad(xy,xi,et);
   detJ=det(J);
   Jinv=(inv(J))';
   
   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
   
   if j==1 || j==2
   v=sqrt(J(1,1)^2+J(2,1)^2);
   else
   v=sqrt(J(1,2)^2+J(2,2)^2);
   end
   [xp,yp]=getcoord(xy,xi,et);


    fun=probdata.neufun;
    arg=probdata.neufunarg;
    index=lbc(j);
    val=fun(xp,yp,index,nm,arg);
    for k=1:esize
     rhsel(k)=rhsel(k)+(val*ph(k)*w(p)*v);
     end
     
   end
 end
end

% apply source term
for p=1:nip
for q=1:nip

ph=sbasish1(:,q+(nip*(p-1)));

J=mapquad(xy,x(p),x(q));
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
[xp,yp]=getcoord(xy,x(p),x(q));
    fun=probdata.srcfun;
    arg=probdata.srcfunarg;
    val=fun(xp,yp,arg);

srho=val;

for i=1:esize
rhsel(i)=rhsel(i)+(w(p)*w(q)*detJ*ph(i)*srho);
end


end
end
