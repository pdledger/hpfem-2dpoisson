function ph=basisqh1(order,x,y);

% Evaluate the H1 basis for quadrilaterials at the point (x,y)
% output:
% ph: set of basis functions appropriate to degree order

% NB order for H^1 basis is 1 greater than that for H(curl)
ph=zeros((order+1)^2,1);

% p=1 Vertex Functions
ph(1)=(1-x)*(1-y);
ph(2)=x*(1-y);
ph(3)=x*y;
ph(4)=(1-x)*y;

% p> 0 Edge Functions
if order > 1
ep=[2.*x-1;
   1-2.*x;
   2.*y-1;
   1-2.*y];

dep=[2 0;
     -2 0;
      0 2;
      0 -2];

la=[1-y;
    y;
    x;
    1-x];
dla=[0 -1;
     0 1;
     1 0;
    -1 0];

for p=1:order-1
for e=1:4
ph(4*p+e)=legi(ep(e),p-1+2)*la(e);
end
end    

end

edgefun=4*(order);

% Interior Functions
if order > 1
basno=edgefun;

% Type 1
for i=0:1:order-2
for j=0:1:order-2
basno=basno+1;
ph(basno)=legi(2*y-1,j+2)*legi(2*x-1,i+2);
end
end


end

[m,n]=size(ph);
if m~=(order+1)^2
m , (order+1)^2
  disp('error wrong number of basis functions')
end
  
