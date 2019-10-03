function ph=divgbasisqh1(order,x,y);

% Evaluate the gradient of the H1 basis for quadrilaterials at the point (x,y)
% output:
% ph: the gradient of the basis functions appropriate to degree order

% NB order for H^1 basis is 1 greater than that for H(curl)

% p=1 Vertex Functions
ph(1,1)=0;
ph(1,2)=1;
ph(1,3)=1;
ph(1,4)=0;

ph(2,1)=0;
ph(2,2)=-1;
ph(2,3)=-1;
ph(2,4)=0;

ph(3,1)=0;
ph(3,2)=1;
ph(3,3)=1;
ph(3,4)=0;

ph(4,1)=0;
ph(4,2)=-1;
ph(4,3)=-1;
ph(4,4)=0;


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

ph(4*p+e,1)=ddlegi(ep(e),p-1+2)*dep(e,1)*la(e)*dep(e,1)+...
            dlegi(ep(e),p-1+2)*dep(e,1)*dla(e,1)+...
            dlegi(ep(e),p-1+2)*dla(e,1)*dep(e,1);

ph(4*p+e,2)=ddlegi(ep(e),p-1+2)*dep(e,1)*la(e)*dep(e,2)+...
            dlegi(ep(e),p-1+2)*dep(e,1)*dla(e,2)+...
            dlegi(ep(e),p-1+2)*dla(e,1)*dep(e,2);
	    

ph(4*p+e,3)=ddlegi(ep(e),p-1+2)*dep(e,2)*la(e)*dep(e,1)+...
            dlegi(ep(e),p-1+2)*dep(e,2)*dla(e,1)+...
	    dlegi(ep(e),p-1+2)*dla(e,2)*dep(e,1);  
	    
ph(4*p+e,4)=ddlegi(ep(e),p-1+2)*dep(e,2)*la(e)*dep(e,2)+...
            dlegi(ep(e),p-1+2)*dep(e,2)*dla(e,2)+...
	    dlegi(ep(e),p-1+2)*dla(e,2)*dep(e,2);  

end
end    

end
edgefun=4*order;

% Interior Functions
if order > 1
basno=edgefun;

% Type 1
for i=0:1:order-2
for j=0:1:order-2
basno=basno+1;

ph(basno,1)=legi(2*y-1,j+2)*ddlegi(2*x-1,i+2)*2*2;
ph(basno,2)=2*dlegi(2*y-1,j+2)*dlegi(2*x-1,i+2)*2;

ph(basno,3)=2*dlegi(2*x-1,i+2)*dlegi(2*y-1,j+2)*2;
ph(basno,4)=legi(2*x-1,i+2)*ddlegi(2*y-1,j+2)*2*2;
end
end


end

[m,n]=size(ph);
if m~=(order+1)^2
m , (order+1)^2
  disp('error wrong number of basis functions')
end
  
