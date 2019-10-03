function ph=basisth1(order,x,y);
esizeh1=3+3*(order-1)+(order-2)*(order-1)/2;

if x < -1
x=-1;
end
if x > 1
x=1;
end
if y < 0
y=0;
end
if y > sqrt(3)
y=sqrt(3);
end

ph=zeros(esizeh1,1);
% NB order for H^1 basis is 1 greater than that for H(curl)
% Vertex functions
ph(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
ph(2)=y/sqrt(3);
ph(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);

% Edge functions
if order >1
s=[(1/2*y*sqrt(3))-1/2-(1/2*x);
   1/2-(1/2*x)-1/2*y*sqrt(3);
   x];
t=[1/2+1/2*x+1/6*y*sqrt(3);
   1/2-1/2*x+1/6*y*sqrt(3);
   1-1/3*y*sqrt(3)];
for p=0:order-2
for e=1:3
if t(e)~=0
a=legi(s(e)/t(e),p+2)*t(e)^(p+2);
else
a=0;
end
ph(3*(p+1)+e)=a;
end
end


l(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
l(2)=y/sqrt(3);
l(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);
s=l(2)-l(1);
t=l(1)+l(2);
basno=3*(order);
for i=0:1:order-3
for j=0:1:order-3
if i+j <= order-3
basno=basno+1;
ui=legi(s/t,i+2)*t^(i+2);
vj=l(3)*leg(l(3)-l(2)-l(1),j); 
ph(basno)=ui*vj;
end
end
end

end
