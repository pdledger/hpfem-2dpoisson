function ph=basistl2(order,x,y)

esizel2=1+((order-1)*(order)/2+order-1);
%order=order+1;
if order < 0
esizel2=1;
end
ph=zeros(esizel2,1);
% constant function
%ph(1)=1;


% Vertex functions
vph(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
vph(2)=y/sqrt(3);
vph(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);
gph(1,1)=1/2;
gph(1,2)=-1/(2*sqrt(3));

gph(2,1)=0;
gph(2,2)=1/sqrt(3);

gph(3,1)=-1/2;
gph(3,2)=-1/(2*sqrt(3));

for i=1:3
curl(i,1)=gph(i,2);
curl(i,2)=-gph(i,1);
end

edges=[ 1 2; 2 3; 3 1];

ph(1)=0;
for i=1:3
ph(1)=ph(1)+(curl(edges(i,1),:)*gph(edges(i,2),:)')-(curl(edges(i,2),:)*gph(edges(i,1),:)');
end


% changed to match H(div functions)

% div ( H(div)) (interior functions only!!!)

% type 1

if order > 1
basno=1;
l(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
l(2)=y/sqrt(3);
l(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);
dl=[1/2 -1/(2*sqrt(3));
    0 1/sqrt(3);
    -1/2 -1/(2*sqrt(3))];    

s=l(2)-l(1);
t=l(1)+l(2);
ds=[dl(2,1)-dl(1,1), dl(2,2)-dl(1,2)];
dt=[dl(2,1)+dl(1,1), dl(2,2)+dl(1,2)];

% Type 1
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
%basno=basno+1;

if t==0
  if i>= 1
    duix=0;
    duiy=0;
    ui=0;
  else
    duix=s*ds(1);
    duiy=s*ds(2);
    ui=s^2/2;
  end
else
  ui=legi(s/t,i+2)*t^(i+2);
  duix=-leg(s/t,i)*t^i*t*dt(1)+leg(s/t,i+1)*t^(i+1)*ds(1);
  duiy=-leg(s/t,i)*t^i*t*dt(2)+leg(s/t,i+1)*t^(i+1)*ds(2);
end

vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));
% this function is zero!!!
%ph(basno)=(duiy*dvjx-duix*dvjy)+(duix*dvjy-duiy*dvjx);

end
end
end

% Type 2
for i=0:1:order-2
for j=0:1:order-2
if i+j <= order-2
basno=basno+1;
if t==0
  if i>= 1
    duix=0;
    duiy=0;
    ui=0;
  else
    duix=s*ds(1);
    duiy=s*ds(2);
    ui=s^2/2;
  end
else
  ui=legi(s/t,i+2)*t^(i+2);
  duix=-leg(s/t,i)*t^i*t*dt(1)+leg(s/t,i+1)*t^(i+1)*ds(1);
  duiy=-leg(s/t,i)*t^i*t*dt(2)+leg(s/t,i+1)*t^(i+1)*ds(2);
end

vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));


ph(basno)=(duiy*dvjx-duix*dvjy)-(duix*dvjy-duiy*dvjx);

end
end
end

% Type 3
for j=0:1:order-2

basno=basno+1;
vj=l(3)*leg(l(3)-l(2)-l(1),j);
dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));
db1x=vj*dl(2,1)+dvjx*l(2);
db1y=vj*dl(2,2)+dvjy*l(2);

da1x=dl(1,1);
da1y=dl(1,2);

db2x=vj*dl(1,1)+dvjx*l(1);
db2y=vj*dl(1,2)+dvjy*l(1);

da2x=dl(2,1);
da2y=dl(2,2);

ph(basno)=(db1x*da1y-db1y*da1x)-(db2x*da2y-db2y*da2x);
end
end



[m]=length(ph);
if m~=esizel2
  disp('error wrong number of basis functions l2basis')
  order
  pause
end
