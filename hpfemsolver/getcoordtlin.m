function [x,y]=getcoordtlin(xy,xi,eta)
la(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*xi-eta);
la(2)=eta/sqrt(3);
la(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*xi-eta);
x=0;
y=0;
for j=1:3
  x=x+(xy(j,1)*la(j));
  y=y+(xy(j,2)*la(j));
end

