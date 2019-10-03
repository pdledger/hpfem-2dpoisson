function J=mapquad(xy,x,y)

dlam=(1/sqrt(3))*[sqrt(3)/2 -1/2 ; 0 1; -sqrt(3)/2 -1/2];

J=zeros(2);
for i=1:2
for j=1:2
for k=1:3
J(i,j)=J(i,j)+(dlam(k,j)*xy(k,i));
end
end
end
