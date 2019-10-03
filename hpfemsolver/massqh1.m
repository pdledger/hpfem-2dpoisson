function massel=massqh1(xy,order,x,w,sbasish1)

% Determine the elemental mass matrix for quadrilaterials
% output:
% massel: elemental mass matrix for quadrilaterials

nip=length(x);
esize=(order+1+1)^2;
massel=zeros(esize);

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

for i=1:esize
for j=1:esize
massel(i,j)=massel(i,j)+(w(p)*w(q)*detJ*ph(i)*ph(j));
end
end

end
end
