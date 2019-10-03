function stiffel=stiffqh1(xy,order,x,w,basisxh1,basisyh1,lmat)

% Determine the elemental stiffness matrix for quadrilaterials
% output:
% stiffel: elemental stiffness matrix for quadrilaterials

nip=length(x);

esize=(order+1+1)^2;
stiffel=zeros(esize);
for p=1:nip
for q=1:nip

phx=basisxh1(:,q+(nip*(p-1)));
phy=basisyh1(:,q+(nip*(p-1)));
ph=[phx phy];
J=mapquad(xy,x(p),x(q));
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
phxy=Jinv*ph';
phxy=phxy';

for i=1:esize
for j=1:esize
for k=1:2
stiffel(i,j)=stiffel(i,j)+(w(p)*w(q)*detJ*phxy(i,k)*phxy(j,k)*lmat);
end
end
end

end
end
