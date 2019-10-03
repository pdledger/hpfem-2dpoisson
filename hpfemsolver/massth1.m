function massel=massth1(xy,order,intxi,inteta,intw,sbasisth1,localblendt,geomflag)

% Determine the elemental mass matrix for triangles
% output:
% massel: elemental mass matrix for triangles

nip=length(intxi);
esize=3+3*order+(order-1)*order/2;
massel=zeros(esize);

if geomflag==0
% linear geometry
for p=1:nip

ph=sbasisth1(:,p);

J=maptrilin(xy,intxi(p),inteta(p));
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';

%for i=1:esize
%for j=1:esize
%massel(i,j)=massel(i,j)+(intw(p)*detJ*ph(i)*ph(j));
%end
%end

massel=massel+(intw(p)*detJ*ph*ph');


end

else

for p=1:nip

ph=sbasisth1(:,p);

J=maptri(xy,intxi(p),inteta(p),localblendt);
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';

%for i=1:esize
%for j=1:esize
%massel(i,j)=massel(i,j)+(intw(p)*detJ*ph(i)*ph(j));
%end
%end

massel=massel+(intw(p)*detJ*ph*ph');


end

end
