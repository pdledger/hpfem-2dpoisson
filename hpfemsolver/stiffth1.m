function stiffel=stiffth1(xy,order,intxi,inteta,intw,basisxth1,basisyth1,lmat,localblendt,geomflag)

% Determine the elemental stiffness matrix for triangles
% output:
% stiffel: elemental stiffness matrix for triangles

nip=length(intxi);
esize=3+3*order+(order-1)*order/2;
stiffel=zeros(esize);

if geomflag==0

for p=1:nip

phx=basisxth1(:,p);
phy=basisyth1(:,p);
ph=[phx phy];
J=maptrilin(xy,intxi(p),inteta(p));
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
phxy=Jinv*ph';
phxy=phxy';

kappa=lmat;


%for i=1:esize
%for j=1:esize
%for k=1:2
%stiffel(i,j)=stiffel(i,j)+(kappa*intw(p)*detJ*phxy(i,k)*phxy(j,k));
%end
%end

stiffel=stiffel+(kappa*intw(p)*detJ*phxy*phxy');

%end

%end
end

else

for p=1:nip

phx=basisxth1(:,p);
phy=basisyth1(:,p);
ph=[phx phy];
J=maptri(xy,intxi(p),inteta(p),localblendt);
detJ=det(J);
if detJ==0
 disp('Error cannot continue')
 return
end
Jinv=(inv(J))';
phxy=Jinv*ph';
phxy=phxy';

kappa=lmat;


%for i=1:esize
%for j=1:esize
%for k=1:2
%stiffel(i,j)=stiffel(i,j)+(kappa*intw(p)*detJ*phxy(i,k)*phxy(j,k));
%end
%end

stiffel=stiffel+(kappa*intw(p)*detJ*phxy*phxy');

%end

%end
end

end
