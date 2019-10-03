function [basisxthdiv,basisythdiv,sdivbasisthdiv,rbasisxthdiv,rbasisythdiv,l2basist,rbasistl2,gl2basistx,gl2basisty]=myevalthdiv(order,intxi,inteta,xt)

nipt=length(intxi);
nip=length(xt);


basisxthdiv=[];
basisythdiv=[];
for p=1:nipt
ph=basisthdiv(order,intxi(p),inteta(p));
basisxthdiv=[basisxthdiv, ph(:,1)];
basisythdiv=[basisythdiv, ph(:,2)];
end

sdivbasisthdiv=[];
for p=1:nipt
sph=divbasisthdiv(order,intxi(p),inteta(p));
sdivbasisthdiv=[sdivbasisthdiv, sph(:)];
end

rbasisxthdiv=[];
rbasisythdiv=[];
for p=1:nip
% edge 1
eta=sqrt(3.)*0.5*(xt(p)+1);
xi=(-1+xt(p))*(-0.5);
ph=basisthdiv(order,xi,eta);
rbasisxthdiv=[rbasisxthdiv, ph(:,1)];
rbasisythdiv=[rbasisythdiv, ph(:,2)];
end
% Edge 2 y=1
for p=1:nip
eta=sqrt(3.)*(-0.5)*(xt(p)-1);
xi=-(1+xt(p))*0.5;
ph=basisthdiv(order,xi,eta);
rbasisxthdiv=[rbasisxthdiv, ph(:,1)];
rbasisythdiv=[rbasisythdiv, ph(:,2)];
end
% Edge 3 x=1
for p=1:nip
eta=0.;
xi=xt(p);
ph=basisthdiv(order,xi,eta);
rbasisxthdiv=[rbasisxthdiv, ph(:,1)];
rbasisythdiv=[rbasisythdiv, ph(:,2)];
end

l2basist=[];
for p=1:nipt
phl2=basistl2(order,intxi(p),inteta(p));
l2basist=[l2basist, phl2(:)];
end


% compute gradient of L2 basis at integration points
gl2basistx=[];
gl2basisty=[];
%for p=1:nipt
%gphl2=gbasistl2(order,intxi(p),inteta(p));
%gl2basistx=[gl2basistx, gphl2(:,1)];
%gl2basisty=[gl2basisty, gphl2(:,2)];
%end




rbasistl2=[];
for p=1:nip
% edge 1
phlag=basislag(xt(p),order);
rbasistl2=[rbasistl2, phlag(:)];
end
% Edge 2 
for p=1:nip
phlag=basislag(xt(p),order);
rbasistl2=[rbasistl2, phlag(:)];
end
% Edge 3 
for p=1:nip
phlag=basislag(xt(p),order);
rbasistl2=[rbasistl2, phlag(:)];
end
