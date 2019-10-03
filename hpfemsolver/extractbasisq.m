function [lbasisxh1,lbasisyh1,lsbasish1,lrbasish1,ldivbasisxxh1,ldivbasisxyh1,ldivbasisyxh1,ldivbasisyyh1,lbasisxqhdiv,lbasisyqhdiv,lsdivbasisqhdiv,lrbasisxqhdiv,lrbasisyqhdiv,ll2basisq,lrbasisql2]=extractbasisq(x,w,...
basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1,order,maxorder,...
basisxqhdiv,basisyqhdiv,sdivbasisqhdiv,rbasisxqhdiv,rbasisyqhdiv,l2basisq,rbasisql2)

nip=length(x);
for p=1:nip
for q=1:nip

phx=basisxh1(:,q+(nip*(p-1)));
phy=basisyh1(:,q+(nip*(p-1)));

ph=sbasish1(:,q+(nip*(p-1)));

dxx=divbasisxxh1(:,q+(nip*(p-1)));
dxy=divbasisxyh1(:,q+(nip*(p-1)));
dyx=divbasisyxh1(:,q+(nip*(p-1)));
dyy=divbasisyyh1(:,q+(nip*(p-1)));


lphx=phx(1:4);
lphy=phy(1:4);
lph=ph(1:4);

ldxx=dxx(1:4);
ldxy=dxy(1:4);
ldyx=dyx(1:4);
ldyy=dyy(1:4);


for pp=1:order
for e=1:4
lph(4*pp+e)=ph(4*pp+e);
lphx(4*pp+e)=phx(4*pp+e);
lphy(4*pp+e)=phy(4*pp+e);

ldxx(4*pp+e)=dxx(4*pp+e);
ldxy(4*pp+e)=dxy(4*pp+e);
ldyx(4*pp+e)=dyx(4*pp+e);
ldyy(4*pp+e)=dyy(4*pp+e);

end
end
basno=4*(maxorder+1);
lbasno=4*(order+1);
for i=0:1:maxorder-1
for j=0:1:maxorder-1
basno=basno+1;
if i<=order-1 & j<= order-1
lbasno=lbasno+1;
lph(lbasno)=ph(basno);
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

ldxx(lbasno)=dxx(basno);
ldxy(lbasno)=dxy(basno);
ldyx(lbasno)=dyx(basno);
ldyy(lbasno)=dyy(basno);

end
end
end

lbasisxh1(:,q+(nip*(p-1)))=lphx;
lbasisyh1(:,q+(nip*(p-1)))=lphy;

lsbasish1(:,q+(nip*(p-1)))=lph;

ldivbasisxxh1(:,q+(nip*(p-1)))=ldxx;
ldivbasisxyh1(:,q+(nip*(p-1)))=ldxy;
ldivbasisyxh1(:,q+(nip*(p-1)))=ldyx;
ldivbasisyyh1(:,q+(nip*(p-1)))=ldyy;


end
end


for j=1:4
for p=1:nip

   ph=rbasish1(:,p+(nip*(j-1)));

lph=ph(1:4);

for pp=1:order
for e=1:4
lph(4*pp+e)=ph(4*pp+e);
end
end

basno=4*(maxorder+1);
lbasno=4*(order+1);

for i=0:1:maxorder-1
for jj=0:1:maxorder-1
basno=basno+1;
if i<=order-1 & jj<= order-1
lbasno=lbasno+1;
lph(lbasno)=ph(basno);
end
end
end

 lrbasish1(:,p+(nip*(j-1)))=lph;

end
end



orderh=order+1;
maxorderh=maxorder+1;


for p=1:nip
for q=1:nip

phx=basisxqhdiv(:,q+(nip*(p-1)));
phy=basisyqhdiv(:,q+(nip*(p-1)));
sph=sdivbasisqhdiv(:,q+(nip*(p-1)));

phl2=l2basisq(:,q+(nip*(p-1)));

lphx=phx(1:4);
lphy=phy(1:4);
lsph=sph(1:4);

lphl2=phl2(1);

for pp=1:orderh
for e=1:4
lsph(4*pp+e)=sph(4*pp+e);
lphx(4*pp+e)=phx(4*pp+e);
lphy(4*pp+e)=phy(4*pp+e);

end
end
basno=4*(maxorderh+1);
lbasno=4*(orderh+1);
%type 1
for i=0:1:maxorderh-1
for j=0:1:maxorderh-1
basno=basno+1;
if i<=orderh-1 & j<= orderh-1
lbasno=lbasno+1;
lsph(lbasno)=sph(basno);
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
end
end
end

%type 2
basnol2=1;
lbasnol2=1;
for i=0:1:maxorderh-1
for j=0:1:maxorderh-1
basno=basno+1;
basnol2=basnol2+1;
if i<=orderh-1 & j<= orderh-1
lbasno=lbasno+1;
lbasnol2=lbasnol2+1;

lsph(lbasno)=sph(basno);
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

lphl2(lbasnol2)=phl2(basnol2);
end
end
end

%type 3
for i=0:1:maxorderh-1
basno=basno+1;
basnol2=basnol2+1;
if i<=orderh-1 
lbasno=lbasno+1;
lsph(lbasno)=sph(basno);
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
lbasnol2=lbasnol2+1;
lphl2(lbasnol2)=phl2(basnol2);


end
basno=basno+1;
basnol2=basnol2+1;

if i<=orderh-1 

lbasno=lbasno+1;
lsph(lbasno)=sph(basno);
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

lbasnol2=lbasnol2+1;
lphl2(lbasnol2)=phl2(basnol2);

end
end


lbasisxqhdiv(:,q+(nip*(p-1)))=lphx;
lbasisyqhdiv(:,q+(nip*(p-1)))=lphy;
lsdivbasisqhdiv(:,q+(nip*(p-1)))=lsph;

ll2basisq(:,q+(nip*(p-1)))=lphl2;

end
end


for jj=1:4
for p=1:nip

phx=rbasisxqhdiv(:,p+(nip*(jj-1)));
phy=rbasisyqhdiv(:,p+(nip*(jj-1)));

phl2=rbasisql2(:,p+(nip*(jj-1)));

lphx=phx(1:4);
lphy=phy(1:4);

for pp=0:orderh
lphl2(pp+1)=phl2(pp+1);
end

for pp=1:orderh
for e=1:4
lphx(4*pp+e)=phx(4*pp+e);
lphy(4*pp+e)=phy(4*pp+e);

end
end
basno=4*(maxorderh+1);
lbasno=4*(orderh+1);
%type 1
for i=0:1:maxorderh-1
for j=0:1:maxorderh-1
basno=basno+1;
if i<=orderh-1 & j<= orderh-1
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
end
end
end

%type 2
basnol2=1;
lbasnol2=1;
for i=0:1:maxorderh-1
for j=0:1:maxorderh-1
basno=basno+1;
if i<=orderh-1 & j<= orderh-1
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
end
end
end


%type 3
for i=0:1:maxorderh-1
basno=basno+1;
if i<=orderh-1 
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

end
basno=basno+1;
if i<=orderh-1 
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

end
end


lrbasisxqhdiv(:,p+(nip*(jj-1)))=lphx;
lrbasisyqhdiv(:,p+(nip*(jj-1)))=lphy;

lrbasisql2(:,p+(nip*(jj-1)))=lphl2;


end
end

