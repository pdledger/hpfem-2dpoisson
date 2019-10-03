function [lbasisxth1,lbasisyth1,lsbasisth1,lrbasisth1,ldivbasisxxth1,ldivbasisxyth1,ldivbasisyxth1,ldivbasisyyth1,lbasisxthdiv,lbasisythdiv,lsdivbasisthdiv,lrbasisxthdiv,lrbasisythdiv,ll2basist,lrbasistl2]=extractbasis(basisxth1,...
basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1,order,maxorder,intw,intxi,inteta,xt,wt...
,basisxthdiv,basisythdiv,sdivbasisthdiv,rbasisxthdiv,rbasisythdiv,l2basist,rbasistl2);

nipt=length(intw);

% set size of local arrays
esizeh1full=3+3*maxorder+(maxorder-1)*maxorder/2;

lesizeh1full=3+3*order+(order-1)*order/2;

for i=1:nipt

ph=sbasisth1(:,i);

phx=basisxth1(:,i);
phy=basisyth1(:,i);

dxx=divbasisxxth1(:,i);
dxy=divbasisxyth1(:,i);

dyx=divbasisyxth1(:,i);
dyy=divbasisyyth1(:,i);


% lowest order dofs
lph=ph(1:3);
lphx=phx(1:3);
lphy=phy(1:3);

ldxx=dxx(1:3);
ldxy=dxy(1:3);
ldyx=dyx(1:3);
ldyy=dyy(1:3);

% edge dofs
for p=0:order-1
for e=1:3
lph(3+3*p+e)=ph(3+3*p+e);
lphx(3+3*p+e)=phx(3+3*p+e);
lphy(3+3*p+e)=phy(3+3*p+e);

ldxx(3+3*p+e)=dxx(3+3*p+e);
ldxy(3+3*p+e)=dxy(3+3*p+e);
ldyx(3+3*p+e)=dyx(3+3*p+e);
ldyy(3+3*p+e)=dyy(3+3*p+e);
end
end

% higher order dofs
basno=3*(maxorder+1);
lbasno=3*(order+1);
for ii=0:maxorder-2
for jj=0:maxorder-2
if ii+jj <= maxorder-2
basno=basno+1;
if ii<= order-2 & jj<=order-2 & ii+jj<=order-2
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
end

lsbasisth1(:,i)=lph;

lbasisxth1(:,i)=lphx;
lbasisyth1(:,i)=lphy;

ldivbasisxxth1(:,i)=ldxx;
ldivbasisxyth1(:,i)=ldxy;

ldivbasisyxth1(:,i)=ldyx;
ldivbasisyyth1(:,i)=ldyy;

end

nip=length(xt);
for j=1:3
for i=1:nip

   ph=rbasisth1(:,i+(nip*(j-1)));
   % lowest order dofs
   lph=ph(1:3);

% edge dofs
for p=0:order-1
for e=1:3
lph(3+3*p+e)=ph(3+3*p+e);
end
end

% higher order dofs
basno=3*(maxorder+1);
lbasno=3*(order+1);
for ii=0:maxorder-2
for jj=0:maxorder-2
if ii+jj <= maxorder-2
basno=basno+1;
if ii<= order-2 & jj<=order-2 & ii+jj<=order-2
lbasno=lbasno+1;
lph(lbasno)=ph(basno);
end
end
end
end

   lrbasisth1(:,i+(nip*(j-1)))=lph;
end
end

% now prepare the basis functions for  H(div)

% set size of local arrays
maxorderh=maxorder+1;
orderh=order+1;
esizehdiv=3+3*maxorderh;
if maxorderh > 0
esizehdivfull=esizehdiv+maxorderh*(maxorderh-1)+maxorderh-1;
else
esizehdivfull=esizehdiv;
end
maxorderl2=maxorderh-1;
if maxorderh==0
esizel2=1;
else
esizel2=1+(maxorderl2*(maxorderl2+1)/2+maxorderl2);
end

esizelag=3*(maxorderh+1);


for i=1:nipt
phx=basisxthdiv(:,i);
phy=basisythdiv(:,i);
sph=sdivbasisthdiv(:,i);

phl2=l2basist(:,i);

lphx=phx(1:3);
lphy=phy(1:3);
lsph=sph(1:3);

lphl2=phl2(1);
% edge dofs
for p=0:orderh-1
for e=1:3
lphx(3+3*p+e)=phx(3+3*p+e);
lphy(3+3*p+e)=phy(3+3*p+e);
lsph(3+3*p+e)=sph(3+3*p+e);
end
end


% higher order dofs type 1
basno=3*(maxorderh+1);
lbasno=3*(orderh+1);
for ii=0:maxorderh-2
for jj=0:maxorderh-2
if ii+jj <= maxorderh-2
basno=basno+1;
if ii<= orderh-2 & jj<=orderh-2 & ii+jj<=orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
lsph(lbasno)=sph(basno);
end
end
end
end
% type 2
basnol2=1;
lbasnol2=1;
for ii=0:maxorderh-2
for jj=0:maxorderh-2
if ii+jj <= maxorderh-2
basno=basno+1;
basnol2=basnol2+1;
if ii<= orderh-2 & jj<=orderh-2 & ii+jj<=orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
lsph(lbasno)=sph(basno);

lbasnol2=lbasnol2+1;
lphl2(lbasnol2)=phl2(basnol2);
end
end
end
end
% type 3
for ii=0:maxorderh-2
basno=basno+1;
basnol2=basnol2+1;
if ii<= orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
lsph(lbasno)=sph(basno);

lbasnol2=lbasnol2+1;
lphl2(lbasnol2)=phl2(basnol2);
end
end



lbasisxthdiv(:,i)=lphx;
lbasisythdiv(:,i)=lphy;
lsdivbasisthdiv(:,i)=lsph;

ll2basist(:,i)=lphl2;

end


for j=1:3
for i=1:nip

   phx=rbasisxthdiv(:,i+(nip*(j-1)));
   phy=rbasisythdiv(:,i+(nip*(j-1)));

   phl2=rbasistl2(:,i+(nip*(j-1)));
   
   
lphx=phx(1:3);
lphy=phy(1:3);

for p=0:orderh
lphl2(p+1)=phl2(p+1);
end

% edge dofs
for p=0:orderh-1
for e=1:3
lphx(3+3*p+e)=phx(3+3*p+e);
lphy(3+3*p+e)=phy(3+3*p+e);
end
end


% higher order dofs type 1
basno=3*(maxorderh+1);
lbasno=3*(orderh+1);
for ii=0:maxorderh-2
for jj=0:maxorderh-2
if ii+jj <= maxorderh-2
basno=basno+1;
if ii<= orderh-2 & jj<=orderh-2 & ii+jj<=orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);
end
end
end
end
% type 2
basnol2=1;
lbasnol2=1;
for ii=0:maxorderh-2
for jj=0:maxorderh-2
if ii+jj <= maxorderh-2
basno=basno+1;
basnol2=basnol2+1;
if ii<= orderh-2 & jj<=orderh-2 & ii+jj<=orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

end
end
end
end
% type 3
for ii=0:maxorderh-2
basno=basno+1;
basnol2=basnol2+1;
if ii<= orderh-2
lbasno=lbasno+1;
lphx(lbasno)=phx(basno);
lphy(lbasno)=phy(basno);

end
end


   lrbasisxthdiv(:,i+(nip*(j-1)))=lphx;
   lrbasisythdiv(:,i+(nip*(j-1)))=lphy;

   lrbasistl2(:,i+(nip*(j-1)))=lphl2;
end
end
