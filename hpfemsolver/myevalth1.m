function [basisxth1,basisyth1,sbasisth1,rbasisth1,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1]=myevalth1(order,intxi,inteta,xt)

nipt=length(intxi);
nip=length(xt);


basisxth1=[];
basisyth1=[];
for p=1:nipt
ph=gbasisth1(order+1,intxi(p),inteta(p));
basisxth1=[basisxth1, ph(:,1)];
basisyth1=[basisyth1, ph(:,2)];
end

sbasisth1=[];
for p=1:nipt
ph=basisth1(order+1,intxi(p),inteta(p));
sbasisth1=[sbasisth1, ph(:)];
end

divbasisxxth1=[];
divbasisxyth1=[];
divbasisyxth1=[];
divbasisyyth1=[];

for p=1:nipt
divph=divgbasisth1(order+1,intxi(p),inteta(p));
divbasisxxth1=[divbasisxxth1, divph(:,1)];
divbasisxyth1=[divbasisxyth1, divph(:,2)];
divbasisyxth1=[divbasisyxth1, divph(:,3)];
divbasisyyth1=[divbasisyyth1, divph(:,4)];

end




rbasisth1=[];
for p=1:nip
% edge 1
eta=sqrt(3.)*0.5*(xt(p)+1);
xi=(-1+xt(p))*(-0.5);
ph=basisth1(order+1,xi,eta);
rbasisth1=[rbasisth1, ph(:)];
end
% Edge 2 y=1
for p=1:nip
eta=sqrt(3.)*(-0.5)*(xt(p)-1);
xi=-(1+xt(p))*0.5;
ph=basisth1(order+1,xi,eta);
rbasisth1=[rbasisth1, ph(:)];
end
% Edge 3 x=1
for p=1:nip
eta=0.;
xi=xt(p);
ph=basisth1(order+1,xi,eta);
rbasisth1=[rbasisth1, ph(:)];
end



