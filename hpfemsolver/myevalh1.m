function [basisxh1,basisyh1,sbasish1,rbasish1,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1]=myevalh1(order,x)

nip=length(x);


basisxh1=[];
basisyh1=[];
for p=1:nip
for q=1:nip 
ph=gbasisqh1(order+1,x(p),x(q));
basisxh1=[basisxh1, ph(:,1)];
basisyh1=[basisyh1, ph(:,2)];
end
end

sbasish1=[];
for p=1:nip
for q=1:nip
ph=basisqh1(order+1,x(p),x(q));
sbasish1=[sbasish1, ph(:)];
end
end

divbasisxxh1=[];
divbasisxyh1=[];
divbasisyxh1=[];
divbasisyyh1=[];

for p=1:nip
for q=1:nip
divph=divgbasisqh1(order+1,x(p),x(q));
divbasisxxh1=[divbasisxxh1, divph(:,1)];
divbasisxyh1=[divbasisxyh1, divph(:,2)];
divbasisyxh1=[divbasisyxh1, divph(:,3)];
divbasisyyh1=[divbasisyyh1, divph(:,4)];

end
end

% rbasish1 not required at present.
% 4 bc edge possibilities
% Edge 1 y=0
rbasish1=[];
for p=1:nip
ph=basisqh1(order+1,x(p),0);
rbasish1=[rbasish1, ph(:)];
end
% Edge 2 y=1
for p=1:nip
ph=basisqh1(order+1,x(p),1);
rbasish1=[rbasish1, ph(:)];
end
% Edge 3 x=1
for p=1:nip
ph=basisqh1(order+1,1,x(p));
rbasish1=[rbasish1, ph(:)];
end
% Edge 4 x=0
for p=1:nip
ph=basisqh1(order+1,0,x(p));
rbasish1=[rbasish1, ph(:)];
end




