function [basisxqhdiv,basisyqhdiv,sdivbasisqhdiv,rbasisxqhdiv,rbasisyqhdiv,l2basisq,rbasisql2]=myevalqhdiv(order,x,w)

nip=length(x);


basisxqhdiv=[];
basisyqhdiv=[];
for p=1:nip
for q=1:nip
ph=basisqhdiv(order,x(p),x(q));
basisxqhdiv=[basisxqhdiv, ph(:,1)];
basisyqhdiv=[basisyqhdiv, ph(:,2)];

end
end

sdivbasisqhdiv=[];
for p=1:nip
    for q=1:nip
sph=divbasisqhdiv(order,x(p),x(q));
sdivbasisqhdiv=[sdivbasisqhdiv, sph(:)];
    end
end

rbasisxqhdiv=[];
rbasisyqhdiv=[];
for p=1:nip
% edge 1
ph=basisqhdiv(order,x(p),0);
rbasisxqhdiv=[rbasisxqhdiv, ph(:,1)];
rbasisyqhdiv=[rbasisyqhdiv, ph(:,2)];
end
% Edge 2 y=1
for p=1:nip
ph=basisqhdiv(order,x(p),1);
rbasisxqhdiv=[rbasisxqhdiv, ph(:,1)];
rbasisyqhdiv=[rbasisyqhdiv, ph(:,2)];
end
% Edge 3 x=1
for p=1:nip
ph=basisqhdiv(order,1,x(p));
rbasisxqhdiv=[rbasisxqhdiv, ph(:,1)];
rbasisyqhdiv=[rbasisyqhdiv, ph(:,2)];
end
% Edge 3 x=1
for p=1:nip
ph=basisqhdiv(order,0,x(p));
rbasisxqhdiv=[rbasisxqhdiv, ph(:,1)];
rbasisyqhdiv=[rbasisyqhdiv, ph(:,2)];
end



l2basisq=[];
for p=1:nip
for q=1:nip
phl2=basisql2(order,x(p),x(q));
l2basisq=[l2basisq, phl2(:)];
end
end


rbasisql2=[];
for p=1:nip
% edge 1
phlag=basislag(2*x(p)-1,order);
rbasisql2=[rbasisql2, phlag(:)];
end
% Edge 2 y=1
for p=1:nip
phlag=basislag(1-2*x(p),order);
rbasisql2=[rbasisql2, phlag(:)];
end
% Edge 3 x=1
for p=1:nip
phlag=basislag(2*x(p)-1,order);
rbasisql2=[rbasisql2, phlag(:)];
end
% Edge 3 x=0
for p=1:nip
phlag=basislag(1-2*x(p),order);
rbasisql2=[rbasisql2, phlag(:)];
end


