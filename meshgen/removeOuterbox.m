function [pNew,tNew] = removeOuterbox(p,t)

nnodes = size(p,1);
nelem = size(t,1);
help = zeros(nnodes, 1);
help2 = zeros(nnodes,1);
nnodesNew = 0;
pNew = zeros(nnodes-4,2);
tNew = zeros(nelem,3);

for i = 1:nelem
    for j = 1:3
        help(t(i,j)) = 1;
    end
end

for i = 1:nnodes
    if help(i) == 1
        nnodesNew = nnodesNew + 1;
        help2(i) = nnodesNew;
        pNew(nnodesNew,:) = p(i,:);
    end
end

for i = 1:nelem
    for j = 1:3
        tNew(i,j) = help2(t(i,j));
    end
end