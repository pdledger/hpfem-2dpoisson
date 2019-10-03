function [bsido] = myerrorboundaries(p,t)
%scatterer,bcscatter)

% Build data structures used in the solver.
%
% This function builds the mesh parameters and connectivity that will be
% used in tvd_rk2.m. Building this initially saves HEAPS of CPU time later 
% in tvd_rk2.m, but of course, my meshes are fixed for the integration.
%
% Darren Engwirda - 2005-2006.
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.

numn = size(p,1);
numt = size(t,1);
vect = 1:numt;

% DETERMINE UNIQUE EDGES IN MESH
 
e       = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];             % Edges - not unique
vec     = (1:size(e,1))';                                   % List of edge numbers
[e,i,j] = unique(sort(e,2),'rows');                         % Unique edges
vec     = vec(j);                                           % Unique edge numbers
eINt    = [vec(vect), vec(vect+numt), vec(vect+2*numt)];    % Unique edges in each triangle



% DETERMINE EDGE TO TRIANGLE CONNECTIVITY

% Each row has two entries corresponding to the triangle numbers
% associated with each edge. Boundary edges have one entry = 0.
nume = size(e,1);
e2t  = repmat(0,nume,2);
ndx  = repmat(1,nume,1);
for k = 1:numt
    % Edge in kth triangle
    e1 = eINt(k,1); e2 = eINt(k,2); e3 = eINt(k,3);
    % Edge 1
    e2t(e1,ndx(e1)) = k; ndx(e1) = ndx(e1)+1;
    % Edge 2
    e2t(e2,ndx(e2)) = k; ndx(e2) = ndx(e2)+1;
    % Edge 3
    e2t(e3,ndx(e3)) = k; ndx(e3) = ndx(e3)+1;
end



% DETERMINE NODE TO EDGE CONNECTIVITY

% Determine maximum neighbours
j = e(:);
v = repmat(0,max(j),1);
for k = 1:length(j)
    jk = j(k); v(jk) = v(jk)+1;
end
maxN = max(v);

n2e = repmat(0,numn,maxN+1);
ndx = repmat(1,numn,1);
for k = 1:nume
    % End nodes
    n1 = e(k,1); n2 = e(k,2);
    % Connectivity
    n2e(n1,ndx(n1)) = k; ndx(n1) = ndx(n1)+1;
    n2e(n2,ndx(n2)) = k; ndx(n2) = ndx(n2)+1;
end


% DETERMINE NODE TO NODE CONNECTIVITY

n2n = repmat(0,numn,maxN+1);
for k = 1:numn
    next = 1; m = 1;
    while n2e(k,m)>0
        if e(n2e(k,m),1)==k
            n2n(k,next) = e(n2e(k,m),2); next = next+1;
        else
            n2n(k,next) = e(n2e(k,m),1); next = next+1; 
        end
        m = m+1;
    end
end



% FLAG BOUNDARY ELEMENTS

vec     = (1:nume)';            % Edge list
be      = vec(~all(e2t,2));     % Boundary edges
bn      = e(be,:);          
bno     = unique(bn(:));        % Boundary nodes
bnd     = false(numn,1);        
bnd(bn) = true;                 % True for boundary nodes

[nbe dum]=size(bn);

disp(['we have',num2str(nbe),'boundary edges'])

for i=1:nbe
bsido(i,1)=bn(i,1);
bsido(i,2)=bn(i,2);
bsido(i,3)=bn(i,1);
bsido(i,4)=3;
bsido(i,5)=3;
end

