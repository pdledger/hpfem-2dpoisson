function [nunk,nunki,unkv,unkvh,unkvi,unkqih,ndir,helpbc,orderel]=nounk(ne,edges,nelemq,nelemt,Mesht,nboun,bsido,glob,Meshq,order,bctype,orderel)

% create number of vertices, edges and interiors (appropriate to degree of
% elements
% outputs
% nunk: Number of continous DOFS (excluding interiors)
% nunki: Total number of DOFS (including interiors)
% unkv: Global unknown number of each vertex
% unkvh: Global unknown number of higher edge functions
% unkvi: Global unknown number of higher order interior functions (triangle)
% unkqih: Global unknown number of higher order interior functions
% (quadrilaterial)
% ndir: Number of known (Dirichlet values)
% helpbc: Array to assist in application of Dirichlet values for vertex
% functions

intmat=Mesht.Elements;
intmaq=Meshq.Elements;
coord=Mesht.Coordinates;
[npoin dum]=size(coord);

% intialize arrays for higher order numbering
unkh=[];
unkqi=[];
unkti=[];

% set boundary conditions at nodes
% if corner vertex has multiple types, enforce to be type 2.
helpbc=zeros(npoin,1);
for i=1:nboun
    if helpbc(bsido(i,1))==0
        helpbc(bsido(i,1))=bsido(i,4);
    end
    if helpbc(bsido(i,2))==0
        helpbc(bsido(i,2))=bsido(i,4);
    end

end


% fix the order on edges
orderedge=1e6*ones(ne,1);
for i=1:nelemt
    for j=1:3
        if orderedge(glob(i+nelemq,j)) > orderel(i+nelemq)
            orderedge(glob(i+nelemq,j))=orderel(i+nelemq);
        end
    end
end

for i=1:nelemq
    for j=1:4
        if orderedge(glob(i,j)) > orderel(i)
            orderedge(glob(i,j))=orderel(i);
        end
    end
end

% now check the order inside again!
for i=1:nelemt
    maxorderi=0;
    for j=1:3
        if maxorderi < orderedge(glob(i+nelemq,j)) 
            maxorderi=orderedge(glob(i+nelemq,j));
        end
        orderel(i+nelemq)=maxorderi;
    end
end

for i=1:nelemq
    maxorderi=0;
    for j=1:4
        if maxorderi < orderedge(glob(i,j)) 
            maxorderi=orderedge(glob(i,j));
	end
        orderel(i)=maxorderi;
    end
end


% number vertex functions
% note that Dirichlet (known values) are negative
nunk=0;
ndir=0;
for i=1:npoin
    if helpbc(i)==0 | (helpbc(i)~=0 & bctype(helpbc(i))~=2 )
        nunk=nunk+1;
        unkv(i)=nunk;
    else
        ndir=ndir+1;
        unkv(i)=-ndir; % dirichlet nodes
    end
end
disp(['we have',num2str(ndir),'dirchet nodes'])


%Higher order Dof's
unkvh=[];
unkvi=[];
unkqih=[];
nunki=nunk;
if order >=1
    unkvh=zeros(ne,order);
    % Higher order edge functions
    for i=1:ne
        if edges(i,3)==0 | (edges(i,3)~=0 & bctype(edges(i,3)) ~=2 )
            for k=1:orderedge(i)
                nunk=nunk+1;
                unkvh(i,k)=nunk;
            end
        else
            for k=1:orderedge(i)
                ndir=ndir+1;
                unkvh(i,k)=-ndir; % dirichlet edges
            end
        end
    end

    % higher order interior functions (triangles)
    nunki=nunk;
    unkvi=zeros(nelemt,(order-1)*order/2);
    %for i=1:nelemt
    %for j=1:(orderel(i)-1)*orderel(i)/2
    %nunki=nunki+1;
    %unkvi(i,j)=nunki;
    %end
    %end
    % to ensure correct number of interiors!
    
    for ii=1:nelemt
    basno=0;
        for i=0:1:order-2
            for j=0:1:order-2
                if i+j <= order-2
                    basno=basno+1;
                    if i <= orderel(ii+nelemq)-2 & j <= orderel(ii+nelemq)-2 & i+j <= orderel(ii+nelemq)-2
                        nunki=nunki+1;
                        unkvi(ii,basno)=nunki;
                    end
                end
            end
        end
    end


    % higher order interior functions (quadrilaterials)
    unkqih=zeros(nelemq,(order-1)^2);
    for i=1:nelemq
        basno=0;
        for p=0:1:order-1
            for q=0:1:order-1
                basno=basno+1;
                if p <=orderel(i)-1 & q<= orderel(i)-1
                    nunki=nunki+1;
                    unkqih(i,basno)=nunki;
                end
            end
        end
    end

end

disp(['There are',num2str(nunki),' unknowns for this domain']);


