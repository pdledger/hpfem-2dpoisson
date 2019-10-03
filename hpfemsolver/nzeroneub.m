function [nz,map,irn,icn]=nzeroneub(nelemq,nelemt,order,glob,nunki,Mesht,unkv,unkvh,unkvi,unkqih,Meshq,nunk)

intmat=Mesht.Elements;
intmaq=Meshq.Elements;
coord=Mesht.Coordinates;
[npoin dum]=size(coord);

%create the connectivity table for quadrilaterals
connectq=[];
for i=1:nelemq
data=[];
% Nodal dof's
for j=1:4
data=[data unkv(intmaq(i,j))];
end

% Nodal edges
for j=1:4
for p=1:order
data=[data unkvh(glob(i,j),p)];
end
end
% Interior block
for j=1:order*order
%data=[data unkqih(i,j)];
end
connectq=[connectq; data ];
end


ldataq=4*(order+1);%+order*order; DO NOT INCLUDE INTERIORS
%disp('Found all quadrilateral connections')


%create the connectivity table for triangles
connectt=[];
for i=1:nelemt
data=[];
% Nodal dof's
for j=1:3
data=[data unkv(intmat(i,j))];
end

% Nodal edges
for j=1:3
for p=1:order
data=[data unkvh(glob(i+nelemq,j),p)];
end
end
% Interior block
%for j=1:order*(order-1)/2
%data=[data unkvi(i,j)];
%end

connectt=[connectt; data ];
end
 
% note that ldatat increased by 1 to include possible extra dof for
% circulation
if order >0
ldatat=3+3*order;%+(order*(order-1)/2); %DO NOT INCLUDE INTERIORS
else
ldatat=3;
end
%if order==0
%ldatat=3;
%else
%ldatat=(order+1)*(order+2);
%end
%disp('Found all triangular connections')

% to here !!!


% create the structure ready for building the nonzero enteries
help1=zeros(nunk,1); % help(i)=Number of Unknowns connected to Unknown i
help2=zeros(nunk,max(ldatat,ldataq)); % help(i,j)=List of Unknowns connected to Unknown i
help3=zeros(nunk,1);
for i=1:nelemq
   for j=1:ldataq
      Row=connectq(i,j);
      if Row > 0 
         for k=1:ldataq
            Col=connectq(i,k);
            if Col > 0
% check to see if Columne Node is already stored
               flag=0;
               if( help1(Row) > 0)
                  for p=1:help1(Row)
                     if (help2(Row,p) == Col) 
                        flag=1;
                     end
                  end
               end
               if(flag==0)
                  help1(Row)=help1(Row)+1;
                  help2(Row,help1(Row))=Col;
               end
	    end
         end
      end
   end
end


for i=1:nelemt
   for j=1:ldatat
      Row=connectt(i,j);
      if Row > 0 
         for k=1:ldatat
            Col=connectt(i,k);
            if Col > 0
% check to see if Columne Node is already stored
               flag=0;
               if( help1(Row) > 0)
                  for p=1:help1(Row)
                     if (help2(Row,p) == Col) 
                        flag=1;
                     end
                  end
               end
               if(flag==0)
                  help1(Row)=help1(Row)+1;
                  help2(Row,help1(Row))=Col;
               end
	    end
         end
      end
   end
end

%disp('finished part1')
% create the list of Nonzero enteries in the matrix
nz=0;
map=zeros(nunk+1,1); %map(i)=First non-zero entry on row i
for i=1:nunk
%store the nonzero entry on which each row begins
   map(i)=nz+1;
   for j=1:help1(i)
      nz=nz+1;
      irn(nz)=i;
      icn(nz)=help2(i,j);
   end
end
map(nunk+1)=nz+1;

%disp(['We have found'])
%nz

clear help1
clear help2
clear connectq
clear connectt
