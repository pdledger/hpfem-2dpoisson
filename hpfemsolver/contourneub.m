function contourneub(nelemt,nelemq,orders,Meshst,Meshsq,order,splitq,dir,glob,Meshq,Mesht,sol,splitt,unkv,...
unkvh,unkvi,unkqih,edges,coeff,coefft);

% this function plots the solution on the split mesh.
% outputs:
% plots of the contours of the solution and its gradient

% extract relevent data from structures
coordn=Meshsq.Coordinates;
intman=Meshst.Elements;
[npoinn,i]=size(coordn);
intmaq=Meshq.Elements;
intmat=Mesht.Elements;
coord=Meshq.Coordinates;

% Quantities for triangular splittings
nintpt=0;
for i=1:orders
   nintpt=nintpt+i;
end
%  number of elements per triangle
nelpt=4;
for i=1:orders
   nelpt=nelpt+3+(2*i);
end
% number of points per triangle
nptri=6;
for i=1:orders
   nptri=nptri+(3+i);
end

% determine the location of the nodes on a split reference element
epetq=splitquad(orders);
[epett,layer]=splittri(orders,nptri);

% evakyate the basis functions at the nodes on the split elements
phxh1store=[];
phyh1store=[];
phh1store=[];
for q=1:(orders+3)^2            
   xi=epetq(q,1);
   eta=epetq(q,2);
   ph=gbasisqh1(order+1,xi,eta);
   phxh1store=[phxh1store, ph(:,1)];
   phyh1store=[phyh1store, ph(:,2)];
   phh1=basisqh1(order+1,xi,eta);
   phh1store=[phh1store, phh1];
end
phxth1store=[];
phyth1store=[];
phth1store=[];
for q=1:nptri          
   xi=0.98*epett(q,1);
   eta=0.98*epett(q,2);
   ph=gbasisth1(order+1,xi,eta);
   phxth1store=[phxth1store, ph(:,1)];
   phyth1store=[phyth1store, ph(:,2)];
   phh1=basisth1(order+1,xi,eta);
   phth1store=[phth1store, phh1];
end

out=zeros(npoinn,2);
help=zeros(npoinn,1);
outc=zeros(npoinn,1);
outh1=zeros(npoinn,1);

% loop over quad elements...      
for j=1:nelemq     
   esize=4*(order+1)+2*order*(order+1);
   esizeh1=4+4*order+order*order;

   for q=1:(orders+3)^2
       for p=1:4
       for pp=1:2
          xy(p,pp)=coord(intmaq(j,p),pp);
       end
       end
       
      for p=1:4
      for pp=1:2
      localblend(p,pp)=coeff(j,p,pp);
      end
      end

            
       xi=epetq(q,1);
       eta=epetq(q,2);
       J=mapquad(xy,xi,eta);
       detJ=det(J);
       Jinvt=(inv(J))';
       phx=phxh1store(:,q);
       phy=phyh1store(:,q);
       ph=[phx phy];
       phxyh1=Jinvt*ph';
       phxyh1=phxyh1';
       phh1=phh1store(:,q);
% lowest order block
      for l=1:4
	 lunkv(l)=unkv(intmaq(j,l));
	 ldrv(l)=1;
      end
      if order > 0
 % Higher order edges Vertex
         for l=1:4
         for p=0:order-1
            lunkv(4+l+4*(p))=unkvh(glob(j,l),p+1);
            ldrv(4+l+4*(p))=dir(j,l)^(p+2);
         end
         end
% Interior block Vertex
         for k=1:order*order
            lunkv(4+4*order+k)=unkqih(j,k);
            ldrv(4+4*(order)+k)=1;    
         end
       end
      

%     determine the solution and the gradient of the solution at this node      
      etil=zeros(2,1);
      phiv=0;
      for k=1:esizeh1
      if lunkv(k) >0
        etil(1)=etil(1)+(phxyh1(k,1)'.*(sol(lunkv(k)).*ldrv(k)));   
        etil(2)=etil(2)+(phxyh1(k,2)'.*(sol(lunkv(k)).*ldrv(k)));   
      end        
      end
      for k=1:esizeh1
      if lunkv(k) >0
        phiv=phiv+(phh1(k)'.*(sol(lunkv(k)).*ldrv(k)));
      end        
      end

%     store these values at relevent nodes in the split mesh      
      out(splitq(j,q),1)=out(splitq(j,q),1)+(etil(1));
      out(splitq(j,q),2)=out(splitq(j,q),2)+(etil(2));
      outh1(splitq(j,q),1)=outh1(splitq(j,q),1)+(phiv);      
      help(splitq(j,q))=help(splitq(j,q))+1;
   end
end

% loop over triangular elements...      
for j=1:nelemt   
   esize=3;
   if order>0   
   esize=(order+1)*(order+2);
   end
   esizeh1=3+3*order+(order-1)*order/2;

       for p=1:3
       for pp=1:2
          xyt(p,pp)=coord(intmat(j,p),pp);
       end
       end
       
       for p=1:3
       for pp=1:2
       localblendt(p,pp)=coefft(i,p,pp);
       end
       end

       
       for p=1:3
       lbc(p)=edges(glob(j,p),3);
       end

 
   for q=1:nptri  

            
       xi=epett(q,1);
       eta=epett(q,2);
       [xp,yp]=getcoordt(xyt,xi,eta,localblendt);
       J=maptri(xyt,xi,eta,localblendt);
       detJ=det(J);
       Jinvt=(inv(J))';
       
       phx=phxth1store(:,q);
       phy=phyth1store(:,q);
       ph=[phx phy];
       phxyh1=Jinvt*ph';
       phxyh1=phxyh1';
       phh1=phth1store(:,q);
       
% lowest order block
      for l=1:3
	   lunkv(l)=unkv(intmat(j,l));
	   ldrv(l)=1;
      end
      if order > 0
 % Higher order edges Vertex
         for l=1:3
         for p=0:order-1
            lunkv(3+l+3*(p))=unkvh(glob(j+nelemq,l),p+1);
            ldrv(3+l+3*(p))=dir(j+nelemq,l)^(p+2);
         end
         end
% Interior block Vertex
         for k=1:order*(order-1)/2
           lunkv(3+3*order+k)=unkvi(j,k);
           ldrv(3+3*(order)+k)=1;    
         end
      end 
      
%     determine the solution and the gradient of the solution at this node      
      etil=zeros(2,1);
      phiv=0;
      for k=1:esizeh1
      if lunkv(k) >0
        etil(1)=etil(1)+(phxyh1(k,1)'.*(sol(lunkv(k)).*ldrv(k)));   
        etil(2)=etil(2)+(phxyh1(k,2)'.*(sol(lunkv(k)).*ldrv(k)));   
      end        
      end
      for k=1:esizeh1
      if lunkv(k) >0
        phiv=phiv+(phh1(k)'.*(sol(lunkv(k)).*ldrv(k)));
      end        
      end
      
       rp=sqrt(xp^2+yp^2);
 thvertaxis=pi/2-atan2(yp,xp);

% phiv=rp^(2/3)*sin((2*thvertaxis+pi)/3);
% phiv=rp^(2/3)*(sin(2*thvertaxis/3)*cos(pi/3)+cos(2*thvertaxis/3)*sin(pi/3));     
      out(splitt(j,q),1)=out(splitt(j,q),1)+(etil(1));
      out(splitt(j,q),2)=out(splitt(j,q),2)+(etil(2));
      outh1(splitt(j,q))=outh1(splitt(j,q))+phiv;
      help(splitt(j,q))=help(splitt(j,q))+1;
   end
end


for i=1:npoinn
if help(i) ~= 0
%minus since E=-grad phi
out(i,1)=out(i,1)./help(i);
out(i,2)=out(i,2)./help(i);
outh1(i)=outh1(i)/help(i);
else 
out(i,1)=0;
out(i,2)=0;
disp('point not found')
end
end

% plot out the solution
figure
plot_LFE(outh1,Meshst);
plot_LFE(outh1,Meshsq);
colorbar
title('plot of the  potential')
figure
quiver(coordn(:,1),coordn(:,2),out(:,1),out(:,2))
title('quiver plot of the gradient of the potential')


% plot out the gradient of the solution
for i=1:2
figure
hold on;
plot_LFE(real(out(:,i)),Meshsq);
plot_LFE(real(out(:,i)),Meshst);
colorbar
hold off;
title(['plot of grad phi(',num2str(i),')']);
end

