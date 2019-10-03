function pressurecoef(nelemt,nelemq,order,splitq,dir,glob,Meshq,Mesht,sol,splitt,unkv,xt,wt,...
intxi,inteta,intw,unkvh,unkvi,sbasisth1,edges,probdata,coeff,coefft)

intmaq=Meshq.Elements;
intmat=Mesht.Elements;
coord=Meshq.Coordinates;

[nelem dum]=size(intmat);

esize=3;
if order>0   
esize=(order+1)*(order+2);
end
esizeh1=3+3*order+(order-1)*order/2;
if probdata.lingeom==0
    [xt,wt]=gaussquad(-1.,1,10);
else
    [xt,wt]=gaussquad(-1.,1,1);
end

nip=length(xt);
data=[];

for i=1:nelem

for j=1:3

if edges(glob(i,j),3)==probdata.presscoefbc
% surface of object

       for p=1:3
       for pp=1:2
          xy(p,pp)=coord(intmat(i,p),pp);
       end
       end
       
      for l=1:3
      for k=1:2
      localblendt(l,k)=coefft(i,l,k);
      end
      end
              
      for p=1:nip
      %sam as fortran code...
      if j==1
      eta=sqrt(3.)*0.5*(xt(p)+1);
      xi=(-1+xt(p))*(-0.5);
      elseif j==2
      eta=sqrt(3.)*(-0.5)*(xt(p)-1);
      xi=-1*(1+xt(p))*0.5;
      else
      eta=0.;
      xi=xt(p);
      end
      
      ph=gbasisth1(order+1,xi,eta);
      J=maptri(xy,xi,eta,localblendt);
      Jinv=(inv(J))';
      phxyh1=Jinv*ph';
      phxyh1=phxyh1';


% lowest order block
      for l=1:3
	 lunkv(l)=unkv(intmat(i,l));
	 ldrv(l)=1;
      end
      if order > 0
 % Higher order edges Vertex
         for l=1:3
         for p=0:order-1
            lunkv(3+l+3*(p))=unkvh(glob(i+nelemq,l),p+1);
            ldrv(3+l+3*(p))=dir(i+nelemq,l)^(p+2);
         end
         end
% Interior block Vertex
         for k=1:order*(order-1)/2
            lunkv(3+3*order+k)=unkvi(i,k);
            ldrv(3+3*(order)+k)=1;    
         end
      end
       
      etil=zeros(2,1);
      for k=1:esizeh1
      if lunkv(k) >0
        etil(1)=etil(1)+(phxyh1(k,1)'.*(sol(lunkv(k)).*ldrv(k)));   
        etil(2)=etil(2)+(phxyh1(k,2)'.*(sol(lunkv(k)).*ldrv(k)));   
      end        
      end
      uinf=probdata.vinf;
      [x,y]=getcoordt(xy,xi,eta,localblendt);
      if probdata.presscoefout==1
      coeff=1.-((etil(1)^2+etil(2)^2)/(uinf(1)^2+uinf(2)^2));
      else
      coeff=((etil(1)^2+etil(2)^2)/(uinf(1)^2+uinf(2)^2));
      end
      % compare with
      
      data=[data; x y coeff];
      end
   end
   end
end


[xsdata,ind]=sort(data(:,1));
xdata=[xsdata,data(ind(:),3)];
if probdata.presscoefout==1
figure;
plot(xdata(:,1), xdata(:,2),'r');
title('Pressue Coefficent with x')
xlabel('x')
ylabel('Cp')
[ysdata,ind]=sort(data(:,2));
ydata=[ysdata,data(ind(:),3)];
figure;
plot(ydata(:,1), ydata(:,2),'r');
title('Pressue Coefficent with y')
xlabel('x')
ylabel('Cp')
else
figure;
plot(xdata(:,1), xdata(:,2),'r');
title('Pressue Coefficent with x')
xlabel('x')
ylabel('1-Cp')
end
