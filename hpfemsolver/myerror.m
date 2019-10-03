function errengy=myerror(sol,sbasish1,nelemq,Meshq,order,unkv,unkvh,unkqih,unkvi,...
glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,basisxth1,basisyth1,...
sbasisth1,nelemt,basisxh1,basisyh1,probdata,coeff,coefft)

% This function compute the relative L2 norm of the error, relative H1 norm and
% energy norm of the error

%extract relevent data	 
coord=Meshq.Coordinates;
intmaq=Meshq.Elements;
intmat=Mesht.Elements;

%intialize values
nip=length(x);
esize=(order+1+1)^2;
errn=0;
errd=0;
errengy=0;
mat=probdata.mat;

fnum=Meshq.Fnum;
for i=1:nelemq

   for l=1:4
      for k=1:2
         xyq(l,k)=coord(intmaq(i,l),k);
      end
   end
   
   for l=1:4
      for k=1:2
      localblend(l,k)=coeff(i,l,k);
      end
   end

% Vertex based dof's  
  for l=1:4
     lunkv(l)=unkv(intmaq(i,l));
     ldrv(l)=1;
  end
  % Higher order edges
  if order >0
  for l=1:4
  for p=0:order-1
  lunkv(4+l+4*(p))=unkvh(glob(i,l),p+1);
  ldrv(4+l+4*(p))=dir(i,l)^(p+2);
  end
  end
  % Interior block
  for l=1:order^2
  lunkv(4+4*order+l)=unkqih(i,l);
  ldrv(4+4*(order)+l)=1;    
  end
  end  
   
   for p=1:nip
   for q=1:nip
      
      [xp,yp]=getcoord(xyq,x(p),x(q));
      
      J=mapquad(xyq,x(p),x(q));
      detJ=det(J);
      Jinvt=inv(J)';
      ph=sbasish1(:,q+(nip*(p-1)));

      v=0;      
      for k=1:esize
      if lunkv(k) ~=0
        v=v+(ph(k).*sol(lunkv(k)).*ldrv(k));      
      end
      end

      phx=basisxh1(:,q+(nip*(p-1)));
      phy=basisyh1(:,q+(nip*(p-1)));
      gph=[phx phy];
      gphxy=Jinvt*gph';
      gphxy=gphxy';

      dv=[0;0];      
      for k=1:esize
      if lunkv(k) ~=0
        dv(1)=dv(1)+(gphxy(k,1).*sol(lunkv(k)).*ldrv(k));      
        dv(2)=dv(2)+(gphxy(k,2).*sol(lunkv(k)).*ldrv(k));  
      end
      end
      
      fun=probdata.exactfun;
      arg=probdata.exactfunarg;
      vex=fun(xp,yp,arg);

      fun=probdata.dexactfun;
      arg=probdata.dexactfunarg;
      dvex=fun(xp,yp,arg);
      
%     perform integration
      errn=errn+(w(p)*w(q)*detJ*(norm(v-vex)^2));
      errd=errd+(w(p)*w(q)*detJ*(norm(vex)^2));

%     perform integration
      errengy=errengy+(w(p)*w(q)*detJ*(norm(dv-dvex)^2))*mat(fnum(i));
   end
   end
end
coord=Mesht.Coordinates;

nipt=length(intxi);
esize=3+3*order+(order-1)*order/2;

fnum=Mesht.Fnum;
for i=1:nelemt

   for l=1:3
      for k=1:2
         xyt(l,k)=coord(intmat(i,l),k);
      end
   end


   for l=1:3
      for k=1:2
      localblendt(l,k)=coefft(i,l,k);
      end
   end
%% Vertex based dof's  
  for l=1:3
     lunkv(l)=unkv(intmat(i,l));
     ldrv(l)=1;
  end
  % Higher order edges
  if order >0
  for l=1:3
  for p=0:order-1
  lunkv(3+l+3*(p))=unkvh(glob(i+nelemq,l),p+1);
  ldrv(3+l+3*(p))=dir(i+nelemq,l)^(p+2);
  end
  end
  % Interior block
  for p=1:order*(order-1)/2
  lunkv(3+3*order+p)=unkvi(i,p);
  ldrv(3+3*(order)+p)=1;    
  end
  end  
   
   for p=1:nipt
      
      [xp,yp]=getcoordt(xyt,intxi(p),inteta(p),localblendt);
      
      J=maptri(xyt,intxi(p),inteta(p),localblendt);
      detJ=det(J);
      Jinvt=inv(J)';
      
      ph=sbasisth1(:,p);
      
      phx=basisxth1(:,p);
      phy=basisyth1(:,p);
      gph=[phx phy];
      phxy=Jinvt*gph';
      phxy=phxy';
      
      
      v=0;      
      for k=1:esize
      if lunkv(k) ~=0
        v=v+(ph(k).*sol(lunkv(k)).*ldrv(k));      
      end
      end
      
      dv=[0;0];      
      for k=1:esize
      if lunkv(k) ~=0
        dv(1)=dv(1)+(phxy(k,1).*sol(lunkv(k)).*ldrv(k));      
        dv(2)=dv(2)+(phxy(k,2).*sol(lunkv(k)).*ldrv(k));      
      end
      end
      
      fun=probdata.exactfun;
      arg=probdata.exactfunarg;
      vex=fun(xp,yp,arg);

      fun=probdata.dexactfun;
      arg=probdata.dexactfunarg; 
      dvex=fun(xp,yp,arg);
      
      
 
%     perform integration
      errn=errn+(intw(p)*detJ*(norm(v-vex)^2));
      errd=errd+(intw(p)*detJ*(norm(vex)^2));

%     perform integration
      errengy=errengy+(intw(p)*detJ*(norm(dv-dvex)^2)*mat(fnum(i)));

   end
end

errn=sqrt(errn);
errengy=sqrt(errengy)

disp(['The L2 norm of the error is ',num2str(errn)]);

disp(['The energy norm of the error is ',num2str(errengy)]);

errn=errn/sqrt(errd);

disp(['The normalised L2 norm of the error is ',num2str(errn)]);
   

%pause
