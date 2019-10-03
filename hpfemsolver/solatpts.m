function [gradsoltx,gradsolty,jump,divsolt,rgradsolt,divsol,rgradsol]=solatpts(sol,sbasish1,nelemq,Meshq,order,unkv,unkvh,unkqih,unkvi,...
glob,nunk,dir,nunki,x,w,Mesht,intxi,inteta,intw,basisxth1,basisyth1,...
sbasisth1,nelemt,basisxh1,basisyh1,xt,wt,ne,helpbc,edges,divbasisxxth1,divbasisxyth1,divbasisyxth1,divbasisyyth1...
,divbasisxxh1,divbasisxyh1,divbasisyxh1,divbasisyyh1,probdata)

% This function compute the solution, its gradient and div (a grad u)
% at integration points inside an element.
% also computes n.a.gradu on boundaries

%extract relevent data	 
coord=Meshq.Coordinates;
[npoin dum]=size(coord);
intmaq=Meshq.Elements;
intmat=Mesht.Elements;

%intialize values
nip=length(x);
esize=(order+1+1)^2;
errn=0;
errd=0;

fnum=Meshq.Fnum;


% Material data
mat=probdata.mat;

disp('Computing data for error estimator');

divsol=[];
for i=1:nelemq

   for l=1:4
      for k=1:2
         xyq(l,k)=coord(intmaq(i,l),k);
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
      phx=basisxh1(:,q+(nip*(p-1)));
      phy=basisyh1(:,q+(nip*(p-1)));
      ph=[phx phy];
      phxy=Jinvt*ph';
      phxy=phxy';

      gradsolqx(i,q+(nip*(p-1)))=0;      
      gradsolqy(i,q+(nip*(p-1)))=0;      

%     Gradient of solution evaluated at integration point
      for k=1:esize
      if lunkv(k) ~=0
      gradsolqx(i,q+(nip*(p-1)))=gradsolqx(i,q+(nip*(p-1)))+(phxy(k,1)*sol(lunkv(k)).*ldrv(k));     
      gradsolqy(i,q+(nip*(p-1)))=gradsolqy(i,q+(nip*(p-1)))+(phxy(k,2)*sol(lunkv(k)).*ldrv(k));     
      end
      end
      
      
      mdivphxx=divbasisxxh1(:,q+(nip*(p-1)));
      mdivphxy=divbasisxyh1(:,q+(nip*(p-1)));
      
      mdivphyx=divbasisyxh1(:,q+(nip*(p-1)));
      mdivphyy=divbasisyyh1(:,q+(nip*(p-1)));
      
      % inverse jacobian derivatives
      [DJ11,DJ12,DJ21,DJ22]=dmapquad(xyq,x(p),x(q));
      dJ11dx=DJ11(1);
      dJ11dy=DJ11(2);
      
      dJ12dx=DJ12(1);
      dJ12dy=DJ12(2);
            
      dJ21dx=DJ21(1);
      dJ21dy=DJ21(2);
      
      dJ22dx=DJ22(1);
      dJ22dy=DJ22(2);
      
      for k=1:esize
      df1dx=dJ11dx*phx(k)+dJ12dx*phy(k)+Jinvt(1,1)*mdivphxx(k)+Jinvt(1,2)*mdivphyx(k);
      df1dy=dJ11dy*phx(k)+dJ12dy*phy(k)+Jinvt(1,1)*mdivphxy(k)+Jinvt(1,2)*mdivphyy(k);
      
      df2dx=dJ21dx*phx(k)+dJ22dx*phy(k)+Jinvt(2,1)*mdivphxx(k)+Jinvt(2,2)*mdivphyx(k);
      df2dy=dJ21dy*phx(k)+dJ22dy*phy(k)+Jinvt(2,1)*mdivphxy(k)+Jinvt(2,2)*mdivphyy(k);
      
      divph(k)=(Jinvt(1,1)*df1dx+Jinvt(1,2)*df1dy)+(Jinvt(2,1)*df2dx+Jinvt(2,2)*df2dy);
      end
      
      
      
      
      divphxy=divph;
      divsol(i,q+(nip*(p-1)))=0;
      for k=1:esize
      if lunkv(k) ~=0
      divsol(i,q+(nip*(p-1)))=divsol(i,q+(nip*(p-1)))+(divphxy(k)*sol(lunkv(k)).*ldrv(k)*mat(fnum(i)));     
      end
      end
 
   end
   end
end




% intialisations
norml=[0 -1; 0 1; 1 0; -1 0]';

% compute gradient on edges
rbasisxh1=[];
rbasisyh1=[];

for j=1:4
    for p=1:nip
        if j==1
            xi=x(p);
            et=0;
        elseif j==2
            xi=x(p);
            et=1;
        elseif j==3
            xi=1;
            et=x(p);
        else
            xi=0;
            et=x(p);
        end

        gph=gbasisqh1(order+1,xi,et);
        rbasisxh1=[rbasisxh1, gph(:,1)];
        rbasisyh1=[rbasisyh1, gph(:,2)];

    end
end


jump=zeros(ne,nip);
rgradsol=[];
bctype=probdata.bctype;
for i=1:nelemq

    for l=1:4
        for k=1:2
            xyq(l,k)=coord(intmaq(i,l),k);
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

    for j=1:4
        for p=1:nip
            if j==1
                xi=x(p);
                et=0;
            elseif j==2
                xi=x(p);
                et=1;
            elseif j==3
                xi=1;
                et=x(p);
            else
                xi=0;
                et=x(p);
            end
            J=mapquad(xyq,xi,et);
            detJ=det(J);
            Jinvt=inv(J)';


            phx=rbasisxh1(:,p+(nip*(j-1)));
            phy=rbasisyh1(:,p+(nip*(j-1)));
            ph=[phx phy];

            phxy=Jinvt*ph';
            phxy=phxy';

            nm=inv(J)'*norml(:,j);
            nm=nm./norm(nm);

            %     Gradient of solution evaluated at integration point
            rgradsolx(i,j,p)=0;
            rgradsoly(i,j,p)=0;
            rgradsol(i,j,p)=0;


            for k=1:esize
                if lunkv(k) ~=0
                    rgradsolx(i,j,p)=rgradsolx(i,j,p)+(phxy(k,1)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));
                    rgradsoly(i,j,p)=rgradsoly(i,j,p)+(phxy(k,2)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));
                    rgradsol(i,j,p)=rgradsol(i,j,p)+(nm(1)*phxy(k,1)+nm(2)*phxy(k,2))*sol(lunkv(k)).*ldrv(k);%*mat(fnum(i));
                end
            end

            % correct this term if on a Neumann edge
            if edges(glob(i,j),3)~=0 & bctype(edges(glob(i,j),3))==1
                [xp,yp]=getcoord(xyq,xi,et);
                fun=probdata.neufun;
                arg=probdata.neufunarg;
                index=edges(glob(i,j),3);
                val=fun(xp,yp,index,nm,arg);
%                gradu=[cos(xp)*sin(yp); sin(xp)*cos(yp)];
%                rgradsol(i,j,p)=rgradsol(i,j,p)-(nm(1)*gradu(1)+nm(2)*gradu(2));
                rgradsol(i,j,p)=rgradsol(i,j,p)-val;
            end


            jump(glob(i,j),p)=jump(glob(i,j),p)+(nm(1)*rgradsolx(i,j,p)+nm(2)*rgradsoly(i,j,p));


        end
    end
end

coord=Mesht.Coordinates;

nipt=length(intxi);
esize=3+3*order+(order-1)*order/2;

fnum=Mesht.Fnum;

% Material data
mat=probdata.mat;

divsolt=[];
gradsoltx=[];
gradsolty=[];
for i=1:nelemt

   for l=1:3
      for k=1:2
         xyt(l,k)=coord(intmat(i,l),k);
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
      
      [xp,yp]=getcoordt(xyt,intxi(p),inteta(p));
      
      J=maptri(xyt,intxi(p),inteta(p));
      detJ=det(J);
      Jinvt=inv(J)';
      
      
      phx=basisxth1(:,p);
      phy=basisyth1(:,p);
      ph=[phx phy];
      phxy=Jinvt*ph';
      phxy=phxy';

      mdivphxx=divbasisxxth1(:,p);
      mdivphxy=divbasisxyth1(:,p);
      
      mdivphyx=divbasisyxth1(:,p);
      mdivphyy=divbasisyyth1(:,p);
      
      
      for k=1:esize
      df1dx=Jinvt(1,1)*mdivphxx(k)+Jinvt(1,2)*mdivphyx(k);
      df1dy=Jinvt(1,1)*mdivphxy(k)+Jinvt(1,2)*mdivphyy(k);
      
      df2dx=Jinvt(2,1)*mdivphxx(k)+Jinvt(2,2)*mdivphyx(k);
      df2dy=Jinvt(2,1)*mdivphxy(k)+Jinvt(2,2)*mdivphyy(k);
      
      divph(k)=(Jinvt(1,1)*df1dx+Jinvt(1,2)*df1dy)+(Jinvt(2,1)*df2dx+Jinvt(2,2)*df2dy);
      end
      
      %for k=1:esize
      %divph(k)=mdivphxx(k)+mdivphyy(k);
      %Jinvt(1,1)*mdivphxx(k)+Jinvt(1,2)*mdivphyx(k)+Jinvt(2,1)*mdivphxy(k)+Jinvt(2,2)*mdivphyy(k);
      %end
      %divphxy=divph/detJ;
      divphxy=divph;
      
      gradsoltx(i,p)=0;      
      gradsolty(i,p)=0;      

%     Gradient of solution evaluated at integration point
      for k=1:esize
      if lunkv(k) ~=0
      gradsoltx(i,p)=gradsoltx(i,p)+(phxy(k,1)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));     
      gradsolty(i,p)=gradsolty(i,p)+(phxy(k,2)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));     
      end
      end
      
      divsolt(i,p)=0;
      for k=1:esize
      if lunkv(k) ~=0
      divsolt(i,p)=divsolt(i,p)+(divphxy(k)*sol(lunkv(k)).*ldrv(k));
      end
      end
      divsolt(i,p)=divsolt(i,p)*mat(fnum(i));
      
      
 
   end
end




% compute normal component of solution on edges

% intialisations
nip=length(xt);
esizeh1=3+3*order+(order-1)*order/2;
norml=[0.5*sqrt(3.) 0.5; -0.5*sqrt(3.) 0.5; 0 -1.]';

% compute gradient on edges
rbasisxth1=[];
rbasisyth1=[];

for j=1:3
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
   rbasisxth1=[rbasisxth1, ph(:,1)];
   rbasisyth1=[rbasisyth1, ph(:,2)];

end
end

rgradsolt=[];
bctype=probdata.bctype;
for i=1:nelemt

   for l=1:3
      for k=1:2
         xyt(l,k)=coord(intmat(i,l),k);
      end
   end

% Vertex based dof's  
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

for j=1:3
% This is a boundary edge and a boundary integral is to be applied   

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
   
   J=maptri(xyt,xi,eta);
   detJ=det(J);
   Jinvt=(inv(J))';
   
   phx=rbasisxth1(:,p+(nip*(j-1)));
   phy=rbasisyth1(:,p+(nip*(j-1)));
   ph=[phx phy];

   phxy=Jinvt*ph';
   phxy=phxy';

   nm=inv(J)'*norml(:,j);
   nm=nm./norm(nm);
         
%     Gradient of solution evaluated at integration point
      rgradsoltx(i,j,p)=0;
      rgradsolty(i,j,p)=0;
      rgradsolt(i,j,p)=0;

      % compute a grad u
      for k=1:esize
      if lunkv(k) ~=0
      rgradsoltx(i,j,p)=rgradsoltx(i,j,p)+(phxy(k,1)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));     
      rgradsolty(i,j,p)=rgradsolty(i,j,p)+(phxy(k,2)*sol(lunkv(k)).*ldrv(k))*mat(fnum(i));     
      rgradsolt(i,j,p)=rgradsolt(i,j,p)+(nm(1)*phxy(k,1)+nm(2)*phxy(k,2))*sol(lunkv(k)).*ldrv(k); %
      end
      end
       rgradsolt(i,j,p)= rgradsolt(i,j,p)*mat(fnum(i));

% correct this term if on a Neumann edge
      if edges(glob(i,j),3)~=0 & bctype(edges(glob(i,j),3))==1
      [xp,yp]=getcoordt(xyt,xi,eta);
%      gradu=[cos(xp)*sin(yp); sin(xp)*cos(yp)];
%      rgradsolt(i,j,p)=rgradsolt(i,j,p)-(nm(1)*gradu(1)+nm(2)*gradu(2));

                fun=probdata.neufun;
                arg=probdata.neufunarg;
                index=edges(glob(i,j),3);
                val=fun(xp,yp,index,nm,arg);
               rgradsolt(i,j,p)=rgradsolt(i,j,p)-val;
%
      end


jump(glob(i,j),p)=jump(glob(i,j),p)+(nm(1)*rgradsoltx(i,j,p)+nm(2)*rgradsolty(i,j,p));


 end
end       
end
disp('done');

% consider vectices

for i=1:1:0%npoin

 if    helpbc(i)~=0
%     bctype(helpbc(i))
 else 
%     helpbc(i)
 end
    testint=0;
    esize=3+3*order+(order-1)*order/2;

    for ii=1:nelemt
        nodes=intmat(ii,1:3);

        flag=0;
        for j=1:3
            if nodes(j)==i
                flag=j;
            end
        end

        for l=1:3
            for k=1:2
                xyt(l,k)=coord(intmat(ii,l),k);
            end
        end
        if flag~=0
            % this element contains this point

            % work out contribution from element
            nip=length(intw);
            for p=1:nip

                [xp,yp]=getcoordt(xyt,intxi(p),inteta(p));
                J=maptri(xyt,intxi(p),inteta(p));
                detJ=det(J);
                Jinvt=inv(J)';
                
                phat(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*intxi(p)-inteta(p));
                phat(2)=inteta(p)/sqrt(3);
                phat(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*intxi(p)-inteta(p));

                gphat(1,1)=1/2;
                gphat(1,2)=-1/(2*sqrt(3));

                gphat(2,1)=0;
                gphat(2,2)=1/sqrt(3);

                gphat(3,1)=-1/2;
                gphat(3,2)=-1/(2*sqrt(3));

                phxyhat=Jinvt*gphat';
                phxyhat=phxyhat';

		 fun=probdata.srcfun;
                 arg=probdata.srcfunarg;
                 val=fun(xp,yp,arg);
		 srho=val;
                testint=testint+detJ*srho*phat(flag)*intw(p);

                %testint=testint-detJ*(phxyhat(flag,1)*gradsoltx(ii,p)+phxyhat(flag,2)*gradsolty(ii,p))*intw(p);

                testint=testint+detJ*phat(flag)*divsolt(ii,p)*intw(p);
                % work out contribution from edge
            end
            nip=length(xt);

            for j=1:3

                
                if edges(glob(ii,j),3)==0 | (edges(glob(ii,j),3)~=0 & bctype(edges(glob(ii,j),3))~=2)
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

                        J=maptri(xyt,xi,eta);
                        detJ=det(J);
                        if detJ==0
                            disp('Error cannot continue')
                            return
                        end
                        Jinv=(inv(J))';
                        norml=[0.5*sqrt(3.) 0.5; -0.5*sqrt(3.) 0.5; 0 -1.]';

                        nm=inv(J)'*norml(:,j);
                        nm=nm./norm(nm);
                        %nm
                        %pause
                        if j~=3
                            dx=xyt(j+1,1)-xyt(j,1);
                            dy=xyt(j+1,2)-xyt(j,2);
                            len=sqrt(dx^2+dy^2);
                        else
                            dx=xyt(1,1)-xyt(j,1);
                            dy=xyt(1,2)-xyt(j,2);
                            len=sqrt(dx^2+dy^2);
                        end
                        %nm=[dy/len; -dx/len];
                        %pause


                        %   %NB diffferent from fortran code!!!
                        if j==1
                            ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*0.5));
                            ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*0.5));
                            v=sqrt((ddgx^2)+(ddgy^2));
                        elseif j==2
                            ddgx=(J(1,1)*(-0.5))+(J(1,2)*(sqrt(3.)*-0.5));
                            ddgy=(J(2,1)*(-0.5))+(J(2,2)*(sqrt(3.)*-0.5));
                            v=sqrt((ddgx^2)+(ddgy^2));
                        else
                            ddgx=J(1,1);
                            ddgy=J(2,1);
                            v=sqrt((ddgx^2)+(ddgy^2));
                        end

                        approx=rgradsolt(ii,j,p);%nm(1)*rgradsoltx(ii,j,p)+nm(2)*rgradsolty(ii,j,p);

                        phat(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*xi-eta);
                        phat(2)=eta/sqrt(3);
                        phat(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*xi-eta);

                        testint=testint-(v*wt(p)*approx*phat(flag)) ;

                    end
                end
            end

        end
    end
    
        esize=(order+1+1)^2;
       for ii=1:nelemq
        nodes=intmaq(ii,1:4);

        flag=0;
        for j=1:4
            if nodes(j)==i
                flag=j;
            end
        end

        for l=1:4
            for k=1:2
                xyq(l,k)=coord(intmaq(ii,l),k);
            end
        end
        if flag~=0
            % this element contains this point

            % work out contribution from element
            nip=length(x);
            for p=1:nip
                for q=1:nip
                [xp,yp]=getcoord(xyq,x(p),x(q));
                J=mapquad(xyq,x(p),x(q));
                detJ=det(J);
                Jinvt=inv(J)';
                srho=2*sin(xp)*sin(yp);
                phat(1)=(1-x(p))*(1-x(q));
                phat(2)=x(p)*(1-x(q));
                phat(3)=x(p)*x(q);
                phat(4)=(1-x(p))*x(q);

                gphat(1,1)=-(1-x(q));
                gphat(1,2)=-(1-x(p));

                gphat(2,1)=(1-x(q));
                gphat(2,2)=-x(p);

                gphat(3,1)=x(q);
                gphat(3,2)=x(p);

                gphat(4,1)=-x(q);
                gphat(4,2)=1-x(p);

                phxyhat=Jinvt*gphat';
                phxyhat=phxyhat';

		 fun=probdata.srcfun;
                 arg=probdata.srcfunarg;
                 val=fun(xp,yp,arg);
		 srho=val;
                testint=testint+detJ*srho*phat(flag)*w(p)*w(q);
 
                %testint=testint-detJ*(phxyhat(flag,1)*gradsolqx(ii,q+(nip*(p-1)))+phxyhat(flag,2)*gradsolqy(ii,q+(nip*(p-1))))*w(p)*w(q);

                testint=testint+detJ*phat(flag)*divsol(ii,q+(nip*(p-1)))*w(p)*w(q);
                % work out contribution from edge
                end
            end
            
            for j=1:4
                if edges(glob(ii,j),3)~=2
                    for p=1:nip

                        if j==1
                            xi=x(p);
                            et=0;
                        elseif j==2
                            xi=x(p);
                            et=1;
                        elseif j==3
                            xi=1;
                            et=x(p);
                        else
                            xi=0;
                            et=x(p);
                        end

                        J=mapquad(xyq,xi,et);
                        detJ=det(J);
                        if detJ==0
                            disp('Error cannot continue')
                            return
                        end
                        Jinv=(inv(J))';
                        norml=[0 -1; 0 1; 1 0; -1 0]';

                        nm=inv(J)'*norml(:,j);
                        nm=nm./norm(nm);
                        
%                         if j~=3
%                             dx=xyt(j+1,1)-xyt(j,1);
%                             dy=xyt(j+1,2)-xyt(j,2);
%                             len=sqrt(dx^2+dy^2);
%                         else
%                             dx=xyt(1,1)-xyt(j,1);
%                             dy=xyt(1,2)-xyt(j,2);
%                             len=sqrt(dx^2+dy^2);
%                         end
%                         nm=[dy/len; -dx/len];
                        %pause


                       
                        if j==1 || j==2
                            v=sqrt(J(1,1)^2+J(2,1)^2);
                        else
                            v=sqrt(J(1,2)^2+J(2,2)^2);
                        end
                        
                        approx=rgradsol(ii,j,p);%nm(1)*rgradsoltx(ii,j,p)+nm(2)*rgradsolty(ii,j,p);
                        phat(1)=(1-xi)*(1-et);
                        phat(2)=xi*(1-et);
                        phat(3)=xi*et;
                        phat(4)=(1-xi)*et;
                        
                        %[xp,yp]=getcoord(xyq,xi,et);
                        %gradu=[cos(xp)*sin(yp); sin(xp)*cos(yp)];
                        %approx=(nm(1)*gradu(1))+nm(2)*gradu(2);
                        %testint=testint+(v*w(p)*approx*phat(flag)) ;
 
                        testint=testint-(v*w(p)*approx*phat(flag)) ;
                        

                    end
                end
            end

        end
    end

    %testint
    %pause
end
