function ph=divbasisqhdiv(order,x,y);
esizehdivfull=4*(order+1)+2*order*(order+1);

ph=zeros(esizehdivfull,1);

% constant function
sigdiff=[2*x-1; 1-2*x; 2*y-1; 1-2*y];
lamplus=[1-y; y; x ;1-x];

glamplus=[0 -1; 0 1; 1 0; -1 0];

csigdiff=[0 -2; 0 2; 2 0; -2 0];

for m=1:4
ph(m)=1/2*(glamplus(m,1)*csigdiff(m,1)+glamplus(m,2)*csigdiff(m,2));
end


% divergence of higher order edge based functions vanish
for p=1:order
for e=1:4
ph(4*p+e)=0;
end
end    

edgefun=4*(order+1);

if order > 0
basno=edgefun;

% Type 1 vanish
% Type 1
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno)=0;
end
end


% Type 2
for i=0:1:order-1
for j=0:1:order-1
basno=basno+1;
ph(basno)=dlegi(2*y-1,j+2)*dlegi(2*x-1,i+2)*2+dlegi(2*x-1,i+2)*dlegi(2*y-1,j+2)*2;
end
end

% Type 3
for i=0:1:order-1
basno=basno+1;
ph(basno)=2*dlegi(2*y-1,i+2);
basno=basno+1;
ph(basno)=2*dlegi(2*x-1,i+2);
end


end

[m,n]=size(ph);
if m~=4*(order+1)+2*order*(order+1)
  disp('error wrong number of basis functions')
end
  
