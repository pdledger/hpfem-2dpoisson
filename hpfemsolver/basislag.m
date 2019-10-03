function ph=basislag(x,order)

n=order+1;
ph=zeros(n,1);
for i=0:order
ph(i+1)=leg(x,i);
end
