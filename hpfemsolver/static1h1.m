function [kcc,rc]=static1h1(massel,rhsel,order,esize,cont,stiffel)

if esize~=cont
acc=stiffel(1:cont,1:cont);
aci=stiffel(1:cont,cont+1:esize);
aic=stiffel(cont+1:esize,1:cont);
aii=stiffel(cont+1:esize,cont+1:esize);

lc=rhsel(1:cont);
li=rhsel(cont+1:esize);

kcc=acc-(aci*(aii\aic));
rc=lc-(aci*(aii\li));


else
kcc=stiffel;
rc=rhsel;
end
