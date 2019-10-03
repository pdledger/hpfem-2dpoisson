function [sol,unkv,nunki,unkvh]=streat(sol,dirval,unkv,unkvh,nunki,order)

% thus function places the Diriclet DOFs in the solution vector
% output:
% sol: solution vector updated with Dirichlet values
% unkv: vertex unknowns updated to refer to position of Diriclet values in sol
% nunki: number of unknowns increased to include Dirichlet values
% unkvh: edge unknowns updated to refer to position of Dirichlet values in sol

% update vertex functions
[npoin]=length(unkv);
for i=1:npoin
if unkv(i) < 0
nunki=nunki+1;
sol(nunki)=dirval(abs(unkv(i)));
unkv(i)=nunki;
end
end

% update edge functions
[ne,dum]=size(unkvh);
for i=1:ne
for j=1:order
if unkvh(i,j) < 0
nunki=nunki+1;
sol(nunki)=dirval(abs(unkvh(i,j)));
unkvh(i,j)=nunki;
end
end
end
