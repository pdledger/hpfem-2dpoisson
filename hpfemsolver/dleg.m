function dln=dleg(x,n)

% evaluate the derivative of the Legendre polynomial d (L_n(x))/dx
% output:
% dln: dL_n(x)/dx

if n <=0  
dln=0.;
else
a=((-2.*n.*n.*x.*leg(x,n))+(2.*n.*n.*leg(x,n-1)));
b=(2.*n.*(1-(x.^2)));
if a==0 && b==0
dln=0; %0???
tol=1e-12;
% cannot evaluate derivatives at end points using above formula!
if abs(x+1) < 1e-9
dln=dleg(x+tol,n);
else
dln=dleg(x-tol,n);
end


elseif a~=0 && b==0
disp('error in gettting dln');
else
dln=a./b;
end
end
