function dli=dlegi(x,n)

% evaluate the derivative of the integrated Legendre polynomial d (l_n(x))/dx
% output:
% dlin: dl_n(x)/dx

dli=(dleg(x,n)-dleg(x,n-2))/(2*n-1);
