function ddli=ddlegi(x,n)

% evaluate the second derivative of the integrated Legendre polynomial d^2
% (l_n(x))/dx^2
% output:
% ddlin: d^2l_n(x)/dx^2

ddli=(ddleg(x,n)-ddleg(x,n-2))/(2*n-1);
