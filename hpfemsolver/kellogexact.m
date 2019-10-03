function kellog=kellogexact(x,y)
R=161.4476387975881;
tau=0.1;
rho=pi/4;
sigma=-14.92256510455152;

r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
 if abs(atan2(y,x)-pi) < 1e-10 
 theta=pi;
 elseif abs(atan2(y,x)+pi) < 1e-10
 theta=pi;
 end

if theta < 0
theta=theta+2*pi;
end
if theta <= pi/2
mu=cos((pi/2-sigma)*tau)*cos((theta-pi/2+rho)*tau);
elseif theta <= pi
		mu=cos(rho*tau)*cos((theta-pi+sigma)*tau);
elseif theta <= 3*pi/2
			mu=cos(sigma*tau)*cos((theta-pi-rho)*tau);
else
			mu=cos((pi/2-rho)*tau)*cos((theta-3*pi/2-sigma)*tau);
end
kellog=(r^tau)*mu;
