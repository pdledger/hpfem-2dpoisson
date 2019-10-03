function ph=divgbasisth1(order,x,y);

ph(1,1:4)=0;
ph(2,1:4)=0;
ph(3,1:4)=0;

if order > 1
    s=[(1/2*y*sqrt(3))-1/2-(1/2*x);
        1/2-(1/2*x)-1/2*y*sqrt(3);
        x];
    t=[1/2+1/2*x+1/6*y*sqrt(3);
        1/2-1/2*x+1/6*y*sqrt(3);
        1-1/3*y*sqrt(3)];

    ds=[-1/2 1/2*sqrt(3);
        -1/2 -1/2*sqrt(3);
        1 0];

    dt=[1/2 1/6*sqrt(3);
        -1/2 1/6*sqrt(3);
        0 -1/3*sqrt(3)];
    %for i=1:3
    %  dx(i)=(t(i)*ds(i,1)-s(i)*dt(i,1))/t(i)^2;
    %  dy(i)=(t(i)*ds(i,2)-s(i)*dt(i,2))/t(i)^2;
    %end

    for p=0:order-2
        for e=1:3

            if t(e)==0
                a=0; %the limit of ls(s,t,n) =0 for n>=2 but a=t(e)*ls(s,t,n) for n>=0 and so a=0 for all n when t(e)=0
                if p>=1
                    b=0;  % limit of ls(s,t,n) =0 for n >=2
                else
                    b=s(e)/2 % limit only true for n>=2 but for n=1 ls(s,t,1)=s
                end
            else
                a=t(e)^(p+1)*leg(s(e)/t(e),p);
                b=t(e)^(p+1)*leg(s(e)/t(e),p+1);

                dargdx=(ds(e,1)*t(e)-s(e)*dt(e,1))/t(e)^2;
                dadx=(p+1)*t(e)^p*dt(e,1)*leg(s(e)/t(e),p)+t(e)^(p+1)*dleg(s(e)/t(e),p)*dargdx;

                dargdy=(ds(e,2)*t(e)-s(e)*dt(e,2))/t(e)^2;
                dady=(p+1)*t(e)^p*dt(e,2)*leg(s(e)/t(e),p)+t(e)^(p+1)*dleg(s(e)/t(e),p)*dargdy;

                dbdx=(p+1)*t(e)^p*dt(e,1)*leg(s(e)/t(e),p+1)+t(e)^(p+1)*dleg(s(e)/t(e),p+1)*dargdx;
                dbdy=(p+1)*t(e)^p*dt(e,2)*leg(s(e)/t(e),p+1)+t(e)^(p+1)*dleg(s(e)/t(e),p+1)*dargdy;


            end

%            ph(3*(p+1)+e)=-dadx*dt(e,1)+dbdx*ds(e,1)-dady*dt(e,2)+dbdy*ds(e,2);

            ph(3*(p+1)+e,1)=-dadx*dt(e,1)+dbdx*ds(e,1);
            ph(3*(p+1)+e,2)=-dady*dt(e,1)+dbdy*ds(e,1);

            ph(3*(p+1)+e,3)=-dadx*dt(e,2)+dbdx*ds(e,2);


	    ph(3*(p+1)+e,4)=-dady*dt(e,2)+dbdy*ds(e,2);
	    



            %p*t(e)^(p-1)*dt(e,1)*legi(s(e)/t(e),p)+t(e)^p*dlegi(s(e)/t(e),p)*dx(e);
            %ph(3*(p+1)+e,2)=p*t(e)^(p-1)*dt(e,2)*legi(s(e)/t(e),p)+t(e)^p*dlegi(s(e)/t(e),p)*dy(e);
        end
    end
end


    edgefun=3*(order);

    % Interior Functions
    if order > 1
        basno=edgefun;
        l(1)=1/(2*sqrt(3))*(sqrt(3)+sqrt(3)*x-y);
        l(2)=y/sqrt(3);
        l(3)=1/(2*sqrt(3))*(sqrt(3)-sqrt(3)*x-y);
        dl=[1/2 -1/(2*sqrt(3));
            0 1/sqrt(3);
            -1/2 -1/(2*sqrt(3))];

        s=l(2)-l(1);
        t=l(1)+l(2);
        ds=[dl(2,1)-dl(1,1), dl(2,2)-dl(1,2)];
        dt=[dl(2,1)+dl(1,1), dl(2,2)+dl(1,2)];

        % Type 1
        for i=0:1:order-3
            for j=0:1:order-3
                if i+j <= order-3
                    basno=basno+1;

                    if t==0
                        if i>= 1
                            duix=0;
                            duiy=0;
                            ui=0;
                        else
                            duix=s*ds(1);
                            duiy=s*ds(2);
                            ui=s^2/2;
                        end
                    else
                        ui=legi(s/t,i+2)*t^(i+2);
                        duix=-leg(s/t,i)*t^i*t*dt(1)+leg(s/t,i+1)*t^(i+1)*ds(1);
                        duiy=-leg(s/t,i)*t^i*t*dt(2)+leg(s/t,i+1)*t^(i+1)*ds(2);

                        dargx=(ds(1)*t-s*dt(1))/t^2;
                        dduixx=-dleg(s/t,i)*dargx*t^i*t*dt(1)-leg(s/t,i)*(i+1)*t^i*dt(1)*dt(1)+dleg(s/t,i+1)*dargx*t^(i+1)*ds(1)+leg(s/t,i+1)*(i+1)*t^(i)*dt(1)*ds(1);

                        dargy=(ds(2)*t-s*dt(2))/t^2;
                        dduixy=-dleg(s/t,i)*dargy*t^i*t*dt(1)-leg(s/t,i)*(i+1)*t^i*dt(1)*dt(2)+dleg(s/t,i+1)*dargy*t^(i+1)*ds(1)+leg(s/t,i+1)*(i+1)*t^(i)*dt(2)*ds(1);


                        dduiyy=-dleg(s/t,i)*dargy*t^i*t*dt(2)-leg(s/t,i)*(i+1)*t^i*dt(2)*dt(2)+dleg(s/t,i+1)*dargy*t^(i+1)*ds(2)+leg(s/t,i+1)*(i+1)*t^(i)*dt(2)*ds(2);

                        dduiyx=-dleg(s/t,i)*dargx*t^i*t*dt(2)-leg(s/t,i)*(i+1)*t^i*dt(1)*dt(2)+dleg(s/t,i+1)*dargx*t^(i+1)*ds(2)+leg(s/t,i+1)*(i+1)*t^(i)*dt(1)*ds(2);

                    end

                    vj=l(3)*leg(l(3)-l(2)-l(1),j);
                    dvjx=dl(3,1)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1));
                    dvjy=dl(3,2)*leg(l(3)-l(2)-l(1),j)+l(3)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2));

                    dargx=dl(3,1)-dl(2,1)-dl(1,1);

                    dargy=dl(3,2)-dl(2,2)-dl(1,2);

                    ddvjxx=dl(3,1)*dleg(l(3)-l(2)-l(1),j)*dargx+dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))+l(3)*ddleg(l(3)-l(2)-l(1),j)*dargx^2;

                   
                    ddvjxy=dl(3,1)*dleg(l(3)-l(2)-l(1),j)*dargy+dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,1)-dl(2,1)-dl(1,1))+l(3)*ddleg(l(3)-l(2)-l(1),j)*dargx*dargy;


                    ddvjyy=dl(3,2)*dleg(l(3)-l(2)-l(1),j)*dargy+dl(3,2)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))+l(3)*ddleg(l(3)-l(2)-l(1),j)*dargy^2;
		    
		    ddvjyx=dl(3,2)*dleg(l(3)-l(2)-l(1),j)*dargx+dl(3,1)*dleg(l(3)-l(2)-l(1),j)*(dl(3,2)-dl(2,2)-dl(1,2))+l(3)*ddleg(l(3)-l(2)-l(1),j)*dargy*dargx;


                    %ph(basno,1)=vj*duix+ui*dvjx;
                    %ph(basno,2)=vj*duiy+ui*dvjy;

                    ph(basno,1)=dvjx*duix+vj*dduixx+duix*dvjx+ui*ddvjxx;
		    
		    ph(basno,2)=dvjy*duix+vj*dduixy+duiy*dvjx+ui*ddvjxy;

                    ph(basno,3)=dvjx*duiy+vj*dduiyx+duix*dvjy+ui*ddvjyx;	    
		    
                    ph(basno,4)=dvjy*duiy+vj*dduiyy+duiy*dvjy+ui*ddvjyy;
                end
            end
        end

    end
end
