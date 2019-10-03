function [IDJ11,IDJ12,IDJ21,IDJ22]=dmapquad(xy,x,y)

% compute derivatives of the Inverse Jacobian Matrix
dlam=[-(1-y) -(1-x); (1-y) -x; y x; -y (1-x)];

J=zeros(2);
for i=1:2
for j=1:2
for k=1:4
J(i,j)=J(i,j)+(dlam(k,j)*xy(k,i));
end
end
end

DETJ=det(J);

% Compute derivative of Jacobian Matirx

%J11
J11=J(1,1);
DJ11(1)=0;
DJ11(2)=xy(1,1)-xy(2,1)+xy(3,1)-xy(4,1);

%J12
J12=J(1,2);
DJ12(1)=xy(1,1)-xy(2,1)+xy(3,1)-xy(4,1);
DJ12(2)=0;

%J21
J21=J(2,1);
DJ21(1)=0;
DJ21(2)=xy(1,2)-xy(2,2)+xy(3,2)-xy(4,2);

%J22
J22=J(2,2);
DJ22(1)=xy(1,2)-xy(2,2)+xy(3,2)-xy(4,2);
DJ22(2)=0;

DDETJ(1)=DJ11(1)*J22+J11*DJ22(1)-(DJ12(1)*J21+J12*DJ21(1));
DDETJ(2)=DJ11(2)*J22+J11*DJ22(2)-(DJ12(2)*J21+J12*DJ21(2));

% Note that
%Jinv = 1/detJ [ J22 -J12; -J21 J11]
%Jinvt= 1/detJ [ J22 -J21; -J12 J11]

% Compute derivatives of the inverse transpose Jacobian matrix
%Jinvt_11 = J_22/DET J
IDJ11(1)= (DJ22(1)*DETJ- J22*DDETJ(1))/DETJ^2;
IDJ11(2)= (DJ22(2)*DETJ- J22*DDETJ(2))/DETJ^2;

%Jinvt_12 = -J_21/DET J
IDJ12(1)= -(DJ21(1)*DETJ- J21*DDETJ(1))/DETJ^2;
IDJ12(2)= -(DJ21(2)*DETJ- J21*DDETJ(2))/DETJ^2;

%Jinvt_21 = -J_12/DET J
IDJ21(1)= -(DJ12(1)*DETJ- J12*DDETJ(1))/DETJ^2;
IDJ21(2)= -(DJ12(2)*DETJ- J12*DDETJ(2))/DETJ^2;

%Jinvt_22 = J_11/DET J
IDJ22(1)= (DJ11(1)*DETJ- J11*DDETJ(1))/DETJ^2;
IDJ22(2)= (DJ11(2)*DETJ- J11*DDETJ(2))/DETJ^2;



%IDJ11=zeros(2);
%IDJ12=zeros(2);
%IDJ21=zeros(2);
%IDJ22=zeros(2);
