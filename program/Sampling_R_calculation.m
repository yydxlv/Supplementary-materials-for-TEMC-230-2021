function [azimuth,elevation]=Sampling_R_calculation(NO_rotation,XX,YY,ZZ,roll_x,roll_y,roll_z)

[px,py]=size(XX);

XX_R=zeros(px,py);
YY_R=zeros(px,py);
ZZ_R=zeros(px,py);
azimuth=zeros(px,py,NO_rotation);
elevation=zeros(px,py,NO_rotation);

for kk=1:NO_rotation 
Rx=[1 0 0 0;0 cos(roll_x(kk)) -sin(roll_x(kk)) 0;0 sin(roll_x(kk)) cos(roll_x(kk)) 0;0 0 0 1];
Ry=[cos(roll_y(kk)) 0 sin(roll_y(kk)) 0;0 1 0 0;-sin(roll_y(kk)) 0 cos(roll_y(kk)) 0;0 0 0 1];
Rz=[cos(roll_z(kk)) -sin(roll_z(kk)) 0 0;sin(roll_z(kk)) cos(roll_z(kk)) 0 0;0 0 1 0;0 0 0 1];

for i=1:px
    for j=1:py
        TTT=Rz*Ry*Rx*[XX(i,j),YY(i,j),ZZ(i,j),1]';
        XX_R(i,j)=TTT(1);
        YY_R(i,j)=TTT(2);
        ZZ_R(i,j)=TTT(3);
    end
end

[azimuth_tt,elevation_tt,r] = cart2sph(XX_R,YY_R,ZZ_R);
azimuth(:,:,kk)=rad2deg(azimuth_tt);
elevation(:,:,kk)=rad2deg(pi/2-elevation_tt);
end