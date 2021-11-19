clear
clc
%% composite beam peak search
phi_3dB_element=260;
theta_3dB_element=130; 
Am=30; 
SLAv=30;
Ge_max=1.5; 
NO_rotation=50000;

N_H_sweep = [4,8]; 

for k = 1:2
    N_H=N_H_sweep(k);
    N_V=2; 
    d_H=0.5; 
    d_V=0.5;
    stepsize_sweep = [2.5,3,3.6,4,4.5,5,6,7.5,9,10,12,15,18,20,30];
    for g = 1:15
        stepsize=stepsize_sweep(g);
        phi=deg2rad(-180:stepsize:180); 
        theta=deg2rad(0:stepsize:180);
        [P_PHI,T_THETA]=meshgrid(phi,theta);

        theta_etilt=deg2rad([45 0 -45]);
        phi_escan=deg2rad([-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]); 
        
        level_db = 60;

        roll_x=360*rand(NO_rotation,1);
        roll_y=180*rand(NO_rotation,1);
        roll_z=360*rand(NO_rotation,1);

        roll_x=deg2rad(roll_x);
        roll_y=deg2rad(roll_y);
        roll_z=deg2rad(roll_z);

        XX = sin(theta')*cos(phi).*(ones(size(P_PHI)));
        YY = sin(theta')*sin(phi).*(ones(size(P_PHI)));
        ZZ = cos(theta')*ones(1,length(phi)).*(ones(size(P_PHI)));

        [azimuth,elevation]=Sampling_R_calculation(NO_rotation,XX,YY,ZZ,roll_x,roll_y,roll_z);

        iPattern_P=zeros(size(XX));
        iPattern_N=zeros(size(XX));
        Combined_pattern_P_N=zeros(size(XX));
        peak_EIRP=zeros(NO_rotation,1);

        for i=1:NO_rotation
            [Ae_H_P,Ae_V_P,Ae_P,iPattern_P]=Combined_Pattern_P(phi_3dB_element,theta_3dB_element,Am,SLAv,Ge_max,N_H,N_V,d_H,d_V,azimuth(:,:,i),elevation(:,:,i),theta_etilt,phi_escan,level_db);
            [Ae_H_N,Ae_V_N,Ae_N,iPattern_N]=Combined_Pattern_N(phi_3dB_element,theta_3dB_element,Am,SLAv,Ge_max,N_H,N_V,d_H,d_V,azimuth(:,:,i),elevation(:,:,i),theta_etilt,phi_escan,level_db);

            Combined_pattern_P_N=max(iPattern_P,(iPattern_N-5));

            XXXXX(:,:,i)=Combined_pattern_P_N;
            peak_EIRP(i)=max(max(Combined_pattern_P_N)); 
        end
        peak_H(15*(k-1)+g,:) = peak_EIRP;
    end
end
save('peak_horization.mat','peak_H');


N_V_sweep = [4,8];
for k = 1:2
    N_V=N_V_sweep(k)
    N_H=2; 
    d_H=0.5; 
    d_V=0.5;
    stepsize_sweep = [2.5,3,3.6,4,4.5,5,6,7.5,9,10,12,15,18,20,30];
    for g = 1:15
        stepsize=stepsize_sweep(g);
        phi=deg2rad(-180:stepsize:180);
        theta=deg2rad(0:stepsize:180); 
        [P_PHI,T_THETA]=meshgrid(phi,theta);

        phi_etilt=deg2rad([45 0 -45]);
        theta_escan=deg2rad([-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]); 
        level_db = 60;

        roll_x=360*rand(NO_rotation,1);
        roll_y=180*rand(NO_rotation,1);
        roll_z=360*rand(NO_rotation,1);

        roll_x=deg2rad(roll_x);
        roll_y=deg2rad(roll_y);
        roll_z=deg2rad(roll_z);

        XX = sin(theta')*cos(phi).*(ones(size(P_PHI)));
        YY = sin(theta')*sin(phi).*(ones(size(P_PHI)));
        ZZ = cos(theta')*ones(1,length(phi)).*(ones(size(P_PHI)));

        [azimuth,elevation]=Sampling_R_calculation(NO_rotation,XX,YY,ZZ,roll_x,roll_y,roll_z);

        iPattern_P=zeros(size(XX));
        iPattern_N=zeros(size(XX));
        Combined_pattern_P_N=zeros(size(XX));
        peak_EIRP=zeros(NO_rotation,1);

        for i=1:NO_rotation
            [Ae_H_P,Ae_V_P,Ae_P,iPattern_P]=Combined_Pattern_P(phi_3dB_element,theta_3dB_element,Am,SLAv,Ge_max,N_H,N_V,d_H,d_V,azimuth(:,:,i),elevation(:,:,i),theta_etilt,phi_escan,level_db);

            [Ae_H_N,Ae_V_N,Ae_N,iPattern_N]=Combined_Pattern_N(phi_3dB_element,theta_3dB_element,Am,SLAv,Ge_max,N_H,N_V,d_H,d_V,azimuth(:,:,i),elevation(:,:,i),theta_etilt,phi_escan,level_db);

            Combined_pattern_P_N=max(iPattern_P,(iPattern_N-5));
            peak_EIRP(i)=max(max(Combined_pattern_P_N));
        end
        peak_V(15*(k-1)+g,:) = peak_EIRP;
    end
end
save('peak_vertical.mat','peak_V');

%% data calculation
load peak_horization.mat;
peak1 = peak_H;
load peak_vertical.mat;
peak2 = peak_V;

for i = 1:15
    p(i,1) = mean(peak1(i,:))-70.5309;
    p(i,2) = mean(peak2(i,:))-70.5309;
    p(i,3) = mean(peak1(15+i,:))-73.5412;
    p(i,4) = mean(peak2(15+i,:))-73.5412;  
    
    p(i,5) = std(peak1(i,:)-70.5309);
    p(i,6) = std(peak2(i,:)-70.5309);
    p(i,7) = std(peak1(15+i,:)-73.5412);
    p(i,8) = std(peak2(15+i,:)-73.5412); 

    error_peak = 70.5309 - peak1(i,:); 
    [YCDF,XCDF,n] = cdfcalc(error_peak);
    p(i,9) = XCDF(0.95*50000); 
    
    error_peak = 70.5309 - peak2(i,:); 
    [YCDF,XCDF,n] = cdfcalc(error_peak);
    p(i,10) = XCDF(0.95*50000);
    
    error_peak = 73.5412 - peak1(i+15,:); 
    [YCDF,XCDF,n] = cdfcalc(error_peak);
    p(i,11) = XCDF(0.95*50000);
    
    error_peak = 73.5412 - peak2(i+15,:); 
    [YCDF,XCDF,n] = cdfcalc(error_peak);
    p(i,12) = XCDF(0.95*50000);  
end
data_result_composite = abs(p);
