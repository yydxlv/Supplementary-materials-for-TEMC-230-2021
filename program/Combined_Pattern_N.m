

function [Ae_H,Ae_V,Ae,MagE]=Combined_Pattern_N(phi_3dB_element,theta_3dB_element,Am,SLAv,Ge_max,N_H,N_V,d_H,d_V,P_PHI,T_THETA,theta_etilt,phi_escan,level_db)


Ae_H=min(12*((P_PHI)./phi_3dB_element).^2, Am);  
Ae_H=Ae_H-max(max(Ae_H));  

Ae_V=-min(12*((T_THETA-90)./theta_3dB_element).^2, SLAv); 

Ae=Ge_max-min(-(Ae_H+Ae_V),Am);


TTTT=0;
Aarrar_Beam_linear=0;
for i=1:length(theta_etilt)
    for j=1:length(phi_escan)
SSS=0;
for m=1:N_H
    for n=1:N_V
        v=exp(1i*2*pi*((n-1)*d_V*cos(deg2rad(T_THETA))+(m-1)*d_H*sin(deg2rad(T_THETA)).*sin(deg2rad(P_PHI))));%% super position vector
        W(n,m)=(1/sqrt(N_H*N_V))*exp(1i*2*pi*((n-1)*d_V*sin(theta_etilt(i))-(m-1)*d_H*cos(theta_etilt(i))*sin(phi_escan(j)))); %%weighting
        SSS=SSS+W(n,m)*v;        
    end
end
TTTT=Ae+10*log10(abs(SSS).*abs(SSS));
TTTT_linear=10.^(TTTT./10);
Aarrar_Beam_linear=max(Aarrar_Beam_linear,TTTT_linear);
    end
end

MagE=max((10*log10(Aarrar_Beam_linear) + level_db),0);
