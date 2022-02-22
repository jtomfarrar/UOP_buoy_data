clear; close all

%  Compute Callies-Wu GM result with d'Asaro mixed-layer enhancement.

E0=6.3*10^-5;
b=1300;
% b=500;
D=4500;
g=9.81;
% xlat=15;
% xlat=24;
% xlat=30;
xlat=45;
fcor=2.*7.29*10^-5*sind(xlat)
% 3 cph = 5.2 x 10^-3 s^-1 (Munk, 1981, p. 285)
% N0=5.2*10^-3;
%  convert from cph to s^-1
% N0=3*(2*pi/3600)
N0=3.5*(2*pi/3600)
% N0=4*(2*pi/3600)

%  d'Asaro mixed-layer parameters
% Hml=150
% Hml=125
% Hml=100
% Hml=75
% Hml=55
Hml=40
% Hml=20
% drho_rho0=2.00*10^-3
% drho_rho0=1.50*10^-3
% drho_rho0=1.25*10^-3
drho_rho0=1.00*10^-3
% drho_rho0=0.50*10^-3
gamma=g*drho_rho0
Cml=sqrt(gamma*Hml)

zGM=[0];
% zGM=[200:100:2*b];
nGM=length(zGM)

for jGM=1:nGM
    
%  set GM "observation" depth
z_star=zGM(jGM)
%
%  set N at GM depth
Nstar=N0*exp(-z_star/b)
%


% Normalized vertical mode spectrum
jstar=3;
jj=[1:50];
nj=length(jj);
H=(jj.^2+jstar^2).^-1/sum((jj.^2+jstar^2).^-1);

% Normalized frequency spectrum
domg=0.0001*(Nstar-fcor)
omg=[fcor:domg:Nstar]';
nomg=length(omg);
Bomg=(2/pi)*fcor./(omg.*sqrt(omg.^2-fcor^2));
Bint=sum(Bomg(2:nomg))*domg
Bomg=Bomg/Bint;
Bint=sum(Bomg(2:nomg))*domg

%  Kinetic energy spectrum at z=z_star (Callies-Wu eq. 7)
F_z_omg_j=0.5*b^2*N0*Nstar*(omg.^2+fcor^2).*(omg.^-2).*Bomg*H*E0;
% %  Vertical displacement spectrum at z=z_star
% F_z_omg_j=b^2*N0*Nstar^-1*(omg.^2-fcor^2).*(omg.^-2).*Bomg*H*E0;

% Compute horizontal and vertical wavenumber (m^-1)
% WKB phase speed (Callies-Wu text after eq. 7)
cj=N0*b./(jj*pi);
cj(1:5)
% Then get horizontal wave number from modal dispersion relation
%  (Callies-Wu eq. 3)
Kcap=sqrt(omg.^2-fcor^2)./cj;
% Compute domega/dk for conversion of E(omega,j) to E(k,j),
%  from derivative of omega^2=c_j^2*Kcap^2+f^2
domgdk=((cj.^2)./omg).*Kcap;

%  Convert E(omega,j) to E(k,j)
F_z_Kcap_j=F_z_omg_j.*domgdk;

% convert energy spectrum to SSH displacement spectrum (Callies-Wu eq. 6)
F_eta_Kcap_j=2*((omg.^2-fcor^2).^2./(g^2*Kcap.^2.*(omg.^2+fcor^2))).*F_z_Kcap_j;
%  convert from m^-1 and m^2/m^-1 to cpkm and cm^2/cpkm
Kcap_cpkm=(10^3/(2*pi))*Kcap;
F_eta_Kcap_j_cm2cpkm=10^4*(2*pi/10^3)*F_eta_Kcap_j;

% d'Asaro factor (d'Asaro eq. 2)
% d'Asaro factor (using omg.^2 = fcor^2 + (cj.*Kcap).^2 from above)
R_Uml_Uint=1./((Nstar^2-omg.^2).*cj.^-2*Hml^2 + (1-Cml^2./cj.^2).^2);
% R_Uml_Uint=1./(((Nstar^2-omg.^2)./(omg.^2-fcor^2)).*Kcap.^2*Hml^2 ...
%                + (1-Cml^2*Kcap.^2./(omg.^2-fcor^2)).^2);
F_eta_Kcap_j_cm2cpkm_ml=R_Uml_Uint.*F_eta_Kcap_j_cm2cpkm;

%  sum over vertical modes as function of K
%   interpolate to fixed grid in K, K < 10 cpkm (L_K > 100 m)
%  Here and in the following, the 'kk' variables correspond
%   to the radial wavenumber K = (k^2+l^2)^(1/2).
dkk=0.0005;
nkk=round(10/dkk);
kk=dkk*[1:nkk];
F_eta_Kcap=zeros(1,nkk);
F_z_Kcap=zeros(1,nkk);
F_eta_Kcap_ml=zeros(1,nkk);
for jj=1:nj
    F_eta_Kcapj=interp1(Kcap_cpkm(:,jj),F_eta_Kcap_j_cm2cpkm(:,jj),kk);
    F_eta_Kcapj(isnan(F_eta_Kcapj)==1)=0.;
    F_eta_Kcap=F_eta_Kcap+F_eta_Kcapj;
    F_z_Kcapj=interp1(Kcap_cpkm(:,jj),F_z_Kcap_j(:,jj),kk);
    F_z_Kcapj(isnan(F_z_Kcapj)==1)=0.;
    F_z_Kcap=F_z_Kcap+F_z_Kcapj;
    F_eta_Kcapj_ml=interp1(Kcap_cpkm(:,jj),F_eta_Kcap_j_cm2cpkm_ml(:,jj),kk);
    F_eta_Kcapj_ml(isnan(F_eta_Kcapj_ml)==1)=0.;
    F_eta_Kcap_ml=F_eta_Kcap_ml+F_eta_Kcapj_ml;
end

%  Integrate radial wavenumber spectrum F_eta_Kcap(kk) over l
%   Use same wavenumber grid for 1-d wavenumber
nk15=floor((1/15)/dkk);
F_eta_k=zeros(1,nkk);
F_z_k=zeros(1,nkk);
F_eta_k_ml=zeros(1,nkk);
%  cycle over 1-d wavenumber; stop at k=1/15 cpkm
for jkk=1:nk15
    kkj=kk(jkk);
    %  integrate over 2-d wavenumber; stop at l=1/15 cpkm
    for jKcap=jkk+1:nkk
        llj=sqrt(kk(jKcap)^2-kkj^2);
        if(llj < 1/15)
            F_eta_k(jkk)=F_eta_k(jkk)+F_eta_Kcap(jKcap)/llj;
            F_z_k(jkk)=F_z_k(jkk)+F_z_Kcap(jKcap)/llj;
            F_eta_k_ml(jkk)=F_eta_k_ml(jkk)+F_eta_Kcap_ml(jKcap)/llj;
        end
    end
end
F_eta_k=4.*(1/(2.*pi))*F_eta_k*dkk;
F_z_k=4.*(1/(2.*pi))*F_z_k*dkk;
F_eta_k_ml=4.*(1/(2.*pi))*F_eta_k_ml*dkk;

% SWOT noise spectrum vs. cpkm
nSWOT=251;
dkSWOT=(1/15-1/1000)/(nSWOT-1);
kSWOT=[1/1000:dkSWOT:1/15];
E_SWOT=2+1.25*10^-3*kSWOT.^-2;
% Jason noise spectrum vs. cpkm
kJason=[10^-3 10^-2];
E_Jason=[10^2 10^2];

% plot

% remove zeros for log plots

F_eta_Kcap(F_eta_Kcap==0)=NaN;
F_eta_k(F_eta_k==0)=NaN;
F_eta_Kcap_ml(F_eta_Kcap_ml==0)=NaN;
F_eta_k_ml(F_eta_k_ml==0)=NaN;




fsa=16;


figure
%
%  SSH spectra (cm^2/cpkm)
% figure
% subplot(2,2,1)
loglog(kSWOT,E_SWOT,'-r','LineWidth',[2])
hold on
loglog(kJason,E_Jason,'--g','LineWidth',[1])
%  radial spectrum vs. radial wavenumber
% loglog(kk,F_eta_Kcap,'-k')
% % loglog(kk,F_eta_kk,'-k','LineWidth',[2])
loglog(kk,F_eta_k,'-b','LineWidth',[2])
loglog(kk,F_eta_k_ml,'-g','LineWidth',[2])
% axis(10.^[-3 -1 -2 4])
% axis(10.^[-3 -1 0 4])
xlabel('k (cpkm)','FontSize',fsa)
ylabel('SSH PSD (cm^2/cpkm)','FontSize',fsa)
pltnam=strcat('z_*=',num2str(z_star), ...
    'm, N_*/N_0=',num2str(Nstar/N0,2), ...
    '; (H_m_l,\Delta\rho/\rho_0)=(',num2str(Hml),',',num2str(drho_rho0),')');
title(pltnam,'FontSize',fsa)
set(gca,'FontSize',fsa)
%
set(gca,'Color','none')
box on
% pause;

end

cd plots

exportgraphics(gcf,'GM_Nsfc_CWdA.pdf','BackgroundColor','none','ContentType','vector')

cd ..

% 
% cd plots
% 
% for jGM=1:nGM
%     figure(jGM);
%     pause;
%     if(jGM < 10)
%         fname=strcat('GM_N0_CWdA_f0',num2str(jGM))
%     else
%         fname=strcat('GM_N0_CWdA_f',num2str(jGM))
%     end
%     print(fname,'-depsc');
% end
% 
% cd ..


