clear; close all

%  Solve Levine 1987 JPO equations for surface-reflected IW
%  Integrate GM spectrum
%
Dcap=4500
%  GM parameters
%
E0=6.3*10^-5
N0_cph=3
N0=N0_cph*2*pi/3600
b=1300
g=9.81;
%
%  require N(z) >= N1 = N(z=-2b)
N1=N0*exp(-2);
%  require N(z) = N300 = N(z=-300 m) for z > -300 m
% N300=N0*exp(-300/b);
% %
% Ncap=N0+0*zz;
% N0=5.24*10^-3
%
% depth above which N decreases or increases by fractional distance to surface
d=200
% d=400
% d=600
% d=1000

%
xlat=30;
% xlat=45;
fcor=2.*7.29*10^-5*sind(xlat)

%

dz=1

% vertical mode number
% jm=1
% jm=2
% jm=3
% jm=5
jm=10;

% omega_cph:  internal wave frequency in cph
% zstar:  GM depth at which to match
%
zstarGM=300;
% nzstar=length(zstarn)
% zstarn=[200:10:2500];
% nzstar=length(zstarn)

%  set GM matching depth
z_star=zstarGM
% z_star=zGM(jGM)
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

% omegan_cph=[0.1:0.02:2.5];
% nomega=length(omegan_cph);

Pmag_fac=zeros(nomg,nj);

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

% Levine surface-reflection factor

zstar=zstarGM;
zz=[-zstar:dz:0];
nz=length(zz);

% % Match at fixed depth
kzstar=nz;

zzb=[-Dcap:dz:0];
nzzb=length(zzb);
Nb=N0*exp(zzb/b);
Nb(Nb<N1)=N1;
% Nb(Nb>N300)=N300;

zz(1);
zzb(nzzb-nz+1);
% Z0_fac=max(abs(Zcap));
kzstarb=nzzb-kzstar+1;

N_fac=10^3;
Z_fac=10^3;
Z_facp=7;
P_facp=1;

%  Buoyancy frequency for Levine surface reflection problem
% Ncap=N0*exp(zz/b);
% Ncap=N0*exp(zz/b).*min(1,abs(zz/d));
DeltaNcap=10*2*pi/3600
Ncap=N0*exp(zz/b)+DeltaNcap.*max(0,(d+zz)/d);


for jm=1:nj
    
    jm

    %  WKB full-depth mode
    % First get vertical mode eigenvalue
%     cj=(N0*b+N1*(Dcap-3*b))./(jm*pi);
    cj_m=cj(jm);
    zz_jm=zeros(size(zzb));
    for kzb=2:nzzb
     zz_jm(kzb)=zz_jm(kzb-1)+Nb(kzb)*dz/cj_m;
    end
    zz_jm=zz_jm-zz_jm(nzzb);
    Pcap=sqrt(Nb./N0).*cos(zz_jm);
    Zcap=-sqrt(Nb./N0).*sin(zz_jm).*Nb/cj_m;
    %
    P_fac=max(abs(Pcap));
    Pcap0=Pcap/P_fac;
    Zcap0=Zcap/P_fac;

    % for kzs=1:nzstar
    % 
    %     zstar=zstarn(kzs);

    Z0_fac=Zcap0(kzstarb);
    P0_fac=Pcap0(kzstarb);


    for momg=1:nomg

        omega=omg(momg);

        %  Levine (1987) eq. (7) - GM WKB dispersion relation
        %  c_j = N0*b/(jm*pi), omega << N0
        kh_jm=(jm*pi/b)*sqrt(omega^2-fcor^2)/N0;

        %  use gamma for Levine beta, kh_jm for Levina alpha
        gamma2=kh_jm^2*((Ncap.^2-omega^2)./(omega^2-fcor^2));
        gamma_star=sqrt(gamma2(1));

        %  matrix operator

        Acap=zeros(nz,3);
        Bcap=zeros(nz,1);
        %
        % psi = 0 at z = 0
        Acap(nz,2)=1;
        % psi_zz + gamma^2 psi = 0
        for kz=2:nz-1
            Acap(kz,1)=1/dz^2;
            Acap(kz,2)=-2/dz^2+gamma2(kz);
            Acap(kz,3)=1/dz^2;
        end
        %  psi_z-i*gamma_star*psi=-2*i*Icap*gamma_star*exp(-i*gamma_star*z) at
        %  z=-zstar
        %   with Icap=1
        Acap(1,3)=1/dz-0.5*i*gamma_star;
        Acap(1,2)=-1/dz-0.5*i*gamma_star;
    %     Acap(1,3)=1/dz;
    %     Acap(1,2)=-1/dz-i*gamma_star;
        Bcap(1)=-2*i*gamma_star*exp(i*gamma_star*zstar);

        %  Solve tridiagonal
        %   Eliminate lower diagonal
        for kz=nz-1:-1:2
        %     Acap(kz,1)=Acap(kz,1)-Acap(kz+1,2)*Acap(kz,1)/Acap(kz+1,2);
            Acap(kz,2)=Acap(kz,2)-Acap(kz+1,3)*Acap(kz,1)/Acap(kz+1,2);
        %     Bcap(kz)=Bcap(kz)    -Bcap(kz+1)  *Acap(kz,1)/Acap(kz+1,2);    
        end
        %
        Acap(1,2)=Acap(1,2)-Acap(1+1,3)*Acap(1,1)/Acap(1+1,2);
        % Bcap(1)=Bcap(1)    -Bcap(1+1)  *Acap(1,1)/Acap(1+1,2);    
        Bcap(1)=Bcap(1)/Acap(1,2);
        %   Back-substitute
        for kz=2:nz
            Bcap(kz)=(Bcap(kz)-Acap(kz,3)*Bcap(kz-1))/Acap(kz,2);
        end

        % P(z) such that d2P/dz2=dZcap/dz
        %  Note that this is evaluated at z=-zstar
        Pcap=zeros(size(Bcap));
        Pcap(1)=-((Bcap(2)-Bcap(1))/dz)/gamma_star^2;
        %
    %     % P(z) such that dP/dz=Zcap
    %     %  Note that this is evaluated at z=-zstar
    %     Pcap=zeros(size(Bcap));
    %     Rcap=Bcap(1)*exp(i*gamma_star*zstar)-exp(2*i*gamma_star*zstar);
    %     Pcap(1)=(exp(i*gamma_star*zstar)-Rcap*exp(-i*gamma_star*zstar))*i/gamma_star;
    %
    %     Pcap(1)=(Bcap(1)-2*exp(-i*gamma_star*zstar))/(i*gamma_star);
    %     Pcap(1)=P0_fac;
        for kz=2:nz
            Pcap(kz)=Pcap(kz-1)+0.5*(Bcap(kz)+Bcap(kz-1))*dz;
        end

        %  normalize
        psi_fac=Z0_fac/Bcap(1);
        Zcap=Bcap*psi_fac;
        Pcap=Pcap*psi_fac;
        %
        Pmag_fac(momg,jm)=abs(Pcap(nz));

    end

    F_eta_Kcap_j_cm2cpkm_ml=Pmag_fac.*F_eta_Kcap_j_cm2cpkm;

end

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


%%


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
    'm, N_*/N_0=',num2str(Nstar/N0,2));
%     '; (H_m_l,\Delta\rho/\rho_0)=(',num2str(Hml),',',num2str(drho_rho0),')');
title(pltnam,'FontSize',fsa)
set(gca,'FontSize',fsa)
%
set(gca,'Color','none')
box on
% pause;

figure
plot(Nb*3600/(2*pi),zzb,'-k')
hold on
plot(Ncap*3600/(2*pi),zz,'-b','LineWidth',[2])
xlabel('N (cph)','FontSize',fsa)
ylabel('z (m)','FontSize',fsa)
set(gca,'FontSize',fsa)

