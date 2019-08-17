%% Enter the constants of the small molecule, other molecular and system properties:%%

T= 296;
R=8.314;
NA = 6.023*10^23;
KB = 1.38*10^(-23);
Ka = 2300;
KD=1/Ka;
C0=1;
eta0 = 9.2*10^(-4); %in Pa-s; solvent viscosity at RT i.e. ~23C
x0 = 1.0*10^(-9);
amolecule= 0.6*10^(-9); % Rhodamine6G dye molecule size
aMW = 479*10^(-3);
m = aMW/NA;
Dx_DI = KB*T/(6* pi* eta0* amolecule);
kh = 1/x0 * KB*T*log(C0/KD);

%% Define polymer concentrations, molecular weights and other properties: %%

M0 = 44;
nu =0.588;
GN0 = 1.58*10^6;
b=7.2*10^(-10);
pol_den=1.1*1000;
Ne=4/5*pol_den*R*T/(M0/1000)/GN0;
KN=Ne^0.5*b^3/(22/7/6*b^3); %Kavassalis-Noolandi overlap parameter for linear PEO chain
MW = [20E3 35E3 80E3 100E3 150E3 200E3 500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];
C = 0.5;
C_pol = 0.5/100;
N = MW/M0;
Rg = 0.02* MW.^0.58*10^(-9); %in m
Ccrit=MW./(4/3*22/7*Rg.^3*NA)/1000/1000*100; % wt% 
a1=0.7;
a2=1.5;
points=length(MW); %size of the MW array

%% Estimation of friction factor of PEO using Vogel-Fulcher equation: %%

N100=100*1000/44;%For MW = 100K
l=0.58*10^(-9);
Wl4_inf=28.36; %nm4/ns
Wl4=Wl4_inf * exp(-1090/(T-155));
W=Wl4*(10^(-9))^4/(10^(-9))/l^4;
Zeta_Ne=3*KB*T/(l^2*W);
Zeta_N=Zeta_Ne*(N100/Ne)^3.4;

%% Estimation of D_x in concentration regimes based on critical concentration for each MW:%%

for i = 1:points
    
%%Estimation common to all regions:%%    
A2(i) = 16/3*22/7*(Rg(i)*0.65)^3*NA/(MW(i)^2)*10^6; %in mol cc /g^2
conc_ratio(i)=C/Ccrit(i);
C_Crit_ratio(i)=C_pol/Ccrit(i);
Xicorr(i)=Rg(i)*(conc_ratio(i))^(nu/(1-3*nu));
gamma_micro(i)= exp(a2*(2*amolecule/Xicorr(i))^a1)*eta0;
gamma_macro(i)= exp(a2*(Rg(i)/Xicorr(i))^a1)*eta0;
b_statseg(i)=(6/N(i))^0.5*Rg(i);

%%Dilute region (C_pol<C_crit):%%
if i<= 4
Rt(i)= 0.65*Rg(i); % Assume this to be Rh (hydrodynamic radius)
phi_theoretical(i)= (C_pol*NA/MW(i))*4/3*22/7*(Rt(i))^3*10^6;
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
Dx(i) = Dx_hardsphere(i);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i); % tracking the ratio of eta_theoretical to microviscosity
KN_ratio(i)=(N(i)/Ne)^0.5;% tracking the ratio of overlap parameter to Kavassalis-Noolandi overlap parameter
end

%%Transition region (C_pol~C_crit):%%
if i == 5 || i==6
phi_theoretical(i)=900*(MW(i))^(-0.76);
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
Dx(i) = Dx_hardsphere(i);
KN_ratio(i)=(N(i)/Ne)^0.5;
end

%%Semi dilute region (C_pol>C_crit):%%
if i > 6
phi_theoretical(i)= (1/N(i)*(gamma_macro(i)/eta0 - 1))^(3*nu-1); %volume fraction in early semidilute regime
GNphi(i)=GN0*phi_theoretical(i)^(3*nu/(3*nu-1));
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i);
KN_ratio(i)=(N(i)/Ne)^0.5;
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
atube_theoretical(i) = (4/5* pol_den*R*T*phi_theoretical(i)/GNphi(i)/(M0*10^(-3))*b_statseg(i)^2)^0.5; 
Dx(i) = Dx_hardsphere(i) + 1/(4*Zeta_N*x0)*kh*atube_theoretical(i)^2 - 1/(16*Zeta_N*x0^2)*kh*atube_theoretical(i)^3;
Dx_difference(i)=1/(4*Zeta_N*x0)*kh*atube_theoretical(i)^2 - 1/(16*Zeta_N*x0^2)*kh*atube_theoretical(i)^3; %tracking the difference between 2nd & 3rd terms of diffusion equation

%% Comparison between v0_theory (SE diffusion)to v0_calc (geometric) to probe the atube lengthscale  %%
v0_theory(i)=2*Dx_hardsphere(i)/atube_theoretical(i);
tau_t(i)=(Xicorr(i))^2/Dx_hardsphere(i)/2; %Average intermesh diffusion time
v0_calc(i)=(Xicorr(i)*0.22)/tau_t(i);
v0_ratio(i)=v0_calc(i)/v0_theory(i); % tracking ratio of v0 from geometric consideration to v0 from SE diffusion

%% Comparison between diffusion times & estimation of diffusion anisotropy %%
tauB(i)=(atube_theoretical(i)/2)^2/Dx(i); %Binding time perpendicular to polymer chain i.e. along x
tauD(i)=3*N(i)^3.4/(W*(22/7)^2)*(l/atube_theoretical(i))^2; % Disengagement/ escape time of the polymer chain from the tube
tau_ratio(i)=tauB(i)/tauD(i); % tracking the ratio of small molecule to PEO chain diffusion time
Dx_anisotropy(i)=Dx_hardsphere(i)/Dx(i); % tracking ratio of perpendicular to parallel diffusion w.r.t. the polymer chain along primitive path

end

end

%% Data plotting of selected parameters %%
Dx_exp = 10^(-10) * [4.10867 3.50332 3.30835 3.14145 3.54093 3.60117 3.74643];
MW_exp = [0 MW(1) MW(2) MW(6) MW(9) MW(11) MW(15)];

%Collecting MW and Dx in all 3 regions (i.e. dilute, transition and semidilute)%
Dx_theoretical_all = [Dx_DI Dx(1) Dx(2) Dx(3) Dx(4) Dx(5) Dx(6) Dx(7) Dx(8) Dx(9) Dx(10) Dx(11) Dx(12) Dx(13) Dx(14) Dx(15) Dx(16)];
MW_all = [0 20E3 35E3 80E3 100E3 150E3 200E3 500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];

%Collecting MW, atube, Dx and Dx_difference in semidilute region%
MW_overlap = [500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];
atube_theoretical_overlap=[atube_theoretical(7) atube_theoretical(8) atube_theoretical(9) atube_theoretical(10) atube_theoretical(11) atube_theoretical(12) atube_theoretical(13) atube_theoretical(14) atube_theoretical(15) atube_theoretical(16)];
Dx_theoretical_overlap = [ Dx(7) Dx(8) Dx(9) Dx(10) Dx(11) Dx(12) Dx(13) Dx(14) Dx(15) Dx(16)];
Dx_difference_overlap = [ Dx_difference(7) Dx_difference(8) Dx_difference(9) Dx_difference(10) Dx_difference(11) Dx_difference(12) Dx_difference(13) Dx_difference(14) Dx_difference(15) Dx_difference(16)];


%Plotting graphs of collected parameters%
figure
plot(MW_overlap,Dx_difference_overlap,'r--O')
xlabel('PEO MW in semidilute region')
ylabel('Difference between 2nd and 3rd term of Dx')


figure
plot(MW_overlap,atube_theoretical_overlap,'b--^')
xlabel('PEO MW in semidilute region')
ylabel('atube (m)')

figure 
plot(MW_overlap, eta_theoretical(7:16), 'r--O')
xlabel('PEO MW in semidilute region')
ylabel('Viscosity (pa-s) in semidilute region')

figure
plot(atube_theoretical_overlap, Dx_theoretical_overlap,'m ^')
xlabel('atube (m)')
ylabel('Dx (m2/s)')

figure
plot(Xicorr(7:16), Dx_theoretical_overlap,'r >')
xlabel('PEO mesh size (m)')
ylabel('Dx (m2/s)')

figure
plot(MW_all,Dx_theoretical_all,'g--o')
xlabel('PEO MW')
ylabel('Dx (m2/s)')
ylim([2.5*10^(-10) 4*10^(-10)])
legend('Theoretical')

figure
plot(MW_all,Dx_theoretical_all,'g-o', MW_exp, Dx_exp,'r *')
xlabel('PEO MW')
ylabel('Dx (m2/s)')
legend('Theoretical','Experimental', 'Location', 'Southeast')
ylim([2.5*10^(-10) 4*10^(-10)])


% Collecting data of v0_ratio, tau_ratio, Dx_anisotropy, eta_ratio and ratio of Kavassalis-Noolandi order parameters in semidilute region %
transpose(v0_ratio(7:16))
transpose(tau_ratio(7:16))
transpose(Dx_anisotropy(7:16))
transpose(eta_ratio(7:16))
transpose(KN_ratio(7:16))
