%% Enter the constants of the small molecule, other molecular and system properties:%%
Ka = 9500;  % Fitting parameter at 10 mM KCl conc.
x0 = 0.63*10^(-9); % Fitting parameter at 10 mM KCl conc.
KD=1/Ka;
eta0 =9.15*10^(-4); %(Pa-s; Viscosity at RT i.e. ~24C)
KB = 1.38*10^(-23);
R=8.314;
T= 297;
D0 = 160*10^(-12);
amolecule= KB*T/(6*22/7*eta0*D0); % Ficoll
C0=1;
kh = 1/x0 * KB*T*log(C0/KD);
aMW = 40*1000*10^(-3);
NA = 6.023*10^23;
m = aMW/NA;

%% Define polymer concentrations, molecular weights and other properties: %%
GN0 = 1.58*10^6;
M0 = 44;
C = 0.5;
C_pol = 0.5/100;
nu =0.588;
b=7.2*10^(-10);
pol_den=1.1*1000;
alpha = 0.49;
Vbar=0.84*10^(-4);
Ne=4/5*pol_den*R*T/(M0/1000)/GN0;
Me=Ne*M0;

MW = [20E3 35E3 200E3 1000E3 2000E3 4000E3 8000E3];
N = MW/M0;
Rg = 0.02* MW.^0.58*10^(-9); %in m
Ccrit=MW./(4/3*22/7*Rg.^3*NA)/1000/1000*100; % wt% or g/ml % MW(transition) =  (0.5*4/3*22/7*10^6/100*6.023*10^23*0.02^3*10^(-27))^(-1/0.74) = 2.51E5 => in ln scale MW = 12.43
a1=0.7;
a2=1.5;
points=length(MW);

%% Estimation of friction factor of PEO chains in semi-dilute solution using Vogel-Fulcher equation: %%
N100=100*1000/44;
l=0.58*10^(-9);
Wl4_inf=28.36; %nm4/ns
Wl4=Wl4_inf * exp(-1090/(T-155));
%Wl4=0.1;%nm4/ns
W=Wl4*(10^(-9))^4/(10^(-9))/l^4;
Zeta_N=3*KB*T/(l^2*W);
Zeta_Ne=Zeta_N*N100^3.4/Ne^3.4;
Zeta_Ne=Zeta_N*N100^3.00/Ne^3.4;
segmental_gamma = Zeta_Ne;

%% Estimation of D_x in concentration regimes based on critical concentration for each MW:%%
for i = 1:points
A2(i) = 16/3*22/7*(Rg(i)*0.65)^3*NA/(MW(i)^2)*10^6; %in mol cc /g^2
Ce(i)= Ne^(3*nu-1)/(MW(i) * A2(i)); %gm/cc
C_Ce_ratio(i)=C_pol/Ce(i);
conc_ratio(i)=C/Ccrit(i);
C_Crit_ratio(i)=C_pol/Ccrit(i);
CN_star(i) = N(i)^(-4/5)/ Vbar;
b_statseg(i)=(6/N(i))^0.5*Rg(i);
Dx_DI = KB*T/(6* pi* eta0* amolecule);

%%Dilute region (C_pol<C_crit):%%

if i<= 2

Rt(i)= 0.65*Rg(i); % Assume this to be Rh (hydrodynamic radius)

%Estimating phi_theoretical & gamma_macro
Xicorr(i)=Rg(i)*(conc_ratio(i))^(nu/(1-3*nu)); % Characteristic length scale x(0)
gamma_micro(i)= exp(a2*(2*amolecule/Xicorr(i))^a1)*eta0;
gamma_macro(i)= exp(a2*(Rg(i)/Xicorr(i))^a1)*eta0;
phi_theoretical(i)= (C_pol*NA/MW(i))*4/3*22/7*(Rt(i))^3*10^6;
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
time = 0.4 * 10^(-3); %Defining the time interval for diffussion co-efficient analysis

%Fitting Segmental friction factor
Zeta(i)=2.95*10^(-6); 


tau_c(i) = Zeta(i)*x0/kh; % Characteristic interaction time
Dx(i) = Dx_hardsphere(i) + 1/8*(kh/Zeta(i)/x0*Xicorr(i))^2*time ;
end


%%Transition region (C_pol~C_crit):%%  

if (i == 3)  
%Estimating phi_theoretical & gamma_macro
Xicorr(i)=Rg(i)*(conc_ratio(i))^(nu/(1-3*nu));
gamma_micro(i)= exp(a2*(2*amolecule/Xicorr(i))^a1)*eta0;
gamma_macro(i)= exp(a2*(Rg(i)/Xicorr(i))^a1)*eta0;
phi_theoretical(i)=900*(MW(i))^(-0.77); % Fitting with 0.77 exponent is better
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
Dx(i) = Dx_hardsphere(i);
time = 0.4 * 10^(-3); %Defining the time interval for diffussion co-efficient analysis


%Fitting Segmental friction factor
Zeta(i)=4.0*10^(-6); 

tau_c(i) = Zeta(i)*x0/kh;  % Characteristic interaction time
Dx(i) = Dx_hardsphere(i) + 1/8*(kh/Zeta(i)/x0*Xicorr(i))^2*time ;


end

%%Semi dilute region (C_pol>C_crit):%%

if i >= 4
%Estimating phi_theoretical & gamma_macro
Xicorr(i)=Rg(i)*(conc_ratio(i))^(nu/(1-3*nu));
gamma_micro(i)= exp(a2*(2*amolecule/Xicorr(i))^a1)*eta0;
gamma_macro(i)= exp(a2*(Rg(i)/Xicorr(i))^a1)*eta0;
phi_theoretical(i)= (1/N(i)*(gamma_macro(i)/eta0 - 1))^(3*nu-1); %This is particularly true for semidilute unentangled regime
GNphi(i)=GN0*phi_theoretical(i)^(3*nu/(3*nu-1));
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
eta_ratio(i)=eta_theoretical(i)/gamma_micro(i);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
atube_theoretical(i) = (4/5* pol_den*R*T*phi_theoretical(i)/GNphi(i)/(M0*10^(-3))*b_statseg(i)^2)^0.5; % Calculating tube size of the network compartments

%Fitting Segmental friction factor
Zeta_Ne(i) = segmental_gamma;
Zeta(i) = Zeta_Ne(i);
time = 0.4 * 10^(-3);  %Defining the time interval for diffussion co-efficient analysis

tau_c(i) = Zeta_Ne(i)*x0/kh; % Characteristic interaction time
Dx(i) = Dx_hardsphere(i) + 1/8*(kh/Zeta_Ne(i)/x0*atube_theoretical(i))^2*time ;

end

end

%% Data plotting of selected parameters %%
Dx_exp = 10^(-10) * [1.60 1.68 1.81 1.51 1.86 1.80 1.71 1.68];
MW_exp = [0 MW(1) MW(2) MW(3) MW(4) MW(5) MW(6) MW(7)];

%Collecting MW and Dx in all 3 regions (i.e. dilute, transition and semidilute)%
Dx_theoretical = [Dx_DI Dx(1) Dx(2) Dx(3) Dx(4) Dx(5) Dx(6) Dx(7)];
MW_log = [0 20E3 35E3 200E3 1000E3 2000E3 4000E3 8000E3];

%Plotting graphs of collected parameters%
figure
ylim([2*10^(-10) 4*10^(-10)])
plot(MW_log,Dx_theoretical,'g-.o', MW_exp, Dx_exp,'r *')
xlabel('PEO Molecular Weight (0.5 wt%)')
ylabel('Diffusion co-efficient of Ficoll molecule (m2/s)')
legend('Theoretical','Experimental','Location','southeast')

figure
plot(log(MW_log(2:8)),tau_c,'b  -.O')
xlabel('PEO Molecular Weight (log scale, 0.5 wt%)')
ylabel('Characteristic interaction time (s) of Ficoll in PEO')

