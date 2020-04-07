%% Enter the constants of the small molecule, other molecular and system properties:%%
amolecule= 0.585 *10^(-9); % Rhodamine 6G dye molecule size from experimental measurements
Ka = 2300;
KD=1/Ka;
eta0 =9.15*10^(-4); %(Pa-s; Viscosity at RT i.e. ~24C)
x0 = 1.*10^(-9);
KB = 1.38*10^(-23);
R=8.314;
T= 297;
aMW = 479*10^(-3);
NA = 6.023*10^23;
m = aMW/NA;
C0=1;
kh = 1/x0 * KB*T*log(C0/KD);

%% Define polymer concentrations, molecular weights and other properties: %%
GN0 = 1.58*10^6;
M0 = 44;
C = 0.5; % polymer wt%
C_pol = 0.5/100;
nu =0.588;
b=7.2*10^(-10);
pol_den=1.1*1000;
alpha = 0.49;
Vbar=0.84*10^(-4);
Ne=4/5*pol_den*R*T/(M0/1000)/GN0;
Me=Ne*M0;
MW = [20E3 35E3 80E3 100E3 150E3 200E3 500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];
N = MW/M0;
Rg = 0.02* MW.^0.58*10^(-9); %in m
Ccrit=MW./(4/3*22/7*Rg.^3*NA)/1000/1000*100; % wt% or g/ml % MW(transition) =  (0.5*4/3*22/7*10^6/100*6.023*10^23*0.02^3*10^(-27))^(-1/0.74) = 2.51E5 => in ln scale MW = 12.43
a1=0.7;
a2=1.5;
points=length(MW);


%% Estimation of friction factor of PEO using Vogel-Fulcher equation: %%
N100=100*1000/44; %For MW = 100K
l=0.58*10^(-9);
Wl4_inf=28.36; %nm4/ns
Wl4=Wl4_inf * exp(-1090/(T-155));
W=Wl4*(10^(-9))^4/(10^(-9))/l^4;
Zeta_N=3*KB*T/(l^2*W);
Zeta_Ne=Zeta_N*N100^3/Ne^3.4;
segmental_gamma = Zeta_Ne;

%% Estimation of D_x in concentration regimes based on critical concentration for each MW:%%
for i = 1:points
    
%%Estimation common to all concentration regimes:%%       
conc_ratio(i)=C/Ccrit(i);
C_Crit_ratio(i)=C_pol/Ccrit(i);
Xicorr(i)=Rg(i)*(conc_ratio(i))^(nu/(1-3*nu));
gamma_micro(i)= exp(a2*(2*amolecule/Xicorr(i))^a1)*eta0;
gamma_macro(i)= exp(a2*(Rg(i)/Xicorr(i))^a1)*eta0;
CN_star(i) = N(i)^(-4/5)/ Vbar;
b_statseg(i)=(6/N(i))^0.5*Rg(i);
Dx_DI = KB*T/(6* pi* eta0* amolecule);


%%Dilute region (C_pol<C_crit):%%
if i<= 4

Rt(i)= 0.65*Rg(i); % Assume this to be Rh (hydrodynamic radius)
Rg_to_Xicorr(i) = Rg(i)/Xicorr(i);
phi_theoretical(i)= (C_pol*NA/MW(i))*4/3*22/7*(Rt(i))^3*10^6;
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
Dx(i) = Dx_hardsphere(i);

end

%%Transition region (C_pol~C_crit):%%
if (i == 5) || (i == 6)

Rg_to_Xicorr(i) = Rg(i)/Xicorr(i);
phi_theoretical(i)=900*(MW(i))^(-0.77); 
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
Dx(i) = Dx_hardsphere(i);

end

%%Semi dilute region (C_pol>C_crit):%%
if i > 6
    
Rg_to_Xicorr(i) = Rg(i)/Xicorr(i);
phi_theoretical(i)= (1/N(i)*(gamma_macro(i)/eta0 - 1))^(3*nu-1); %This is particularly true for semidilute unentangled regime
GNphi(i)=GN0*phi_theoretical(i)^(3*nu/(3*nu-1));
eta_theoretical(i) = eta0 * (1+2.5*phi_theoretical(i)+6.2*phi_theoretical(i)^2);
gamma_theoretical(i) = 6* pi* eta_theoretical(i)* amolecule;
Dx_hardsphere(i)= KB*T/gamma_theoretical(i);
atube_theoretical(i) = (4/5* pol_den*R*T*phi_theoretical(i)/GNphi(i)/(M0*10^(-3))*b_statseg(i)^2)^0.5; 
Zeta_Ne(i) = segmental_gamma;
time = 0.4 * 10^(-3);
Dx(i) = Dx_hardsphere(i) + 1/8*(kh/Zeta_Ne(i)/x0*atube_theoretical(i))^2*time ;


%% Calculation of approximate MSD profile with reptation correction in case of 0.5% PEO MW = 1000K %%
   if i == 9
ah = 1.4 * 10^(-10);
zeta0 = 6*22/7* eta0 * ah;
tau_s(i) = (b)^2 * zeta0  / (3 * (22/7)^2 * KB * T); % segmental relaxation time

deltime=10^(-9); % Time interval of simulation
time=0.02/1000; 
time0=time;

    for j = 1 : 0.4*10^6
    
MSD(j)= 6*(Dx_hardsphere(i) + 1/8*(kh/Zeta_Ne(i)/x0*atube_theoretical(i))^2*time + 1/2 *(b_statseg(i))^2 * Ne^0.5 * (tau_s(i))^(-0.25) * time^(-0.75) )*time ; % MSD with Reptation correction

MSD0=MSD(1);
delMSD_rept(j) = MSD(j)- MSD0;
time = time + deltime;
t_actual(j)=time;
delMSD_rept_time(j) =  t_actual(j) - time0;

    end

  end

end

end

%% Data plotting of selected parameters %%
Dx_exp = 10^(-10) * [4.10 3.55 3.43 3.34 4.10 4.07 3.94];
MW_exp = [0 MW(1) MW(2) MW(6) MW(9) MW(11) MW(15)];

%Collecting MW and Dx in all 3 regions (i.e. dilute, transition and semidilute)%
Dx_theoretical = [Dx_DI Dx(1) Dx(2) Dx(3) Dx(4) Dx(5) Dx(6) Dx(7) Dx(8) Dx(9) Dx(10) Dx(11) Dx(12) Dx(13) Dx(14) Dx(15) Dx(16)];
MW_log = [0 20E3 35E3 80E3 100E3 150E3 200E3 500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];


%Collecting MW, atube and Dx in semidilute entangled region%
MW_overlap = [500E3 800E3 1000E3 1500E3 2000E3 3000E3 4000E3 6000E3 8000E3 10000E3];
atube_theoretical_overlap=[atube_theoretical(7) atube_theoretical(8) atube_theoretical(9) atube_theoretical(10) atube_theoretical(11) atube_theoretical(12) atube_theoretical(13) atube_theoretical(14) atube_theoretical(15) atube_theoretical(16) ]
Dx_theoretical_overlap = [ Dx(7) Dx(8) Dx(9) Dx(10) Dx(11) Dx(12) Dx(13) Dx(14) Dx(15) Dx(16)];

%Plotting graphs of collected parameters%
figure
ylim([2*10^(-10) 4*10^(-10)])
plot(MW_log,Dx_theoretical,'g-.o', MW_exp, Dx_exp,'r *')
xlabel('PEO Molecular Weight (0.5 wt%)')
ylabel('Diffusion co-efficient of Rh6G molecule (m2/s)')
legend('Theoretical','Experimental','Location','southeast')

figure
plot(delMSD_rept_time,delMSD_rept,'b  .')
xlabel('Time interval (s)')
ylabel('MSD of Rh6G in 1000k PEO solution (m2)')
