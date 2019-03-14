function [Y T]=rk45func_thalmodwithin(Carr, initval, Tval, ampaval, gabaAval, gabaBval, leakval, T_ret, Cm)
%% THIS FUNCTION IMPLEMENTS THE RK45 ALGORITHM FOR THE KINETIC MODEL FOR
%% THE EXTENDED MODEL BEING WORKED ON BY DAN AND TOM. FEBRUARY 2014.
% Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];
% initval=[Ginit vinit_tcr vinit_trn vinit_in];
% Tval=[Tmax Kp Vp];
% ampaval=[E_ampa alpha_ampa beta_ampa g_ampa_tcr2trn g_ampa_ret2tcr g_ampa_ret2in];
% gabaAval=[alpha_gaba_a beta_gaba_a g_gaba_a_trn2tcr g_gaba_a_trn2trn E_gaba_a_trn2tcr E_gaba_a_trn2trn g_gaba_a_in2tcr g_gaba_a_in2in E_gaba_a_in2tcr E_gaba_a_in2in];
% gabaBval=[alpha1_gaba_b  beta1_gaba_b  alpha2_gaba_b beta2_gaba_b g_gaba_b E_gaba_b Kd_gaba_b n];
% leakval=[E_leak_tcr g_leak_tcr E_leak_trn g_leak_trn E_leak_in
% g_leak_in];

%% set the default values
yvect = zeros(12,length(T_ret));
yvect(1:9,1)=initval(1);%% initialising the r parameters
yvect(10,1)= initval(2); %% initialising V_tcr
yvect(11,1)=initval(3); %% initialising V_trn
yvect(12,1)=initval(4); %% initialising V_in

I_ampa_ret2tcrArr = zeros(1,length(T_ret)-1);
I_ampa_ret2inArr = zeros(1,length(T_ret)-1);
I_ampa_tcr2trnArr = zeros(1,length(T_ret)-1);
I_gaba_a_trn2tcrArr = zeros(1,length(T_ret)-1);
I_gaba_a_trn2trnArr = zeros(1,length(T_ret)-1);
I_gaba_a_in2tcrArr = zeros(1,length(T_ret)-1);
I_gaba_a_in2inArr = zeros(1,length(T_ret)-1);

    I_gaba_bArr = zeros(1,length(T_ret)-1);
    I_leak_tcrArr = zeros(1,length(T_ret)-1);
    I_leak_trnArr = zeros(1,length(T_ret)-1);
    I_leak_inArr = zeros(1,length(T_ret)-1);
    OArr = zeros(1,length(T_ret)-1);
    T_tcrArr = zeros(1,length(T_ret)-1);
    T_trnArr = zeros(1,length(T_ret)-1);
    T_inArr = zeros(1,length(T_ret)-1);

h=0.001;
n=gabaBval(8);
%% RungeKutta45
for i=1:length(T_ret)
r_ampa_ret2tcr=yvect(1,i); %% ampa for the retinal synapse to tcr
r_ampa_ret2in=yvect(2,i); %% ampa for the retinal synapse to in
r_ampa_tcr2trn=yvect(3,i); %% ampa for the tcr synapse to trn
r_gaba_a_trn2tcr=yvect(4,i);%% gabaA synapse 
r_gaba_a_trn2trn=yvect(5,i);%% gabaA synapse 
r_gaba_a_in2tcr=yvect(6,i);%% gabaA synapse 
r_gaba_a_in2in=yvect(7,i);%% gabaA synapse 
%     r_ampar=yvect(1,i); %% ampa for the retinal synapse 
%     r_ampat=yvect(2,i); %% ampa for the tcr synapse 
%     r_gaba_a=yvect(3,i);%% gabaA synapse 
    r_gaba_b=yvect(8,i);%% gabaB synapse 
    G=yvect(9,i); %% The G protein concentration
    V_tcr = yvect(10,i); %% model output: TCR output
    V_trn=yvect(11,i); %% TRN output
    V_in=yvect(12,i); %% IN output
    
% Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];
% initval=[Ginit vinit_tcr vinit_trn vinit_in];
% Tval=[Tmax Kp Vp];
    
    T_tcr=Tval(1)/(1+exp(-(V_tcr-Tval(3))/Tval(2))); %% transmitter concentration of the TCR synapse
    T_trn=Tval(1)/(1+exp(-(V_trn-Tval(3))/Tval(2))); %% transmitter concentration of the TRN synapse
    T_in=Tval(1)/(1+exp(-(V_in-Tval(3))/Tval(2))); %% transmitter concentration of the IN synapse
% ampaval=[E_ampa alpha_ampa beta_ampa g_ampa_tcr2trn g_ampa_ret2tcr g_ampa_ret2in];
% gabaAval=[alpha_gaba_a beta_gaba_a g_gaba_a_trn2tcr g_gaba_a_trn2trn E_gaba_a_trn2tcr E_gaba_a_trn2trn g_gaba_a_in2tcr g_gaba_a_in2in E_gaba_a_in2tcr E_gaba_a_in2in];
% gabaBval=[alpha1_gaba_b beta1_gaba_b alpha2_gaba_b   beta2_gaba_b g_gaba_b E_gaba_b Kd_gaba_b n];
% leakval=[E_leak_tcr g_leak_tcr E_leak_trn g_leak_trn E_leak_in
% g_leak_in];

    I_ampa_ret2tcr =ampaval(5) * r_ampa_ret2tcr * (V_tcr - ampaval(1));
    I_ampa_ret2in =ampaval(6) * r_ampa_ret2in * (V_in - ampaval(1));
    I_ampa_tcr2trn= ampaval(4) * r_ampa_tcr2trn * (V_trn - ampaval(1));
    
    I_gaba_a_trn2tcr=gabaAval(3) * r_gaba_a_trn2tcr * (V_tcr - gabaAval(5));
    I_gaba_a_trn2trn = gabaAval(4) * r_gaba_a_trn2trn * (V_trn - gabaAval(6));
    I_gaba_a_in2tcr=gabaAval(7) * r_gaba_a_in2tcr * (V_tcr - gabaAval(9));
    I_gaba_a_in2in = gabaAval(8) * r_gaba_a_in2in * (V_in - gabaAval(10));

    O = G^n / (G^n + gabaBval(7));
    I_gaba_b= gabaBval(5) * O * (V_tcr - gabaBval(6));
    I_leak_tcr = leakval(2) * (V_tcr - leakval(1));
    I_leak_trn = leakval(4) * (V_trn - leakval(3));
    I_leak_in = leakval(6) * (V_in - leakval(5));

k1_r_ampa_ret2tcr = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2tcr) - ampaval(3) * r_ampa_ret2tcr);
k1_r_ampa_ret2in = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2in) - ampaval(3) * r_ampa_ret2in);
k1_r_ampa_tcr2trn = h * (ampaval(2) * T_tcr * (1-r_ampa_tcr2trn) - ampaval(3) * r_ampa_tcr2trn);
k1_r_gaba_a_trn2tcr = h * (gabaAval(1) * T_trn * (1-r_gaba_a_trn2tcr) - gabaAval(2) * r_gaba_a_trn2tcr);
k1_r_gaba_a_trn2trn = h * (gabaAval(1) * T_trn * (1-r_gaba_a_trn2trn) - gabaAval(2) * r_gaba_a_trn2trn);
k1_r_gaba_a_in2tcr = h * (gabaAval(1) * T_in * (1-r_gaba_a_in2tcr) - gabaAval(2) * r_gaba_a_in2tcr);
k1_r_gaba_a_in2in = h * (gabaAval(1) * T_in * (1-r_gaba_a_in2in) - gabaAval(2) * r_gaba_a_in2in);
k1_r_gaba_b = h * (gabaBval(1) *  T_trn * (1-r_gaba_b) - gabaBval(2)*r_gaba_b);
k1_G = h * (gabaBval(3) * r_gaba_b - gabaBval(4) * G);
k1_V_tcr = h * -(Carr(4) * I_ampa_ret2tcr + Carr(2) * I_gaba_a_trn2tcr + Carr(3) * I_gaba_b + Carr(7) * I_gaba_a_in2tcr + I_leak_tcr)*(1/Cm);
k1_V_trn = h * -(Carr(1)* I_ampa_tcr2trn + Carr(5) * I_gaba_a_trn2trn + I_leak_trn)*(1/Cm);
k1_V_in = h * -(Carr(6)* I_ampa_ret2in + Carr(8) * I_gaba_a_in2in + I_leak_in)*(1/Cm);


% Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];
    r_ampa_ret2tcr_temp=r_ampa_ret2tcr+(1/4)*k1_r_ampa_ret2tcr;
    r_ampa_ret2in_temp=r_ampa_ret2in+(1/4)*k1_r_ampa_ret2in;
    r_ampa_tcr2trn_temp=r_ampa_tcr2trn+(1/4)*k1_r_ampa_tcr2trn;
    r_gaba_a_trn2tcr_temp=r_gaba_a_trn2tcr+(1/4)*k1_r_gaba_a_trn2tcr;
    r_gaba_a_trn2trn_temp=r_gaba_a_trn2trn+(1/4)*k1_r_gaba_a_trn2trn;
    r_gaba_a_in2tcr_temp=r_gaba_a_in2tcr+(1/4)*k1_r_gaba_a_in2tcr;
    r_gaba_a_in2in_temp=r_gaba_a_in2in+(1/4)*k1_r_gaba_a_in2in;
    r_gaba_b_temp=r_gaba_b+(1/4)*k1_r_gaba_b;
    G_temp=G+(1/4)*k1_G;
    
    V_tcr_temp = V_tcr+(1/4)*k1_V_tcr;
    V_trn_temp = V_trn+(1/4)*k1_V_trn;
    V_in_temp = V_in+(1/4)*k1_V_in;
    
    T_tcr_temp = Tval(1)/(1+exp(-(V_tcr_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TCR synapse
    T_trn_temp = Tval(1)/(1+exp(-(V_trn_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TRN synapse
    T_in_temp = Tval(1)/(1+exp(-(V_in_temp - Tval(3))/Tval(2))); %% transmitter concentration of the IN synapse
    

    I_ampa_ret2tcr_temp = ampaval(5) * r_ampa_ret2tcr_temp * (V_tcr_temp - ampaval(1));
    I_ampa_ret2in_temp =ampaval(6) * r_ampa_ret2in_temp * (V_in_temp - ampaval(1));
    I_ampa_tcr2trn_temp= ampaval(4) * r_ampa_tcr2trn_temp * (V_trn_temp - ampaval(1));
    I_gaba_a_trn2tcr_temp =gabaAval(3) * r_gaba_a_trn2tcr_temp * (V_tcr_temp - gabaAval(5));
    I_gaba_a_trn2trn_temp = gabaAval(4) * r_gaba_a_trn2trn_temp * (V_trn_temp - gabaAval(6));
    I_gaba_a_in2tcr_temp=gabaAval(7) * r_gaba_a_in2tcr_temp * (V_tcr_temp - gabaAval(9));
    I_gaba_a_in2in_temp = gabaAval(8) * r_gaba_a_in2in_temp * (V_in_temp - gabaAval(10));
    
    O_temp = G_temp^n / (G_temp^n + gabaBval(7));
    I_gaba_b_temp = gabaBval(5) * O_temp * (V_tcr_temp - gabaBval(6));
    I_leak_tcr_temp = leakval(2) * (V_tcr_temp - leakval(1));
    I_leak_trn_temp = leakval(4) * (V_trn_temp - leakval(3));
    I_leak_in_temp = leakval(6) * (V_in_temp - leakval(5));
    
k2_r_ampa_ret2tcr = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2tcr_temp) - ampaval(3) * r_ampa_ret2tcr_temp);
k2_r_ampa_ret2in =  h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2in_temp) - ampaval(3) * r_ampa_ret2in_temp);
k2_r_ampa_tcr2trn = h * (ampaval(2) * T_tcr_temp * (1-r_ampa_tcr2trn_temp) - ampaval(3) * r_ampa_tcr2trn_temp);
k2_r_gaba_a_trn2tcr = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2tcr_temp) - gabaAval(2) * r_gaba_a_trn2tcr_temp);
k2_r_gaba_a_trn2trn = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2trn_temp) - gabaAval(2) * r_gaba_a_trn2trn_temp);
k2_r_gaba_a_in2tcr = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2tcr_temp) - gabaAval(2) * r_gaba_a_in2tcr_temp);
k2_r_gaba_a_in2in = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2in_temp) - gabaAval(2) * r_gaba_a_in2in_temp);
k2_r_gaba_b = h * (gabaBval(1) *  T_trn_temp * (1-r_gaba_b_temp) - gabaBval(2)*r_gaba_b_temp);
k2_G = h * (gabaBval(3) * r_gaba_b_temp - gabaBval(4) * G_temp);
k2_V_tcr = h * -(Carr(4) * I_ampa_ret2tcr_temp + Carr(2) * I_gaba_a_trn2tcr_temp + Carr(3) * I_gaba_b_temp + Carr(7) * I_gaba_a_in2tcr_temp + I_leak_tcr_temp)*(1/Cm);
k2_V_trn = h * -(Carr(1) * I_ampa_tcr2trn_temp + Carr(5) * I_gaba_a_trn2trn_temp + I_leak_trn_temp)*(1/Cm);
k2_V_in =  h * -(Carr(6) * I_ampa_ret2in_temp + Carr(8) * I_gaba_a_in2in_temp + I_leak_in_temp)*(1/Cm);    
 
% k2_r_ampar = h * (ampaval(2) * T_ret(1,i) * (1-r_ampar_temp) - ampaval(3) * r_ampar_temp);
% k2_r_ampat = h * (ampaval(2) * T_tcr_temp * (1-r_ampat_temp) - ampaval(3) * r_ampat_temp);
% k2_r_gaba_a = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_temp) - gabaAval(2) * r_gaba_a_temp);
% k2_r_gaba_b = h * (gabaBval(1) *  T_trn_temp * (1-r_gaba_b_temp) - gabaBval(2)*r_gaba_b_temp);
% k2_G = h * (gabaBval(3) * r_gaba_b_temp - gabaBval(4) * G_temp);
% k2_V_tcr = h * -(Carr(4) * I_ampar_temp + Carr(2) * I_gaba_a_trn2tcr_temp + Carr(3) * I_gaba_b_temp + I_leak_tcr_temp)*(1/Cm);
% k2_V_trn = h * -(Carr(1)* I_ampat_temp + Carr(5) * I_gaba_a_trn2trn_temp + I_leak_trn_temp)*(1/Cm);

r_ampa_ret2tcr_temp=r_ampa_ret2tcr+(3/32)*k1_r_ampa_ret2tcr+(9/32)*k2_r_ampa_ret2tcr;
r_ampa_ret2in_temp=r_ampa_ret2in+(3/32)*k1_r_ampa_ret2in+(9/32)*k2_r_ampa_ret2in;
r_ampa_tcr2trn_temp=r_ampa_tcr2trn+(3/32)*k1_r_ampa_tcr2trn+(9/32)*k2_r_ampa_tcr2trn;
r_gaba_a_trn2tcr_temp=r_gaba_a_trn2tcr+(3/32)*k1_r_gaba_a_trn2tcr+(9/32)*k2_r_gaba_a_trn2tcr;
r_gaba_a_trn2trn_temp=r_gaba_a_trn2trn+(3/32)*k1_r_gaba_a_trn2trn+(9/32)*k2_r_gaba_a_trn2trn;
r_gaba_a_in2tcr_temp=r_gaba_a_in2tcr+(3/32)*k1_r_gaba_a_in2tcr+(9/32)*k2_r_gaba_a_in2tcr;
r_gaba_a_in2in_temp=r_gaba_a_in2in+(3/32)*k1_r_gaba_a_in2in+(9/32)*k2_r_gaba_a_in2in;
r_gaba_b_temp=r_gaba_b+(3/32)*k1_r_gaba_b+(9/32)*k2_r_gaba_b;
G_temp=G+(3/32)*k1_G+(9/32)*k2_G;
V_tcr_temp = V_tcr+(3/32)*k1_V_tcr+(9/32)*k2_V_tcr;
V_trn_temp = V_trn+(3/32)*k1_V_trn+(9/32)*k2_V_trn;
V_in_temp = V_in+(3/32)*k1_V_in+(9/32)*k2_V_in;


    T_tcr_temp = Tval(1)/(1+exp(-(V_tcr_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TCR synapse
    T_trn_temp = Tval(1)/(1+exp(-(V_trn_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TRN synapse
    T_in_temp = Tval(1)/(1+exp(-(V_in_temp - Tval(3))/Tval(2))); %% transmitter concentration of the IN synapse

    I_ampa_ret2tcr_temp = ampaval(5) * r_ampa_ret2tcr_temp * (V_tcr_temp - ampaval(1));
    I_ampa_ret2in_temp =ampaval(6) * r_ampa_ret2in_temp * (V_in_temp - ampaval(1));
    I_ampa_tcr2trn_temp= ampaval(4) * r_ampa_tcr2trn_temp * (V_trn_temp - ampaval(1));
    I_gaba_a_trn2tcr_temp =gabaAval(3) * r_gaba_a_trn2tcr_temp * (V_tcr_temp - gabaAval(5));
    I_gaba_a_trn2trn_temp = gabaAval(4) * r_gaba_a_trn2trn_temp * (V_trn_temp - gabaAval(6));
    I_gaba_a_in2tcr_temp=gabaAval(7) * r_gaba_a_in2tcr_temp * (V_tcr_temp - gabaAval(9));
    I_gaba_a_in2in_temp = gabaAval(8) * r_gaba_a_in2in_temp * (V_in_temp - gabaAval(10));
    
    O_temp = G_temp^n / (G_temp^n + gabaBval(7));
    I_gaba_b_temp = gabaBval(5) * O_temp * (V_tcr_temp - gabaBval(6));
    I_leak_tcr_temp = leakval(2) * (V_tcr_temp - leakval(1));
    I_leak_trn_temp = leakval(4) * (V_trn_temp - leakval(3));
    I_leak_in_temp = leakval(6) * (V_in_temp - leakval(5));

    
    
k3_r_ampa_ret2tcr = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2tcr_temp) - ampaval(3) * r_ampa_ret2tcr_temp);
k3_r_ampa_ret2in =  h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2in_temp) - ampaval(3) * r_ampa_ret2in_temp);
k3_r_ampa_tcr2trn = h * (ampaval(2) * T_tcr_temp * (1-r_ampa_tcr2trn_temp) - ampaval(3) * r_ampa_tcr2trn_temp);
k3_r_gaba_a_trn2tcr = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2tcr_temp) - gabaAval(2) * r_gaba_a_trn2tcr_temp);
k3_r_gaba_a_trn2trn = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2trn_temp) - gabaAval(2) * r_gaba_a_trn2trn_temp);
k3_r_gaba_a_in2tcr = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2tcr_temp) - gabaAval(2) * r_gaba_a_in2tcr_temp);
k3_r_gaba_a_in2in = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2in_temp) - gabaAval(2) * r_gaba_a_in2in_temp);
k3_r_gaba_b = h * (gabaBval(1) *  T_trn_temp * (1-r_gaba_b_temp) - gabaBval(2)*r_gaba_b_temp);
k3_G = h * (gabaBval(3) * r_gaba_b_temp - gabaBval(4) * G_temp);
k3_V_tcr = h * -(Carr(4) * I_ampa_ret2tcr_temp + Carr(2) * I_gaba_a_trn2tcr_temp + Carr(3) * I_gaba_b_temp + Carr(7) * I_gaba_a_in2tcr_temp + I_leak_tcr_temp)*(1/Cm);
k3_V_trn = h * -(Carr(1) * I_ampa_tcr2trn_temp + Carr(5) * I_gaba_a_trn2trn_temp + I_leak_trn_temp)*(1/Cm);
k3_V_in =  h * -(Carr(6) * I_ampa_ret2in_temp + Carr(8) * I_gaba_a_in2in_temp + I_leak_in_temp)*(1/Cm); 


    r_ampa_ret2tcr_temp=r_ampa_ret2tcr+(1932/2197)*k1_r_ampa_ret2tcr-(7200/2197)*k2_r_ampa_ret2tcr+(7296/2197)*k3_r_ampa_ret2tcr;
    r_ampa_ret2in_temp=r_ampa_ret2in+(1932/2197)*k1_r_ampa_ret2in-(7200/2197)*k2_r_ampa_ret2in+(7296/2197)*k3_r_ampa_ret2in;
    r_ampa_tcr2trn_temp=r_ampa_tcr2trn+(1932/2197)*k1_r_ampa_tcr2trn-(7200/2197)*k2_r_ampa_tcr2trn+(7296/2197)*k3_r_ampa_tcr2trn;
    r_gaba_a_trn2tcr_temp=r_gaba_a_trn2tcr+(1932/2197)*k1_r_gaba_a_trn2tcr-(7200/2197)*k2_r_gaba_a_trn2tcr+(7296/2197)*k3_r_gaba_a_trn2tcr;
    r_gaba_a_trn2trn_temp=r_gaba_a_trn2trn+(1932/2197)*k1_r_gaba_a_trn2trn-(7200/2197)*k2_r_gaba_a_trn2trn+(7296/2197)*k3_r_gaba_a_trn2trn;
    r_gaba_a_in2tcr_temp=r_gaba_a_in2tcr+(1932/2197)*k1_r_gaba_a_in2tcr-(7200/2197)*k2_r_gaba_a_in2tcr+(7296/2197)*k3_r_gaba_a_in2tcr;
    r_gaba_a_in2in_temp=r_gaba_a_in2in+(1932/2197)*k1_r_gaba_a_in2in-(7200/2197)*k2_r_gaba_a_in2in+(7296/2197)*k3_r_gaba_a_in2in;
    r_gaba_b_temp=r_gaba_b+(1932/2197)*k1_r_gaba_b-(7200/2197)*k2_r_gaba_b+(7296/2197)*k3_r_gaba_b;
    G_temp=G+(1932/2197)*k1_G - (7200/2197)*k2_G + (7296/2197)*k3_G;
    V_tcr_temp=V_tcr+(1932/2197)*k1_V_tcr-(7200/2197)*k2_V_tcr+(7296/2197)*k3_V_tcr;
    V_trn_temp=V_trn+(1932/2197)*k1_V_trn-(7200/2197)*k2_V_trn+(7296/2197)*k3_V_trn;
    V_in_temp=V_in+(1932/2197)*k1_V_in-(7200/2197)*k2_V_in+(7296/2197)*k3_V_in;
    
    T_tcr_temp = Tval(1)/(1+exp(-(V_tcr_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TCR synapse
    T_trn_temp = Tval(1)/(1+exp(-(V_trn_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TRN synapse
    Tin_temp = Tval(1)/(1+exp(-(V_in_temp - Tval(3))/Tval(2))); %% transmitter concentration of the IN synapse

    I_ampa_ret2tcr_temp = ampaval(5) * r_ampa_ret2tcr_temp * (V_tcr_temp - ampaval(1));
    I_ampa_ret2in_temp =ampaval(6) * r_ampa_ret2in_temp * (V_in_temp - ampaval(1));
    I_ampa_tcr2trn_temp= ampaval(4) * r_ampa_tcr2trn_temp * (V_trn_temp - ampaval(1));
    I_gaba_a_trn2tcr_temp =gabaAval(3) * r_gaba_a_trn2tcr_temp * (V_tcr_temp - gabaAval(5));
    I_gaba_a_trn2trn_temp = gabaAval(4) * r_gaba_a_trn2trn_temp * (V_trn_temp - gabaAval(6));
    I_gaba_a_in2tcr_temp=gabaAval(7) * r_gaba_a_in2tcr_temp * (V_tcr_temp - gabaAval(9));
    I_gaba_a_in2in_temp = gabaAval(8) * r_gaba_a_in2in_temp * (V_in_temp - gabaAval(10));
    O_temp = G_temp^n / (G_temp^n + gabaBval(7));
    I_gaba_b_temp = gabaBval(5) * O_temp * (V_tcr_temp - gabaBval(6));
    I_leak_tcr_temp = leakval(2) * (V_tcr_temp - leakval(1));
    I_leak_trn_temp = leakval(4) * (V_trn_temp - leakval(3));
    I_leak_in_temp = leakval(6) * (V_in_temp - leakval(5));


k4_r_ampa_ret2tcr = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2tcr_temp) - ampaval(3) * r_ampa_ret2tcr_temp);
k4_r_ampa_ret2in =  h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2in_temp) - ampaval(3) * r_ampa_ret2in_temp);
k4_r_ampa_tcr2trn = h * (ampaval(2) * T_tcr_temp * (1-r_ampa_tcr2trn_temp) - ampaval(3) * r_ampa_tcr2trn_temp);
k4_r_gaba_a_trn2tcr = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2tcr_temp) - gabaAval(2) * r_gaba_a_trn2tcr_temp);
k4_r_gaba_a_trn2trn = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2trn_temp) - gabaAval(2) * r_gaba_a_trn2trn_temp);
k4_r_gaba_a_in2tcr = h * (gabaAval(1) * Tin_temp * (1-r_gaba_a_in2tcr_temp) - gabaAval(2) * r_gaba_a_in2tcr_temp);
k4_r_gaba_a_in2in = h * (gabaAval(1) * Tin_temp * (1-r_gaba_a_in2in_temp) - gabaAval(2) * r_gaba_a_in2in_temp);
k4_r_gaba_b = h * (gabaBval(1) *  T_trn_temp * (1-r_gaba_b_temp) - gabaBval(2)*r_gaba_b_temp);
k4_G = h * (gabaBval(3) * r_gaba_b_temp - gabaBval(4) * G_temp);
k4_V_tcr = h * -(Carr(4) * I_ampa_ret2tcr_temp + Carr(2) * I_gaba_a_trn2tcr_temp + Carr(3) * I_gaba_b_temp + Carr(7) * I_gaba_a_in2tcr_temp + I_leak_tcr_temp)*(1/Cm);
k4_V_trn = h * -(Carr(1) * I_ampa_tcr2trn_temp + Carr(5) * I_gaba_a_trn2trn_temp + I_leak_trn_temp)*(1/Cm);
k4_V_in =  h * -(Carr(6) * I_ampa_ret2in_temp + Carr(8) * I_gaba_a_in2in_temp + I_leak_in_temp)*(1/Cm); 

r_ampa_ret2tcr_temp = r_ampa_ret2tcr+(439/216)*k1_r_ampa_ret2tcr-8*k2_r_ampa_ret2tcr+(3680/513)*k3_r_ampa_ret2tcr-(845/4104)*k4_r_ampa_ret2tcr;
r_ampa_ret2in_temp = r_ampa_ret2in+(439/216)*k1_r_ampa_ret2in-8*k2_r_ampa_ret2in+(3680/513)*k3_r_ampa_ret2in-(845/4104)*k4_r_ampa_ret2in;
r_ampa_tcr2trn_temp = r_ampa_tcr2trn+(439/216)*k1_r_ampa_tcr2trn-8*k2_r_ampa_tcr2trn+(3680/513)*k3_r_ampa_tcr2trn-(845/4104)*k4_r_ampa_tcr2trn;
r_gaba_a_trn2tcr_temp = r_gaba_a_trn2tcr+(439/216)*k1_r_gaba_a_trn2tcr-8*k2_r_gaba_a_trn2tcr+(3680/513)*k3_r_gaba_a_trn2tcr-(845/4104)*k4_r_gaba_a_trn2tcr;
r_gaba_a_trn2trn_temp = r_gaba_a_trn2trn+(439/216)*k1_r_gaba_a_trn2trn-8*k2_r_gaba_a_trn2trn+(3680/513)*k3_r_gaba_a_trn2trn-(845/4104)*k4_r_gaba_a_trn2trn;
r_gaba_a_in2tcr_temp = r_gaba_a_in2tcr+(439/216)*k1_r_gaba_a_in2tcr-8*k2_r_gaba_a_in2tcr+(3680/513)*k3_r_gaba_a_in2tcr-(845/4104)*k4_r_gaba_a_in2tcr;
r_gaba_a_in2in_temp = r_gaba_a_in2in+(439/216)*k1_r_gaba_a_in2in-8*k2_r_gaba_a_in2in+(3680/513)*k3_r_gaba_a_in2in-(845/4104)*k4_r_gaba_a_in2in;
r_gaba_b_temp = r_gaba_b+(439/216)*k1_r_gaba_b-8*k2_r_gaba_b+(3680/513)*k3_r_gaba_b-(845/4104)*k4_r_gaba_b;
G_temp = G+(439/216)*k1_G-8*k2_G+(3680/513)*k3_G-(845/4104)*k4_G;
V_tcr_temp=V_tcr+(439/216)*k1_V_tcr-8*k2_V_tcr+(3680/513)*k3_V_tcr-(845/4104)*k4_V_tcr;
V_trn_temp=V_trn+(439/216)*k1_V_trn-8*k2_V_trn+(3680/513)*k3_V_trn-(845/4104)*k4_V_trn;
V_in_temp=V_in+(439/216)*k1_V_in-8*k2_V_in+(3680/513)*k3_V_in-(845/4104)*k4_V_in;

    T_tcr_temp = Tval(1)/(1+exp(-(V_tcr_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TCR synapse
    T_trn_temp = Tval(1)/(1+exp(-(V_trn_temp - Tval(3))/Tval(2))); %% transmitter concentration of the TRN synapse
    T_in_temp = Tval(1)/(1+exp(-(V_in_temp - Tval(3))/Tval(2))); %% transmitter concentration of the IN synapse

    I_ampa_ret2tcr_temp = ampaval(5) * r_ampa_ret2tcr_temp * (V_tcr_temp - ampaval(1));
    I_ampa_ret2in_temp =ampaval(6) * r_ampa_ret2in_temp * (V_in_temp - ampaval(1));
    I_ampa_tcr2trn_temp= ampaval(4) * r_ampa_tcr2trn_temp * (V_trn_temp - ampaval(1));
    I_gaba_a_trn2tcr_temp =gabaAval(3) * r_gaba_a_trn2tcr_temp * (V_tcr_temp - gabaAval(5));
    I_gaba_a_trn2trn_temp = gabaAval(4) * r_gaba_a_trn2trn_temp * (V_trn_temp - gabaAval(6));
    I_gaba_a_in2tcr_temp=gabaAval(7) * r_gaba_a_in2tcr_temp * (V_tcr_temp - gabaAval(9));
    I_gaba_a_in2in_temp = gabaAval(8) * r_gaba_a_in2in_temp * (V_in_temp - gabaAval(10));
    O_temp = G_temp^n / (G_temp^n + gabaBval(7));
    I_gaba_b_temp = gabaBval(5) * O_temp * (V_tcr_temp - gabaBval(6));
    I_leak_tcr_temp = leakval(2) * (V_tcr_temp - leakval(1));
    I_leak_trn_temp = leakval(4) * (V_trn_temp - leakval(3));
    I_leak_in_temp = leakval(6) * (V_in_temp - leakval(5));


k5_r_ampa_ret2tcr = h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2tcr_temp) - ampaval(3) * r_ampa_ret2tcr_temp);
k5_r_ampa_ret2in =  h * (ampaval(2) * T_ret(i) * (1-r_ampa_ret2in_temp) - ampaval(3) * r_ampa_ret2in_temp);
k5_r_ampa_tcr2trn = h * (ampaval(2) * T_tcr_temp * (1-r_ampa_tcr2trn_temp) - ampaval(3) * r_ampa_tcr2trn_temp);
k5_r_gaba_a_trn2tcr = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2tcr_temp) - gabaAval(2) * r_gaba_a_trn2tcr_temp);
k5_r_gaba_a_trn2trn = h * (gabaAval(1) * T_trn_temp * (1-r_gaba_a_trn2trn_temp) - gabaAval(2) * r_gaba_a_trn2trn_temp);
k5_r_gaba_a_in2tcr = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2tcr_temp) - gabaAval(2) * r_gaba_a_in2tcr_temp);
k5_r_gaba_a_in2in = h * (gabaAval(1) * T_in_temp * (1-r_gaba_a_in2in_temp) - gabaAval(2) * r_gaba_a_in2in_temp);
k5_r_gaba_b = h * (gabaBval(1) *  T_trn_temp * (1-r_gaba_b_temp) - gabaBval(2)*r_gaba_b_temp);
k5_G = h * (gabaBval(3) * r_gaba_b_temp - gabaBval(4) * G_temp);
k5_V_tcr = h * -(Carr(4) * I_ampa_ret2tcr_temp + Carr(2) * I_gaba_a_trn2tcr_temp + Carr(3) * I_gaba_b_temp + Carr(7) * I_gaba_a_in2tcr_temp + I_leak_tcr_temp)*(1/Cm);
k5_V_trn = h * -(Carr(1) * I_ampa_tcr2trn_temp + Carr(5) * I_gaba_a_trn2trn_temp + I_leak_trn_temp)*(1/Cm);
k5_V_in =  h * -(Carr(6) * I_ampa_ret2in_temp + Carr(8) * I_gaba_a_in2in_temp + I_leak_in_temp)*(1/Cm); 

%ORDER 4

      yvect(1,i+1)=yvect(1,i)+(25/216)*k1_r_ampa_ret2tcr+(1408/2565)*k3_r_ampa_ret2tcr+(2197/4101)*k4_r_ampa_ret2tcr-(1/5)*k5_r_ampa_ret2tcr;
      yvect(2,i+1)=yvect(2,i)+(25/216)*k1_r_ampa_ret2in+(1408/2565)*k3_r_ampa_ret2in+(2197/4101)*k4_r_ampa_ret2in-(1/5)*k5_r_ampa_ret2in;
      yvect(3,i+1)=yvect(3,i)+(25/216)*k1_r_ampa_tcr2trn+(1408/2565)*k3_r_ampa_tcr2trn+(2197/4101)*k4_r_ampa_tcr2trn-(1/5)*k5_r_ampa_tcr2trn;
      yvect(4,i+1)=yvect(4,i)+(25/216)*k1_r_gaba_a_trn2tcr+(1408/2565)*k3_r_gaba_a_trn2tcr+(2197/4101)*k4_r_gaba_a_trn2tcr-(1/5)*k5_r_gaba_a_trn2tcr;
      yvect(5,i+1)=yvect(5,i)+(25/216)*k1_r_gaba_a_trn2trn+(1408/2565)*k3_r_gaba_a_trn2trn+(2197/4101)*k4_r_gaba_a_trn2trn-(1/5)*k5_r_gaba_a_trn2trn;
      yvect(6,i+1)=yvect(6,i)+(25/216)*k1_r_gaba_a_in2tcr+(1408/2565)*k3_r_gaba_a_in2tcr+(2197/4101)*k4_r_gaba_a_in2tcr-(1/5)*k5_r_gaba_a_in2tcr;
      yvect(7,i+1)=yvect(7,i)+(25/216)*k1_r_gaba_a_in2in+(1408/2565)*k3_r_gaba_a_in2in+(2197/4101)*k4_r_gaba_a_in2in-(1/5)*k5_r_gaba_a_in2in;
      yvect(8,i+1)=yvect(8,i)+(25/216)*k1_r_gaba_b+(1408/2565)*k3_r_gaba_b+(2197/4101)*k4_r_gaba_b-(1/5)*k5_r_gaba_b;
      yvect(9,i+1)=yvect(9,i)+(25/216)*k1_G+(1408/2565)*k3_G+(2197/4101)*k4_G-(1/5)*k5_G;
      yvect(10,i+1)=yvect(10,i)+(25/216)*k1_V_tcr+(1408/2565)*k3_V_tcr+(2197/4101)*k4_V_tcr-(1/5)*k5_V_tcr;
      yvect(11,i+1)=yvect(11,i)+(25/216)*k1_V_trn+(1408/2565)*k3_V_trn+(2197/4101)*k4_V_trn-(1/5)*k5_V_trn;
      yvect(12,i+1)=yvect(12,i)+(25/216)*k1_V_in+(1408/2565)*k3_V_in+(2197/4101)*k4_V_in-(1/5)*k5_V_in;

      T_trnArr(1,i) = T_trn;
      T_tcrArr(1,i) = T_tcr;
      T_inArr(1,i) = T_in;
    
    I_ampa_ret2tcrArr(1,i) = I_ampa_ret2tcr;
    I_ampa_ret2inArr(1,i) = I_ampa_ret2in;
    I_ampa_tcr2trnArr(1,i) = I_ampa_tcr2trn;
    I_gaba_a_trn2tcrArr(1,i) = I_gaba_a_trn2tcr;
    I_gaba_a_trn2trnArr(1,i) = I_gaba_a_trn2trn;
    I_gaba_a_in2tcrArr(1,i) = I_gaba_a_in2tcr;
    I_gaba_a_in2inArr(1,i) = I_gaba_a_in2in;
    OArr(1,i) = O;
    I_gaba_bArr(1,i) = I_gaba_b;
    I_leak_tcrArr(1,i) = I_leak_tcr;
    I_leak_trnArr(1,i) = I_leak_trn;
    I_leak_inArr(1,i) = I_leak_in;
    


end 
Y=yvect;
T=[T_tcrArr; T_inArr; T_trnArr];
% end   
    
