%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 23.05.2022
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, 
% University of Birmingham, 
% B15 2TT, United Kingdom.
% The Tennant Lab
% Kinetic model with free radical generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy=model(x,~,p,O,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metabolites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atp_c = x(1);
adp_c = x(2);
glu_c = x(3);
g6p_c = x(4);
f6p_c = x(5);
nad_c = x(6);
nadh_c= x(7);
gap_c = x(8);
bpg_c = x(9);
pyr_c = x(10);
lac_c = x(11);
gly_c = x(12);
amp_c = x(13);
fa_c  = x(14);
fac_c = x(15);
tg_c  = x(16);
glr_c = x(17);
cr_c  = x(18);
pcr_c = x(19);
badp_c= x(20);
nad_x = x(21);
nadh_x= x(22);
pyr_x = x(23);
fac_x = x(24);
acoa_x= x(25);
cit_x = x(26);
akg_x = x(27);
sca_x = x(28);
atp_x = x(29);
adp_x = x(30);
suc_x = x(31);
mal_x = x(32);
oaa_x = x(33);
o2_m  = x(34);
co2_m = x(35);
asp_x = x(36);
mal_c = x(37);
oaa_c = x(38);
asp_c = x(39);
akg_c = x(40);
glut_c= x(41);
glut_x= x(42);

H_x    = x(43); % Protons H - Hydrogen with no electrons (Hydrogen Ions).
K_x    = x(44); % Potassium K+ - ions
Mg_x   = x(45); % Magnesium Mg2+ - ions
NADH_x = x(46);% Nicotinamide Dehydrogenase 
QH2    = x(47); % Reduced Ubiquinone - oxidation is loss of electrons & reduction is gain of electrons.  
Cred   = x(48); % Cytochrome reduced
ATP_x  = x(49); % ATP  adenosine tri phosphate
ADP_x  = x(50); % ADP adenosine di phosphate. 
mATP_x = x(51); % Mg2+ bound to ATP
mADP_x = x(52); % Mg2+ bound to ADP
Pi_x   = x(53); % Phosphate
ATP_i  = x(54); % ATP in the inner mitochondrial matrix
ADP_i  = x(55); % ADP in the inner mitohcondrial matrix
AMP_i  = x(56); % AMP in the inner mitochondrial matrix
mATP_i = x(57); % Mg2+ bound to ATP
mADP_i = x(58); % Mg2+ bound to ADP
Pi_i   = x(59); % Phosphate inner mitochondrial 
dPsi   = x(60); % mitochondrial membrane potential
O2     = x(61); %not a variable considered in paper
%%
G = p.glu_e;
F = p.fa_e;
L = p.lac_e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ratios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vc = 1.43;
Vm = 5.43;
psp_c = atp_c/adp_c;
psn_c = adp_c/atp_c;
rsn_c = nad_c/nadh_c;
rsp_c = nadh_c/nad_c;
% 
% psp_m = atp_x/adp_x;
psn_m = adp_x/atp_x;
rsn_m = nad_x/nadh_x;
rsp_m = nadh_x/nad_x;

% Define explicitly solved variables. 
NAD_x = param.NADtot-NADH_x;    % Nicotinamide 
Q = param.Qtot-QH2;             % Ubiquinone
Cox = param.Ctot-Cred;          % Cytochrome
% Description of Mg2+ binding to ATP
% f denote the concentrations of ATP and ADP that are not bound to the Mg2+
% in the matrix and inner membrane space.
fATP_x = ATP_x-mATP_x; % mM or M 
fADP_x = ADP_x-mADP_x;
fATP_i = ATP_i-mATP_i;
fADP_i = ADP_i-mADP_i;

%% Cytoplasm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jglut = (p.Tec_glu*G)/(p.Mec_glu+G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEX1: 1 GLU + 1 ATP -> 1 G6P + 1 ADP + 1 H 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jhex =p.V_glu_g6p*(glu_c /(glu_c + p.Kglu))*(psp_c /(p.KPSp_GLU+psp_c));

% With AMP regulation
Jhex =(336027/24250)*(glu_c /(glu_c + p.Kglu))...
*(psp_c /(p.KPSp_GLU+psp_c)) * (1/(1+(amp_c/0.05)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PGI: 1 G6P -> 1 F6P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jpgi = 362/1425 * (g6p_c /( g6p_c + 0.01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PFK: 1 F6P + 1 ATP -> 1 GAP + 1 ADP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jpfk = p.V_g6p_gap *( f6p_c/( f6p_c + p.Kg6p2))...
    *(psp_c /(p.KPSp_2GAP + psp_c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAPD: 1 GAP + 1 NAD -> 1 BPG + 1 NADH + 1(H+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jgapd =(63327/25000)*(gap_c /( gap_c + p.Kgap))...
    * (rsn_c/(rsn_c + p.KRSn_2GAP));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PYK: 2 BPG + 2 ADP -> 2 PYR + 2 ATP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jpyk =p.V_bpg_pyr * (bpg_c/(bpg_c + p.Kbpg))...
    *(psn_c/(psn_c + p.KPSn_BPG));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PYRtm: Mitochondrial Pyruvate  Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jpyrtm = p.Tcm_pyr*pyr_c/(p.Mcm_pyr+pyr_c) * O;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDH: 2 PYR + 2 NADH + 2 (H) -> 2 LAC + 2 NAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jldhn = p.V_pyr_lac*(pyr_c/(pyr_c + p.Kpyr))...
    *(rsp_c/( rsp_c + p.KRSp_PYR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDH: LAC + NAD -> PYR + NADH + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jldhp = p.V_lac_pyr * (lac_c/(lac_c + p.Klac))...
    * (rsn_c/(p.KRSn_LAC + rsn_c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monocarboxylate Transporter 4 (MCT4 - SLC16A3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jmct4e = p.Tec_lac* L/(p.Mec_lac+ L);
Jmct4c = p.Tce_lac*lac_c/(p.Mce_lac+lac_c); 
Jmct4 = Jmct4e + Jmct4c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Glycogen Storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R1: G6P + ATP -> GLY + ADP + 2Pi
Jglyp = p.V_g6p_gly*(g6p_c/(g6p_c + p.Kg6p))...
    *(psp_c /(p.KPSp_G6P + psp_c));

% R2: GLY + Pi -> G6P
Jglyn = 8/25*(gly_c*p.pi_c/(p.pi_c*gly_c+p.Kgly));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAt: Long Chain FA transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jfat = (10500/21979) * (F - p.s_fa * fa_c);

% Jgpat: TG -> 3 FA + GLR 
Jgpat = 3/25 * (tg_c / (tg_c + p.Ktg));

% Jlpat: 3 FA + GLR + 3 ATP -> TG + 3 ADP + 3 Pi 
Jlpat =(9007199254740992/224955625624568575)*...
    ( fa_c * glr_c/(p.K_fa_tg + fa_c * glr_c))...
    * (psp_c/(p.KPSp_FA + psp_c ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACS: FA + CoA + 1 ATP -> FAC + 1 AMP + PPi
% This inhibits PFK... Think how to integrate this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jfacs = (7777/5000) * (fa_c *p.coa_c/(p.Kffa + fa_c*p.coa_c))...
    *(psp_c/(p.KPSp_FFA + psp_c ));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPT1: Fatty Acyl-CoA Mitochondrial Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jcpt1= (1333/100)*fac_c/(p.Mcm_fac+fac_c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP Creatine Kinase (Cytosolic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CKn: PCr + ADP + (H) -> Cr + ATP 
JCKn = p.V_pcr_cr * (pcr_c /(p.KPcr+pcr_c ))...
    *(psn_c/(p.KPSn_PCr + psn_c));

% CKp: Cr + ATP -> PCr + ADP + (H) 
JCKp = p.V_cr_pcr*( cr_c /(p.KCr + cr_c))...
    *( psp_c /(p.KPSp_Cr + psp_c ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADP Buffer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R24: ADP -> bADP
Jbadp = p.V_adp_badp * adp_c / (p.KfADP + adp_c );

% R25: bADP -> ADP
Jbadpn = p.V_badp_adp * (badp_c / (badp_c + p.KbADP) );

%% Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDH: PYR + CoA + NAD -> ACoA + NADH + CO2 + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jpdh = p.V_pyr_acoa * ( pyr_x * p.coa_m /(pyr_x * p.coa_m + p.Kpyrm))...
    * (rsn_m/(rsn_m + p.KRSn_PYRCoA));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAOX: FAC + 7 FAD + 7CoA + 7 NAD -> 7ACoA + 7 FADH2 + 7 NADH + 7(H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jfaox = (14847/5000)*(fac_x*p.coa_m/(fac_x*p.coa_m + p.Kfacm))...
* (rsn_m)/(p.KRSn_FAC+ rsn_m);

% syms x
% Jfaox = x *(fac_x*p.coa_m/(fac_x*p.coa_m + p.Kfacm))...
% * (rsn_m)/(p.KRSn_FAC+ rsn_m) == 0.14;
% solve(Jfaox,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSm: ACoA + OAA + H2O -> CIT + CoA + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jcitm = p.V_acoa_cit*(acoa_x * oaa_x)/(acoa_x * oaa_x + p.KACoA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICDHm: CIT + NAD -> aKG + NADH + CO2 + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jicdh = p.V_cit_akg * (cit_x /(cit_x + p.KCit))...
    *(rsn_m/(p.KRSn_CIT + rsn_m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AKGD: aKG + CoA + NAD -> SCA + NADH + CO2 + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jakgd = p.V_akg_sca * (akg_x / ( akg_x + p.KaKG ))...
    *(rsn_m /( p.KRSn_aKG + rsn_m ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUCOAS: SCA + ADP + Pi -> SUC + ATP + CoA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jsucoas = p.V_sca_suc * (sca_x /( sca_x + p.KSca ))...
    *(psn_m/(p.KPSn_SCA + psn_m));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUCD: 1 SUC + 1 NAD -> MAL + 1 NADH 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jsucd = p.V_suc_mal* (suc_x/(p.KSuc + suc_x))...
    *(rsn_m/(p.KRSn_SUC + rsn_m));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDH: MAL + NAD -> OAA + NADH + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jmdh = (516/25) * ( mal_x /(p.KMal + mal_x))...
    *(rsn_m / (p.KRSn_MAL + rsn_m));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oxygen Transport (Linear diffusion straight into the mitochondria)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jo2 = (98945/35093)*(p.o2_e-p.s_o2*o2_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATPsynathase: 1/2 O2 + 3 ADP + 3 Pi + NADH + H -> H2O + 3 ATP + NAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jatps = (7183407/40000)*(o2_m/(o2_m + p.KO2))...
    *(psn_m/(p.KPSn_O2 + psn_m))...
    *(rsp_m/(p.KRSp_O2+rsp_m));


% T35: Carbon dioxide (Linear diffusion straight into the mitochondria)
% Careful, it is defined as NEGATIVE! -> Removes respiration waste. 
Jpem_co2 = p.l_co2*(p.co2_e-p.s_co2*co2_m);

%% Malate/Aspartate Shuttle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASPTA: 1 OAAm + 1 GLUTm -> 1 AKGm + 1 ASPm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jasptam = (20534783926028289/27553380242255876) ...
    * ((oaa_x * glut_x)/(oaa_x * glut_x + 0.008));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slc25a13 - Citrin:  1 ASPm + 1 GLUTc + 1 H + 1 Ca2+c -> 1 ASPc + 1 GLUTm + 1 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jcitrin = (1480615901066571940948149343879168 ...
/346397766622065054835823125943925)...
* (asp_x/(asp_x + 0.012))*(glut_c/(glut_c + 0.25));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASPTA: 1 OAAc + 1 GLUTc <- 1 AKGc + 1 ASPc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms x
Jasptac = (657525545596092416/167495870506995375) ...
    * ((asp_c * akg_c)/(asp_c * akg_c + 1.3));
% solve(Jasptac,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDH: MAL + NAD <- OAA + NADH + (H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jmdhc = (1334651773/8100) * ( oaa_c /(p.KMal + oaa_c)) ...
    *(rsp_c / (p.KRSn_MAL + rsp_c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AKGLMALtm =  akg[m] + mal_L[c] <=> akg[c] + mal_L[m] (slc25a11)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jakgmal = (870233/51000) * (akg_x/(akg_x + 0.012))*(mal_c/(mal_c + 0.25));

%% ADP/ATP Expenditure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP:  1 ATP + 1 H2O <- 1 ADP + 1 H + 1 Pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP Hydrolysis: H2O + ATP -> ADP + (H) + Pi
Jwatpp = (62547/400) * atp_c / (p.KATP + atp_c);

% syms x
% Jwatpp = x * atp_c / (p.KATP + atp_c) == DM_ATP + Jant - Jfacs + Jamp;
% solve(Jwatpp, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP:  1 AMP + 1 PPi -> 1 ATP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jamp = (62881/12650) * amp_c / (p.KATP + amp_c);

% syms x
% Jamp = x * amp_c / (p.KATP + amp_c) == 0.14;
% solve(Jamp, x)

% Reactions that need ATP: (CAREFUL WITH THE SIGNS!!)
DM_ATP = -Jhex - Jpfk + Jpyk - Jglyp - 3 * Jlpat -JCKp + JCKn;

% Reactions that produce ADP:
DM_ADP = Jhex + Jpfk - Jpyk + Jglyp + 3 * Jlpat + JCKp - JCKn;

% AMP Reaction
DM_AMP = Jfacs;

% Reactions that need ATP: (CAREFUL WITH THE SIGNS!!)
DM_ATPm = Jsucoas + 3 * Jatps;

% Reactions that produce ADP:
DM_ADPm = - Jsucoas - 3 * Jatps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Adenine nucleotide translocator (ANT - SLC25A4/5/6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jant = p.Tmc_atp * (atp_x / (p.Mmc_atp + atp_x));


% syms x
% Jwatpm = x * atp_x / (p.KATP + atp_x) == 18.08;
% solve(Jwatpm, x)

%% NAD/NADH Expenditure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NADH -> NADH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jredox = Jmdhc;
%(803/900) * (nadh_c / (nadh_c +0.01));

% Reactions that need NAD: (CAREFUL WITH THE SIGNS!!)
DM_NAD = - 2 * Jgapd + Jldhn - Jldhp;

% Reactions that produce NADH: (CAREFUL WITH THE SIGNS!!)
DM_NADH= 2 * Jgapd - Jldhn + Jldhp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reactions that need NAD: (CAREFUL WITH THE SIGNS!!)
DM_NADm = - Jpdh - 7 * Jfaox - Jicdh - Jakgd - Jsucd - Jmdh;
% Reactions that produce NADH: (CAREFUL WITH THE SIGNS!!)
DM_NADHm = Jpdh + 7 * Jfaox + Jicdh + Jakgd + Jsucd + Jmdh;

Jredoxm = (902/75) * (nadh_x / (nadh_x + 0.01));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms x
% Jredoxm = x * (nadh_x / (nadh_x +0.01))==DM_NADHm;
% solve(Jredoxm,x)

% syms x
% Jredox = x * (p.nadh_c / ( p.nadh_c +0.01)) == 0.73;
% solve(Jredox,x)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OXPHOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proton Motive Force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG_H = param.F*dPsi+param.RT*log(param.H_i/H_x); % why fixed H+ concentration in the inner membrane?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dehydrogenase Flux % TCA Cycle CRAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_DH = param.x_DH*(param.r*NAD_x-NADH_x)*(1+Pi_x/param.k_Pi1)/(1+Pi_x/param.k_Pi2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Electron Transfer Chain Complex I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG_C1op = param.dG_C1o-param.RT*log(H_x/1e-7)-param.RT*log(Q/QH2);
J_C1 = param.x_C1*(exp(-(dG_C1op+4*dG_H)/param.RT)*NADH_x-NAD_x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Electron Transpoparam.rt Chain Complex III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG_C3op = param.dG_C3o+2*param.RT*log(H_x/1e-7)-param.RT*log(QH2/Q);
J_C3 = param.x_C3*(1+Pi_x/param.k_Pi3)/(1+Pi_x/param.k_Pi4)*...
    (exp(-(dG_C3op+4*dG_H-2*param.F*dPsi)/(2*param.RT))*Cox-Cred);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Electron Transpoparam.rt Chain Complex IV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG_C4op = param.dG_C4o-2*param.RT*log(H_x/1e-7)-param.RT/2*log(O2/1);
J_C4 = param.x_C4*1/(1+param.k_O2/O2)*Cred/param.Ctot*(exp(-(dG_C4op+2*dG_H)/(2*param.RT))*...
    Cred-Cox*exp(param.F*dPsi/param.RT));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ATP Synthase (Complex V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_F1 = param.x_F1*(exp(-(param.dG_F1o-param.n_A*dG_H)/param.RT)*param.K_DD/param.K_DT*mADP_x*Pi_x-mATP_x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adenine nucleotide translocase (ANT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_x = -0.65*dPsi;
Psi_i = 0.35*dPsi;
mincond = 1e-12;
if  (fADP_i > mincond) || (fATP_i > mincond)
    J_ANT = param.x_ANT*(fADP_i/(fADP_i+fATP_i*exp(-param.F*Psi_i/param.RT))...
        -fADP_x/(fADP_x+fATP_x*exp(-param.F*Psi_x/param.RT)))*fADP_i/(param.k_mADP+fADP_i);
else
    J_ANT = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MgATPx/MgADPx Binding Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_MgATPx = param.x_MgA*(fATP_x*Mg_x-param.K_DT*mATP_x);

J_MgADPx = param.x_MgA*(fADP_x*Mg_x-param.K_DD*mADP_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MgATPi/MgADPi Binding Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_MgATPi = param.x_MgA*(fATP_i*param.Mg_i-param.K_DT*mATP_i);

J_MgADPi = param.x_MgA*(fADP_i*param.Mg_i-param.K_DD*mADP_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ATP/ADP Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_ATP = param.gamma*param.p_A*(param.ATP_e-ATP_i);

J_ADP = param.gamma*param.p_A*(param.ADP_e-ADP_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AMP Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_AMP = param.gamma*param.p_A*(param.AMP_e-AMP_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pi Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_Pi2 = param.gamma*param.x_Pi2*(param.Pi_e-Pi_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2H/Pi Transporter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H2PIi = Pi_i*param.H_i/(param.H_i+param.k_dHPi);
H2PIx = Pi_x*H_x/(H_x+param.k_dHPi);
J_Pi1 = param.x_Pi1*(H_x*H2PIi-param.H_i*H2PIx)/(H2PIi+param.k_PiH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ﻿Adenylate Kinase Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_AKi = param.x_AK*(param.K_AK*ADP_i*ADP_i-AMP_i*ATP_i);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proton Leak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GHK Equation.
J_Hle = param.x_Hle*dPsi*...
    (param.H_i*exp(param.F*dPsi/param.RT)-H_x)/(exp(param.F*dPsi/param.RT)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K+ Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GHK Equation.
J_K = param.x_K*dPsi*(param.K_i*exp(param.F*dPsi/param.RT)-K_x)...
                                /(exp(param.F*dPsi/param.RT)-1);
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H+/K+ Transporter (Antiporter) Na+/H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_KH = param.x_KH*(param.K_i*H_x-K_x*param.H_i);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations of the System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy(43,1) = param.x_buff*H_x*(J_DH-5*J_C1-2*J_C3-4*J_C4+(param.n_A-1)*J_F1+2*J_Pi1+J_Hle-J_KH)/param.W_x;
dy(44,1) = (J_KH+J_K)/param.W_x; % Potassium in the matrix % Change for Na+
dy(45,1) = (-J_MgATPx-J_MgADPx)/param.W_x;
dy(46,1) = (J_DH-J_C1)/param.W_x; % We need all of our NADH from the TCA cycle. 
dy(47,1) = (J_C1-J_C3)/param.W_x;
dy(48,1) = (2*J_C3-2*J_C4)/param.W_i;
dy(49,1) = (J_F1-J_ANT)/param.W_x;
dy(50,1) = (-J_F1+J_ANT)/param.W_x; % - dATPx
dy(51,1) = J_MgATPx/param.W_x;
dy(52,1) = J_MgADPx/param.W_x;
dy(53,1) = (-J_F1+J_Pi1)/param.W_x;
dy(54,1) = (J_ATP+J_ANT+J_AKi)/param.W_i;
dy(55,1) = (J_ADP-J_ANT-2*J_AKi)/param.W_i;
dy(56,1) = (J_AMP+J_AKi)/param.W_i;
dy(57,1) = J_MgATPi/param.W_i;
dy(58,1) = J_MgADPi/param.W_i;
dy(59,1) = (-J_Pi1+J_Pi2)/param.W_i;
dy(60,1) = (4*J_C1+2*J_C3+4*J_C4-param.n_A*J_F1-J_ANT-J_Hle-J_K)/param.C_im;
dy(61,1) = 0;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cytosol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(1,1) = Vc * (DM_ATP + Jant - Jfacs + Jamp - Jwatpp);
dy(2,1) = Vc * (DM_ADP - Jant - Jbadp + Jbadpn + Jwatpp);
dy(3,1) = Vc * (Jglut - Jhex);
dy(4,1) = Vc * (Jhex - Jpgi - Jglyp + Jglyn);
dy(5,1) = Vc * (Jpgi - Jpfk);
dy(6,1) = Vc * (Jredox + DM_NAD);
dy(7,1) = Vc * (DM_NADH - Jredox);
dy(8,1) = Vc * (2 * Jpfk - 2 * Jgapd);
dy(9,1) = Vc * (2 * Jgapd - Jpyk);
dy(10,1)= Vc * (Jpyk - Jpyrtm - Jldhn + Jldhp);
dy(11,1)= Vc * (Jmct4 + Jldhn - Jldhp);
dy(12,1)= Vc * 100 * (Jglyp - Jglyn);
dy(13,1)= Vc * (DM_AMP - Jamp);
dy(14,1)= Vc * (Jfat - Jfacs - 3 * Jlpat + Jgpat);
dy(15,1)= Vc * (Jfacs - Jcpt1);
dy(16,1)= Vc * (3 * Jlpat - Jgpat);
dy(17,1)= Vc * (Jgpat - 3 * Jlpat);
dy(18,1)= Vc * (JCKn - JCKp);     
dy(19,1)= Vc * (JCKp - JCKn);
dy(20,1)= Vc * (Jbadp - Jbadpn); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(21,1) = Vm * (DM_NADm + Jredoxm);
dy(22,1) = Vm * (DM_NADHm - Jredoxm);
dy(23,1) = Vm * (Jpyrtm - Jpdh);
dy(24,1) = Vm * (Jcpt1 - Jfaox);
dy(25,1) = Vm * (8*Jfaox + Jpdh - Jcitm);
dy(26,1) = Vm * (Jcitm - Jicdh);
dy(27,1) = Vm * (Jicdh - Jakgd + Jasptam - Jakgmal);
dy(28,1) = Vm * (Jakgd - Jsucoas);
dy(29,1) = Vm * (DM_ATPm - Jant);
dy(30,1) = Vm * (Jant + DM_ADPm);
dy(31,1) = Vm * (Jsucoas - Jsucd);
dy(32,1) = Vm * (Jsucd - Jmdh + Jakgmal);
dy(33,1) = Vm * (Jmdh - Jcitm - Jasptam);
dy(34,1) = Vm * (Jo2 - Jatps);
dy(35,1) = Vm * (Jpem_co2 + Jpdh + Jicdh + Jakgd);
dy(36,1) = Vm * (Jasptam - Jcitrin);
dy(37,1) = Vc * 100 * (Jmdhc - Jakgmal);
dy(38,1) = Vc * 100 * (Jasptac- Jmdhc);
dy(39,1) = Vc * 100 * (Jcitrin - Jasptac);
dy(40,1) = Vc * 100 * (Jakgmal - Jasptac);
dy(41,1) = Vc * 100 * (Jasptac - Jcitrin);
dy(42,1) = Vm * 100 * (Jcitrin - Jasptam);

end
