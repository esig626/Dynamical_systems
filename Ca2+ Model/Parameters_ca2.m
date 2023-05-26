%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the Ca2+ Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.volCt = 0.85;               % [pL] volume cytosol
par.distance = 100;             % [nm] mean intermembrane distance
par.dist = 1e-07;               % [m]
par.r = 5.8e-07;                % [m]
par.SA = 1.690930829868170e-12;	% [m^2]
par.nMtObject = 200;     	% Number of Mitochondria
par.volMd = 0.033818616597363; % [m^3]
par.volER = 0.1;            % [pL]
par.volMt = 0.05;           % [pL] volume ER
par.Vip3r = 1.591683;     	% [pL] volume mitocondria
par.Vserca = 29.12715;     	% [1/s] max flux of IP3R
par.kserca = 0.1933099;     % [uM/s] max flux of SERCA pump
par.ip3 = 0.20;             % [uM] activation constant for SERCA pump
par.a2 = 0.06051925;     	% [uM] IP3 in the cytosol
par.d1 = 0.037678792;     	% [uM^-1*s^-1] IP3R binding rate ca inhibition sites
par.d2 = 1.3270385685;     	% [uM] IP3R dissociation constant for IP3 sites
par.d3 = 1.741082;          % [uM] IP3R dissociation constant for ca inhibition sites
par.d5 = 0.23896425;     	% [uM] IP3R dissociation constant for IP3 sites
par.Vmcu = 7.53115;     	% [uM] IP3R dissociation constant for ca activation sites
par.kmcu = 1.231125;     	% [uM/s] max rate of ca uptake by MCU
par.Vncx = 119.311875;     	% [uM] half-max rate of ca pumping from c to m
par.kncx = 43.2250920;     	% [uM/s] max rate of ca release through NCX
par.kna = 9.4;              % [uM] activation constant for NCX
par.N = 10;                 % [mM] Na activation constant for MCU
par.N_u = 10;               % [mM] Na in cytosol
par.leak_e_u = 0.0433060425;% [mM] Na in microdomain
par.leak_e_c = 0.0106861775;% [1/s] leak constant from ER to Md
par.leak_u_c = 0.0331703975;% [1/s] leak constant from ER to Ct
par.cI = 0.4864205;     	% [1/s] leak constant from Md to Ct
par.cS = 0.6030505;     	% fraction of IP3R facing microdomain
par.cM = 0.8940475;     	% fraction of SERCA facing microdomain
par.cN = 0.5691925;     	% fraction of MCU facing microdomain
par.bt_c = 153.6134;     	% fraction of mNCX facing microdomain
par.K_c = 11.110898;     	% [uM] total buffer concentration in cytosol
par.bt_e = 11082.390;     	% buffer rate constant ratio (Qi 2015)
par.K_e = 966.8265;     	% [uM] total buffer concentration in ER
par.bt_m = 285123.575;     	% buffer rate constant ratio (Qi 2015)
par.K_m = 697.6137;     	% [uM] total buffer concentration in mitocondria
par.bt_u = 190.73308;     	% buffer rate constant ratio (Qi 2015)
par.K_u = 11.984722;    	% [uM] total buffer concentration in micro-domain
par.Vm = 170;               % buffer rate constant ratio (Qi 2015)
par.Vm0 = 91;               % Inner mitochondrial membrane voltage 
par.Vmcu0 = 0.431001; 
par.F = 96485;              % Faraday's Constant 
par.T = 310;                % Temperature
par.R = 8.31;               % Universal Gas Constant
par.RTF = 0.037453903186988; 
par.kx = 0.008;             % Rate constant of bidirectional Ca2+ mPTP
par.p2 = 0.016;             % Voltage dependence coefficient of Pore Activity
par.Kacg = 0.14;            % [uM] Half Maximal activation of Citrin 
par.q2 = 0.1;               % [uM] Half Maximal indirect deactivation of citrin by mitochondrial Ca2+
par.p4 = 0.01;              % [mV] Voltage dependence coefficient of Ctrin activity
par.Vacg = 25;              % Rate constant of NADH production via 
                            % malate-aspartate shuttle (citrin)
                            
                            
                            


IC = [0.108079226929297; 247.83739835315; 0.166826643134091;...
    0.238938153069428; 0.21005822470814; 0.107295996667707];




clearvars -except par IC
