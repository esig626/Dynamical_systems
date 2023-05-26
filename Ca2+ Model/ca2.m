function dx = ca2(~,x,par,ip3,Vm,mPTP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 23.05.2022
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, 
% University of Birmingham, 
% B15 2TT, United Kingdom.
% The Tennant Lab

% DETAILS OF THE MODEL: It is a closed cell model. Thus, there are no 
% Plasma Membrane IN<->OUT Ca2+ Fluxes. Its also well-stirred 
% (aka its homogenous throughout). It contains 4 domains:
% 1. The Cytoplasm  
% 2. A micro-domain subset of the cytoplasm.
% 3. The Endoplasmic Reticulum (or ER)
% 4. A Mitochondrion
% And we keep track of the concentration of Ca2+,in all these compartments
% proportional to the fluxes amongst the 4 compartments. It uses InsP3 as a
% parameter (ip3) to represent the cellular stimulation with agonists.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The variables of the system are: 
% c = [Ca2+] in the cytosol (uM).
% e = [Ca2+] in the ER (uM).
% m = [Ca2+] in the Mitochondrion (uM).
% u = [Ca2+] in a Microdomain (uM)
% h = Slow inactivation gating variable of the ER InsP3 Receptors 
% facing the Cytosol.
% h_u = Slow inactivation gating variable of the ER InsP3 Receptors 
% facing the microdomain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the system variables:
c = x(1);
e = x(2);
m = x(3);
u = x(4);
h = x(5);
h_u = x(6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluxes of the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondrial Permeability Transition Pore Flux (uM/sec)
Jx = mPTP * par.kx * (m - c)*exp(par.p2*Vm);

% Mitochondrial Permeability Transition Pore Flux (uM/sec) -  towards ud
Jxu = mPTP * par.kx * (m - u)*exp(par.p2*Vm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aspartate/Glutamate Exchannger (Gral - Inhibited by Mitochondrial Ca2+)

% Jacg = Vacg *(c/(c + par.Kacg)) * (par.q2/(par.q2 + m)) * exp(par.p4 * Vm);
% Jacgu = par.Vacg *(u/(u + par.Kacg)) * (par.q2/(par.q2 + m)) * exp(par.p4 * Vm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hill function of [InsP3] a& [Ca2+] 
% expressing the probability of the InsP3 Receptors (facing the cytoplasm) 
% being active:
sact = h.*((ip3/(ip3 + par.d1)) .* c./(c + par.d5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total open probability of InsP3 Receptors (facing the cytoplasm):  
Poip3r = sact.^4 + 4*sact.^3.*(1 - sact);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hill function of [InsP3] a& [Ca2+] 
% expressing the probability of the InsP3 Receptors 
% (facing the micro-domain) being active:
sact_u = h.*((ip3/(ip3 + par.d1)) .* u./(u + par.d5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total open probability of InsP3 Receptors (facing the micro-domain): 
Poip3r_u = sact_u.^4 + 4*sact_u.^3.*(1 - sact_u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InsP3 Receptor Fluxes
Jip3r = (1 - par.cI)*(par.Vip3r*Poip3r).*(e - c);
Jip3r_u = par.cI*(par.Vip3r*Poip3r_u).*(e - u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SERCA ATPase Pumps Facing the Cytoplasm
Jserca = (1 - par.cS)*par.Vserca*c.^2./(par.kserca^2 + c.^2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SERCA_u ATPase Pumps Facing the Micro-domain
Jserca_u = par.cS*par.Vserca*u.^2./(par.kserca^2 + u.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondrial Na/Ca2+ Exchanger (NCX) 
% -> Electrogenic Channel (3:1 stoichiometry) facing the cytoplasm
Jncx = (1 - par.cN)*par.Vncx*(par.N^3/(par.kna^3 + par.N^3)).*(m./(par.kncx + m));    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondrial Na/Ca2+ Exchanger (NCX) facing the micro-domain
% -> Electrogenic Channel
Jncx_u = par.cN*par.Vncx*(par.N_u^3/(par.kna^3 + par.N_u^3)).*(m./(par.kncx + m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum rate of ca uptake by Mitochondrial Ca2+ Uniporter
VMCU = par.RTF*par.Vmcu0*(Vm - par.Vm0) * exp(par.RTF*par.Vmcu0*(Vm - par.Vm0))* sinh(par.RTF*par.Vmcu0*(Vm - par.Vm0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondrial Ca2+ Uniporter (facing the cytosol)
Jmcu = (1 - par.cM)*VMCU*(c.^2./(par.kmcu^2 + c.^2));  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mitochondrial Ca2+ Uniporter (facinng the micrdomain)
Jmcu_u = par.cM*VMCU*(u.^2./(par.kmcu^2 + u.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear leaks coupling all domains
Jleak_u_c = par.leak_u_c.*(u - c);
Jleak_e_u = par.leak_e_u.*(e - u);
Jleak_e_c = par.leak_e_c.*(e - c);        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h Slow inactivation gating variable of the ER InsP3 Receptors 
ah = par.a2*par.d2*(ip3 + par.d1)/(ip3 + par.d3);
bh = par.a2*c;        
% h_u Slow inactivation gating variable of the ER InsP3 Receptors 
bh_u = par.a2*u; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heavy Ca2+ buffering
theta_c = par.bt_c*par.K_c./((par.K_c + c).^2); % buffer factor
theta_e = par.bt_e*par.K_e./((par.K_e + e).^2); % buffer factor
theta_m = par.bt_m*par.K_m./((par.K_m + m).^2); % buffer factor
theta_u = par.bt_u*par.K_u./((par.K_u + u).^2); % buffer factor   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations of the system. 

dx(1) = (Jip3r + Jleak_u_c + Jleak_e_c + Jncx...
    - Jserca - Jmcu + Jx )./(1 + theta_c);
dx(2) = (par.volCt/par.volER*(Jserca + Jserca_u ...
    - Jip3r - Jip3r_u - Jleak_e_u - Jleak_e_c ))./(1 + theta_e);
dx(3) = (par.volCt/par.volMt*(Jmcu + Jmcu_u - Jncx...
    - Jncx_u - Jxu - Jx))./(1 + theta_m);
dx(4) = (par.volCt/par.volMd*(Jip3r_u + Jncx_u + Jleak_e_u ...
    - Jserca_u - Jmcu_u - Jleak_u_c + Jxu))./(1 + theta_u);
dx(5) = ah.*(1 - h) - bh.*h;
dx(6) = ah.*(1 - h_u) - bh_u.*h_u;
dx = dx';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end