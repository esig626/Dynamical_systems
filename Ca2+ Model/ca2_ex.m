%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
run('Parameters_ca2.m')

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the error allowance of the system:
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [InsP3]cytosol
ip = 0.353;
% Constant Membrane Potential (This will likely change.) 
% I suggest we use the electric circuit cell model
Vm = 170;  
% Turn the pore off/on independent of Ca2+
mPTP = 1;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the system:
[t, x] = ode15s(@(t,y) ca2(t,y,par,ip,Vm,mPTP),0:500,IC,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign Variables:
c   = x(:,1); % [Ca2+]cyt
e  = x(:,2);  % [Ca2+]ER
m  = x(:,3);  % [Ca2+]Mitochondria
u = x(:,4);   % [Ca2+]microdomain
h   = x(:,5); 
h_u = x(:,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective Ca2+ Concentration & mPTP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_eff = (par.volCt .* c + par.volMd .* u)/(par.volCt + par.volMd);
Jx = mPTP * par.kx * (e - c)*exp(par.p2*Vm);

Jxu = mPTP * par.kx * (m - u)*exp(par.p2*Vm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the model's results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all

figure(1)
subplot(2,2,1)
ax = gca;
plot(t,c_eff,'LineWidth',2,'Color','k')
ax.LineWidth = 2;
ax.FontSize = 20;
ylabel('[Ca^{2+}]_{cyt}^{eff} \muM')
xlabel('Time (sec)')
title(['[InsP_3] = ',num2str(ip),' \muM'])
subplot(2,2,2)
ax = gca;
plot(t,e,'LineWidth',2,'Color','r')
ax.LineWidth = 2;
ax.FontSize = 20;
ylabel('[Ca^{2+}]_{er} \muM')
xlabel('Time (sec)')
title(['[InsP_3] = ',num2str(ip),' \muM'])
subplot(2,2,3)
ax = gca;
plot(t,m,'LineWidth',2,'Color','b')
ax.LineWidth = 2;
ax.FontSize = 20;
ylabel('[Ca^{2+}]_m \muM')
xlabel('Time (sec)')
title(['[InsP_3] = ',num2str(ip),' \muM'])
subplot(2,2,4)
ax = gca;
plot(t,Jx,'LineWidth',2,'Color','k')
hold on
plot(t,Jxu,'LineWidth',2,'Color','r')
ax.LineWidth = 2;
ax.FontSize = 20;
ylabel('mPTP \muM/sec')
xlabel('Time (sec)')
legend boxoff
legend('mPTP','mPTP_{\mud}')
title(['[InsP_3] = ',num2str(ip),' \muM'])
