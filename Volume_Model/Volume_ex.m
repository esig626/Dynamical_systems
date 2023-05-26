clear;close all; clc; format longg; load('Par_V.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 23.05.2022
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, 
% University of Birmingham, 
% B15 2TT, United Kingdom.
% The Tennant Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f_secretion = @(t,x) Volume(t,x,par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Secretion_double
% Integrate

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

[t,U] = ode15s(f_secretion, [0 500], par.IC, options);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nal 	= U(:,1);
Kl 		= U(:,2);
Cll 	= U(:,3);
w 		= U(:,4);
Na 		= U(:,5);
K 		= U(:,6);
Cl 		= U(:,7);
HCO3 	= U(:,8);
H 		= U(:,9);

%%%%%%%

%%%%%%%

figure(1)
set(gca,'Box','off')
plot(t, H, 'LineWidth',2)
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('Volume (\mum^3)')



