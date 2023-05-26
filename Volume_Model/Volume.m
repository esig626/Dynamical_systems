function dx = Volume(t,x,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 23.05.2022
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, 
% University of Birmingham, 
% B15 2TT, United Kingdom.
% The Tennant Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if t > 250
    ca = 0.18;
else
    ca = 0.075;
end

Nal = x(1);

Kl = x(2);

Cll = x(3);

w = x(4);

Na = x(5);

K = x(6);

Cl = x(7);

HCO3 = x(8);

H = x(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PrCl = 1 ./ ( 1 + ( par.KCaCC ./ ca ).^par.eta1 ); 
PrKb = 1 ./ ( 1 + ( par.KCaKC ./ ca ).^par.eta2 ); 

JNaK = par.Sb * par.aNaK * ( par.r * par.Ke^2 * Na^3 ... 
                                  / ( par.Ke^2 + par.alpha1 * Na^3 ) );                               
VCl = par.RTF * log( Cll / Cl );        % Nernst Potential. (J/C)1e3 = mV

VK = par.RTF * log( par.Ke / K );     % Nernst Potential. (J/C)1e3 = mV

VtNa = par.RTF * log( Nal / par.Nae );% Nernst Potential. (J/C)1e3 = mV

VtK = par.RTF * log( Kl / par.Ke );   % Nernst Potential. (J/C)1e3 = mV

Va = -(par.F*par.GtK*JNaK*par.St ...
    + par.F*par.GtNa*JNaK*par.St ...
    + par.GCl*par.GK*PrCl*PrKb*VCl ...
    + par.GCl*par.GtK*PrCl*par.St*VCl ...
    + par.GCl*par.GtNa*PrCl*par.St*VCl ...
    - par.GK*par.GtK*PrKb*par.St*VK ...
    - par.GK*par.GtNa*PrKb*par.St*VK ...
    - par.GK*par.GtK*PrKb*par.St*VtK ...
    - par.GK*par.GtNa*PrKb*par.St*VtNa)...
    /(par.GCl*par.GK*PrCl*PrKb ...
    + par.GCl*par.GtK*PrCl*par.St ...
    + par.GCl*par.GtNa*PrCl*par.St ...
    + par.GK*par.GtK*PrKb*par.St ...
    + par.GK*par.GtNa*PrKb*par.St);

Vb = -(par.F*par.GCl*JNaK*PrCl ...
    + par.F*par.GtK*JNaK*par.St ...
    + par.F*par.GtNa*JNaK*par.St ...
    - par.GCl*par.GK*PrCl*PrKb*VK ...
    + par.GCl*par.GtK*PrCl*par.St*VCl ...
    + par.GCl*par.GtNa*PrCl*par.St*VCl ...
    - par.GK*par.GtK*PrKb*par.St*VK ...
    - par.GK*par.GtNa*PrKb*par.St*VK ...
    + par.GCl*par.GtK*PrCl*par.St*VtK ...
    + par.GCl*par.GtNa*PrCl*par.St*VtNa)...
    /(par.GCl*par.GK*PrCl*PrKb ...
    + par.GCl*par.GtK*PrCl*par.St ...
    + par.GCl*par.GtNa*PrCl*par.St ...
    + par.GK*par.GtK*PrKb*par.St ...
    + par.GK*par.GtNa*PrKb*par.St);


Vt = Va - Vb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we calculate the open probabilities for the Ca2+ activated channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Ca2+ - activated Cl- channel

JCl = par.GCl * PrCl * ( Va + VCl ) / par.F;          % fS.micro-metres^2.mV.mol.C^-1

JK = par.GK * PrKb * ( Vb - VK ) / par.F;          % fS.micro-metres^2.mV.mol.C^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tight Junction Na+ and K+ currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JtNa = par.GtNa * par.St * ( Vt - VtNa ) / par.F;                 % fS.micro-metres^2.mV.mol.C^-1

JtK = par.GtK * par.St * ( Vt - VtK ) / par.F;                      % fS.micro-metres^2.mV.mol.C^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Osmolarities 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qa = par.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - par.CO20 + par.Ul );     % micro-metres^3.s^-1

Qb = par.B2 * ( 2 * ( Na + K + H ) + par.CO20 - ...
                      ( par.Nae + par.Ke + par.Cle + par.HCO3e ) );

Qt = par.B3 * ( 2 * ( Nal + Kl ) + par.Ul - ....
                      ( par.Nae + par.Ke + par.Cle + par.HCO3e ) ); % micro-metres^3.s^-1

Qtot=(Qa+Qt);                                     % micro-metres^3.s^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3Na+)/(2K+) ATP-ase pump (NaK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ K+ 2Cl- co-transporter (Nkcc1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JNkcc1 = par.aNkcc1 * par.Sb * ( par.a1 - par.a2 * Na * K * Cl^2 ) ...
                                             / ( par.a3 + par.a4 * Na * K * Cl^2 );     
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JAe4 = par.Sb * par.G4 * ( ( par.Cle / ( par.Cle + par.KCl ) ) * ( Na / ( Na + par.KNa ) ) ...
             * ( HCO3 / ( HCO3 + par.KB ) )^2 );       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ / H+ Anion exchanger (Nhe1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JNhe1 = par.Sb * par.G1 * ( ( par.Nae / ( par.Nae + par.KNa ) ) * ( H / ( par.KH + H ) )...
                          - ( Na / ( Na + par.KNa ) ) * ( par.He / ( par.KH + par.He ) ) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicarbonate Buffer (Reaction Term)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This equation is a reaction inside the cell, note how it depends on the
% cellular volume

JBB = w * par.GB * ( par.kp * par.CO20 - par.kn * HCO3 * H );                   
                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dx(1) = ( JtNa - Qtot * Nal )/par.wl;

dx(2) = ( JtK - Qtot * Kl )/par.wl;

dx(3) = ( - JCl - Qtot * Cll )/par.wl;

dx(4) = Qb - Qa;

dx(5) = ( JNkcc1 - 3 * JNaK + JNhe1 - JAe4 - dx(4) * Na ) / w;

dx(6) = ( JNkcc1 + 2 * JNaK - JK - dx(4) * K ) / w;

dx(7) = ( 2 * JNkcc1 + JAe4 + JCl - dx(4) * Cl ) / w;

dx(8) = ( JBB - 2 * JAe4 - dx(4) * HCO3 ) / w;

dx(9) = ( JBB - JNhe1 - dx(4) * H ) / w;

dx = dx';
end