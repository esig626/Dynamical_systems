function dx = Volume2(t,x,p,p2)


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

c_no = 7;


w = x(1);
Na = x(2);
K = x(3);
Cl = x(4);
HCO3 = x(5);
H = x(6);
%%%%%%
%%%%%%
%%%%%%
Nal = ones(1,p.ind{c_no});
for i = 1:p.ind{c_no}
    Nal(i)=x(6+i);
end

Kl = ones(1,p.ind{c_no});
for i = 1:p.ind{c_no}
    Kl(i)=x(6+length(Nal)+i);
end

Cll = ones(1,p.ind{c_no});
for i= 1:length(Nal)
    Cll(i) = Nal(i)+Kl(i);
end

Qb = p2.B2 * ( 2 * ( Na + K + H ) + p2.CO20 - p.Ie);
PrCl = zeros(1,length(p.neigh{c_no}));
VCl = PrCl; VtNa = VCl; VtK = VtNa;
Qa = PrCl; Qt = Qa; Qtot=Qa;

for i = 1:length(p.neigh{c_no})
    neighbour = p.neigh{c_no}(i);
    Sa_p=p.Sa_p{c_no,neighbour};
    
    PCl = 1 ./ ( 1 + ( p2.KCaCC ./ ca ).^p2.eta1 ); 
    PrCl(1,i) = Sa_p * PCl;
    VCl(1,i) = Sa_p * p2.RTF * log( Cll(i) / Cl );
    VtNa(1,i) = Sa_p * p2.RTF * log( Nal(i) / p2.Nae );
    VtK(1,i) = Sa_p * p2.RTF * log( Kl(i) / p2.Ke );
    Ii = 2*(Na + K + H) + p2.CO20;
    Il = 2*Cll(i)+p2.Ul;
    Qa(1,i) = Sa_p * p2.B1 * (Il - Ii);
    Qt(1,i) = Sa_p * p2.B3 * ( Il - p.Ie ); 
    Qtot(1,i) =(Qa(1,i)+ Qt(1,i));
end
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
                         
PKb = 1 ./ ( 1 + ( p2.KCaKC ./ ca ).^p2.eta2 ); 
PrKb = PKb; 
VK = p2.RTF * log( p2.Ke / K );
JNaK = p2.Sb*p2.aNaK*(p2.r*p2.Ke^2*Na^3/(p2.Ke^2+p2.alpha1* Na^3));
iNaK = p2.F * JNaK;

GtX = p2.GtNa * p2.St;
GtY = p2.GtK * p2.St;
GY  = p2.GK  * PrKb;
ND  = p2.St*(p2.GtNa+p2.GtK);
vtx = sum(VtNa);
vty = sum(VtK);
Gz  = p2.GCl * sum(PrCl);
Gzv = p2.GCl * dot(sum(PrCl),sum(VCl));

Va=-(ND*iNaK-GY*(GtX*vtx+GtY*vty+ND*VK)...
    +Gzv*(GY+ND))/(ND*(GY+Gz)+GY*Gz);

Vb=-((ND+Gz)*(iNaK-GY*VK)+Gz*(GtX*vtx+GtY*vty)...
    +ND*Gzv)/(ND*GY+Gz*(GY+ND));
Vt = Va - Vb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JCl = Qa; JtNa = JCl; JtK=JtNa;
for i = 1:length(p.neigh{c_no})
    neighbour = p.neigh{c_no}(i);
    Sa_p=p.Sa_p{c_no,neighbour};
    JCl(1,i) = p2.GCl * PrCl(1,i) * (Va + (VCl(1,i)/Sa_p))/p2.F;
    JtNa(1,i) = Sa_p * p2.GtNa * p2.St * ( Vt - (VtNa(1,i)/Sa_p )) / p2.F;
    JtK(1,i) = Sa_p * p2.GtK * p2.St * ( Vt - (VtK(1,i)/Sa_p ) ) / p2.F;  
end

JK = p2.GK * PrKb * ( Vb - VK ) / p2.F;          


JNkcc1=p2.aNkcc*p2.Sb*(p2.a1-p2.a2*Na*K*Cl^2)/(p2.a3+p2.a4*Na*K*Cl^2);     
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JAe4 = p2.Sb * p2.G4 * ( ( p2.Cle / ( p2.Cle + p.KCl ) ) ...
          * ( Na / ( Na + p2.KNa ) ) ...
                    * ( HCO3 / ( HCO3 + p2.KB ) )^2 );       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ / H+ Anion exchanger (Nhe1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JNhe1 = p2.Sb * p2.G1 * ( ( p2.Nae / ( p2.Nae + p2.KNa ) ) ...
            * ( H / ( p2.KH + H ) )...
                          - ( Na / ( Na + p2.KNa ) ) ...
                          * ( p2.He / ( p2.KH + p2.He ) ) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicarbonate Buffer (Reaction Term)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


JBB = w * p2.GB * ( p2.kp * p2.CO20 - p2.kn * HCO3 * H );                   
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jw = Qb - sum(Qa);

dx(1) = Jw;

dx(2) = ( JNkcc1 - 3 * JNaK + JNhe1 - JAe4 - Jw * Na ) / w;

dx(3) = ( JNkcc1 + 2 * JNaK - JK - Jw * K ) / w;

dx(4) = ( 2 * JNkcc1 + JAe4 + sum(JCl) - Jw * Cl ) / w;

dx(5) = ( JBB - 2 * JAe4 - Jw * HCO3 ) / w;

dx(6) = ( JBB - JNhe1 - Jw * H ) / w;

for i = 1:p.ind{c_no}
    dx(6+i)=JtNa(1,i) - Qtot(1,i) * Nal(i);
    dx(6+length(Nal)+i) = JtK(1,i) - Qtot(1,i) * Kl(i);
end

dx = dx';
end