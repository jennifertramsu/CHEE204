% Week6Q1
% Hvl = integral(cpV - cpL) + Hvl(Tb)

clc
format bank

T1 = 300+273.15; %C
%T2 = 1.22; %C
Tb = 373.15; %K
HvlNormal = 40656.2; %J/mol
X = 0.9;
p1 = 101.3; % kPa

% T2 in K

T2 = fsolve(@(T2)integral(@(T)cpH2Og(T)./T,T1,T2) - 8.314 * log(exp(16.5362-3985.44/(T2 - 38.9974))/p1) ... 
    -(1-X)*(integral(@(T)cpH2Og(T), Tb, T2) - integral(@(T)cpH2Ol(T), Tb, T2) + HvlNormal)/T2, 300)

Hvl = integral(@(T)cpH2Og(T), Tb, T2) - integral(@(T)cpH2Ol(T), Tb, T2) + HvlNormal;
fprintf("The latent heat of vaporization at T2 = 1.22 deg C is %.4e J/mol.\n", Hvl)

% Week6Q2

workShaft = (integral(@(T)cpH2Og(T), T2, T1) + (1 - X)*Hvl)/(1000 * 0.018);
fprintf("The shaft work is %.f kJ/kg.\n", workShaft)

function cp=cpH2Og(T) 
% GVR, J/(mol K), T in K

a = 3.40471e1;
b = -9.65064e-3;
c = 3.29983e-5; 
d = -2.04467e-8;
e = 4.30228e-12;

cp=a+b*T+c*T.^2+d*T.^3 + e*T.^4; 

end

function cp=cpH2Ol(T)

% Liquid heat capacity (J/(molK)) for H2O with T in K % Source: GVR

% Appendix 6

a = 1.82964e1;

b = 4.72118e-1;

c = -1.33878e-3;

d = 1.31424e-6;

cp = a+b*T+c*T.^2+d*T.^3;

end