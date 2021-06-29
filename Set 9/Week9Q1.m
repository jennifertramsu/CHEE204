% Problem Set 9

clc

mw1 = 10; % kg/s
T1 = 60 + 273.15; % K
T2 = 30 + 273.15; % K
dH3 = 57.5; % kJ/kg
dH4 = 137; % kJ/kg
dmw = -0.018; % kg/kg dry air

massAir = fsolve(@(ma)(mw1*integral(@(T)cpH2Ol(T), 273.15, T1)/18 + ma*(dH3 - dH4) - (mw1 + dmw*ma)*integral(@(T)cpH2Ol(T), 273.15, T2)/18), 10) % kg

function cp=cpH2Ol(T) % divide by 18 to convert to kJ/kg

% Liquid heat capacity (J/(molK)) for H2O with T in K % Source: GVR

% Appendix 6

a=1.82964e1;

b=4.72118e-1;

c=-1.33878e-3;

d=1.31424e-6;

cp=a+b*T+c*T.^2+d*T.^3;

end