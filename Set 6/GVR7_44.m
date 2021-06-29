% GVR 7.44

dHvl = 5.83*1000 %cal/gmol
meanCp = 36 % cal/(gmol * K)

Tref = 700 + 273.15; % C
T2 = 400 + 273.15; % C
Tb = 391.4; % K

y = fsolve(@(N2)(N2 * (integral(@(T)AAvap_HC(T), Tref, Tb) - ...
    dHvl + meanCp * (50 + 273.15 - Tb)) - 4*integral(@(T)K_HC(T), Tref, T2) - 43 * integral(@(T)CH4_HC(T), Tref, T2) - ...
    10 * integral(@(T)AAvap_HC(T), Tref, T2) - 43 * integral(@(T)CO2_HC(T), Tref, T2) - ...
    N2 * (integral(@(T)AAvap_HC(T), Tref, T2))), 100)

clear

N2 = 40;

Tb = 391.4; % K
Tref = 700 + 273.15; % C
dHvl = 5.83*1000; %cal/gmol
meanCp = 36; % cal/(gmol * K)

x = fsolve(@(T2)(N2 * (integral(@(T)AAvap_HC(T), Tref, Tb) - ...
    dHvl + meanCp * (50 + 273.15 - Tb)) - 4*integral(@(T)K_HC(T), Tref, T2) - 43 * integral(@(T)CH4_HC(T), Tref, T2) - ...
    10 * integral(@(T)AAvap_HC(T), Tref, T2) - 43 * integral(@(T)CO2_HC(T), Tref, T2) - ...
    N2 * (integral(@(T)AAvap_HC(T), Tref, T2))), 100)

function y=CH4_HC(T) %cal/(gmol * K)

a = 4.750;
b = 1.2e-2;
c = 0.303e-5;
d = -2.63e-9;

y = a + b*T + c*T.^2 + d*T.^3;

end

function y=CO2_HC(T)

a = 6.393;
b = 1.01e-2;
c = -0.3405e-5;

y = a + b*T + c*T.^2;

end

function y=AAvap_HC(T)

a = 8.20;
b = 4.805e-2;
c = -3.056e-5;
d = 8.31e-9;

y = a + b*T + c*T.^2 + d*T.^3;

end

function y=K_HC(T)

a = 4.11;
b = 2.966e-2;
c = -1.793e-5;
d = 4.72e-9;

y = a + b*T + c*T.^2 + d*T.^3;

end