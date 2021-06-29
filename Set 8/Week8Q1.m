% Week8Q1

clc
format bank

T0 = 25 + 273.15;
T1 = T0;

DHrxn = -802.3*1000 % J/mol

T2 = fsolve(@(T2) -DHrxn + integral(@(T)cpCH4g(T),T0,T1) + 2.6*integral(@(T)cpO2g(T),T0,T1) + 9.78*integral(@(T)cpN2g(T), T0, T1) ...
    - 9.78*integral(@(T)cpN2g(T), T0, T2) - 0.6*integral(@(T)cpO2g(T),T0,T2) - integral(@(T)cpCO2g(T),T0,T2) - ...
    2 * integral(@(T)cpH2Og(T),T0,T2), 100);

fprintf("The adiabatic flame temperature is %.2f K per mole of CH4.\n", T2)

function cp=cpCH4g(T)
% Gas heat capacity (J/(molK)) for CH4 with T in K. Source: GVR Appendix 3
a=3.8387e1; b=-7.36639e-2; c=2.90981e-4; d=-2.63849e-7; e=8.00679e-11; cp=a+b*T+c*T.^2+d*T.^3+e*T.^4;
end

function cp=cpH2Og(T)
% Gas heat capacity (J/(molK)) for H2O with T in K. Source: GVR Appendix 3
a=3.40471e1; b=-9.65064e-3; c=3.29983e-5; d=-2.04467e-8; e=4.30228e-12;
cp=a+b*T+c*T.^2+d*T.^3+e*T.^4;
end

function cp=cpO2g(T)
% Gas heat capacity (J/(molK)) for O2 with T in K. Source: GVR Appendix 3
a=2.98832e1; b=-1.13842e-2; c=4.33779e-5; d=-3.70082e-8; e=1.01006e-11;
cp=a+b*T+c*T.^2+d*T.^3+e*T.^4;
end

function cp=cpCO2g(T)
% Gas heat capacity (J/(molK)) for CO2 with T in K. Source: GVR Appendix 3
a=1.90223e1; b=7.96291e-2; c=-7.37067e-5; d=3.74572e-8; e=-8.13304e-12;
cp=a+b*T+c*T.^2+d*T.^3+e*T.^4;
end

function cp=cpN2g(T)
% Gas heat capacity (J/(molK)) for N2 with T in K. Source: GVR Appendix 3
a = 2.94119e1; 
b = -3.00681e-3; 
c = 5.45064e-6; 
d = 5.13186e-9; 
e = -4.25308e-12;
cp=a+b*T+c*T.^2+d*T.^3+e*T.^4;
end