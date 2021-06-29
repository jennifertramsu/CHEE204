% Week7Q1

clc
format bank

global a;
global b;
global R;

a=5.536*0.1; % Jm^3/mol^2 
b=0.03049/1000; % L/mol->m^3/mol c
R=8.314; % J/(mol K) = m^3*Pa/(mol K)

T1 = 100 + 273.15; % K
T2 = 200 + 273.15; % K
p1 = 1.0131; % bar
p2 = 15 + (16 - 15)/(201.4 - 198.3); % bar
vliq = 0.001044; % m^3/kg, specific volume at reference, 100degC

Hvl_ref = 2256.9; % kJ/kg

% Pressure correction

dHl = integral(@(T)cpH2Ol(T), T1, T2)/18 + 100*vliq*(p2 - p1); % * 100 to convert from m^3/kg * bar to kJ/kg (= m^2/s^2)

dHv = integral(@(p) dHdpT(p, T1), p1*100000, 0)/1000 + integral(@(T) cpH2Og(T), T1 - 273.15, T2 - 273.15)*1000/18 + integral(@(p) dHdpT(p, T2), 0, p2)/1000;

Hvl = dHv - dHl + Hvl_ref; % kJ/kg

fprintf("(Pressure Correction) Heat of vaporization at 200degC is %.2f kJ/kg.\n", Hvl)

% Without pressure correction

dHl = integral(@(T)cpH2Ol(T), T1, T2)/18;

dHv = integral(@(T) cpH2Og(T), T1 - 273.15, T2 - 273.15)*1000/18;

Hvl = dHv - dHl + Hvl_ref; % kJ/kg

fprintf("(No pressure correction) Heat of vaporization at 200degC is %.2f kJ/kg.\n", Hvl)

% Watson's correlation

Tc = 647.301; % deg C

Hvl = Hvl_ref * ((Tc - T2)/(Tc - T1)) ^ 0.38;

fprintf("(Watson's equation) Heat of vaporization at 200degC is %.2f kJ/kg.\n", Hvl)

function y=dHdpT(p,T) % T in K, p in Pa, V in m^3
global a; 
global b; 
global R;

V=vdw(p,T); 
y=V-T.*R*V.^2./(3*p.*V.^2-2*(b*p+R*T).*V+a); 

end

function V=vdw(p,T)
global a; 
global b; 
global R;

V=p; % Matlab’s integral function calls this with a vector of p values, so cycle through them all. 
for i=1:length(p) 
    V(i)=fzero(@(V)p(i)*V^3-(b*p(i)+R*T)*V^2+a*V-a*b,R*T/p(i)); % solve the e.o.s for V given p and T! 
end 

end

function cp=cpH2Og(T) 
% Ideal-gas heat capacity (kJ/(molK)) for H2O with T in deg C 
% Source: Appendix B of Felder & Rousseau 

a=33.46e-3; 
b=0.6880e-5; 
c=0.7604e-8; 
d=-3.593e-12; 

cp=a+b*T+c*T.^2+d*T.^3; 

if(any(T>1500)||any(T<0)) 
    
    fprintf('ERROR: T=%e deg C is out of range\n', T); 
    
end

end

function cp=cpH2Ol(T) % divide by 18 to convert to kJ/kg

% Liquid heat capacity (J/(molK)) for H2O with T in K % Source: GVR

% Appendix 6

a=1.82964e1;

b=4.72118e-1;

c=-1.33878e-3;

d=1.31424e-6;

cp=a+b*T+c*T.^2+d*T.^3;

end