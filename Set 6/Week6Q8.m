% GVR 7.38

P1 = 200; %kpa
P2 = 55; %kpa

zE = 0.7;
zW = 0.3;

kE = EsatP(70 + 273.15)/P2;
kW = WsatP(70 + 273.15)/P2;

x=fsolve(@(X) (zE*(1 - kE)/(kE*X + 1 - X) + zW*(1 - kW)/(kW*X + 1 - X)), 0.5)

xW = zW/(kW*x + 1 - x)
xE = 1 - xW
yW = zW*kW/(kW*x + 1 - x)
yE = 1 - yW

dHcondW = -41387.4; % J/mol
dHcondE = -38577.3; %J/mol

Tref = 110 + 273.15; %K

heatTransfer = -((zE * integral(@(T)EvapHC(T), Tref, 110 + 273.15)) + (zW * integral(@(T)(WvapHC(T)), Tref, 110 + 273.15)) - x*(yE * integral(@(T)(EvapHC(T)), Tref, 70 + 273.15) + yW * integral(@(T)(WvapHC(T)), Tref, 70 + 273.15)) - (1 - x) * (xE * (integral(@(T)(EvapHC(T)), Tref, 351.481)) + dHcondE + integral(@(T) (EliqHC(T)), 351.481, 70 + 273.15)) + xW * (integral(@(T)(WvapHC(T)), Tref, 356.85) + dHcondW + integral(@(T) (WliqHC(T)), 356.85, 70 + 273.15)))

function y=EvapHC(T)

a = 1.76907e1;
b = 1.49532e-1;
c = 8.94815e-5;
d = -1.97384e-7;
e = 8.31747e-11;

y = a + b.*T + c.*T.^2 + d.*T.^3 + e.*T.^4;

end

function y=WvapHC(T)

a = 3.40471e1;
b = -9.65064e-3;
c = 3.29983e-5;
d = -2.04467e-8;
e = 4.30228e-12;

y = a + b.*T + c.*T.^2 + d.*T.^3 + e.*T.^4;

end

function y=EliqHC(T)

a = -3.25137e2;
b = 4.13787e0;
c = -1.40307e-2;
d = 1.70354e-5;

y = a + b.*T + c.*T.^2 + d.*T.^3;

end

function y=WliqHC(T)

a = 1.82964e1;
b = 4.72118e-1;
c = -1.33878e-3;
d = 1.31424e-6;

y = a + b.*T + c.*T.^2 + d.*T.^3;

end


function y=EsatP(T)

A = 16.1952;
B = 3423.53;
C = -55.7152;

y = exp(A - B./(T + C))

end

function y=WsatP(T)

A = 16.5362;
B = 3985.44;
C = -38.9974;

y = exp(A - B./(T + C))

end