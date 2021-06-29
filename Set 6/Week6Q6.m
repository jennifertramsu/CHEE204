% GVR 7.34

% Enthalpy difference beteen water at 5 bar and 50degC and steam at 5 bar
% and 400degC

%Part b

% Watson's correlation with n = 0.38

% Choose liquid water at 5 bar and 50deg C as reference

dH = integral(@(T)cpWaterLiq(T), 50 + 273, 151.8 + 273) + 37933.2 + integral(@(T)cpWaterVapour(T), 151.8 + 273, 400 + 273)

function cp=cpWaterLiq(T)

A = 1.82964e1;
B = 4.72118e-1;
C = -1.33878e-3;
D = 1.31424e-6;

cp = A + B*T + C*T.^2 + D*T.^3;

end

function cp=cpWaterVapour(T)

% J/(molK)

A = 3.40471e1;
B = -9.65064e-3;
C = 3.29983e-5
D = -2.04467e-8;
E = 4.30228e-12;

cp = A + B*T + C*T.^2 + D*T.^3 + E*T.^4;

end