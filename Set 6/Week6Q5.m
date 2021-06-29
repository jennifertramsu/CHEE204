% GVR 7.26

% Convert 1 lb CO2 to lbmol

CO2lbmol = 1 / (12 + 2*16);

q = CO2lbmol* integral(@(T)heatcapCO2(T), 60, 400)

fprintf("Heat required is %.2f Btu.\n", q)

function cp=heatcapCO2(T) 

cp = 9.00 + 2.71e-3.*T - 0.256e-6.*T.^2;

end