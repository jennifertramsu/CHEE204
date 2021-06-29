function H = heatCapacityWater(T)

A = 3.40471;
B = -9.65064e-3;
C = 3.29983e-5;
D = -2.04467e-8;
E = 4.30228e-12;

H = A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;

return

