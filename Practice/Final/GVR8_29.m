N1 = fsolve(@(N) (N * integral(@(T)water_cp(T), 888, 755) + 750*1000 - ), 100)

F1 = N1 / 18 * 3600 / 1000 %klb/h

p1 = 3.447e+6; % Pa
p2 = 6.895e+6; % Pa


function y=water_cp(T)

a = 3.40471e1;
b = -9.65064e-3;
c = 3.29983e-5;
d = -2.04467e-8;
e = 4.30228e-12;

y = a + b.*T + c.*T.^2 + d*T.^3 + e*T.^4;

end