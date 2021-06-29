% GVR 4.24
format long
x = fsolve(@(x)speciesBalance(x), [0.5 0.5])

function f=speciesBalance(x) 

N3 = x(1);
r2 = x(2);

f = [
    0.21*20/0.2*0.015*N3 - 0.03*N3 - 5*0.015*N3 - 2*r2;
    0.79*20/0.2*0.015*N3 - (N3 - 0.03*N3 - 0.06*N3 - 0.09*N3 - 0.015*N3 - 2*r2) - r2;
    ];

end