% GVR 3.17

y = fsolve(@(x)balance(x), [1 1 1 1 1 1 1 1 1 1])

function f=balance(x)

N2 = x(1);
nW10 = x(2);
r1 = x(3);
nCH3OH4 = x(4);
nW4 = x(5);
nHI8 = x(6);
nW8 = x(7)
N8 = x(8);
N6 = x(9);
N7 = x(10);

f = [
    nW10 - 0.05*N2;
    0.475*N2 - r1;
    1.425*N2 - nCH3OH4;
    0.525*N2 - nW4;
    nHI8 - 0.95*N2;
    0.0525*N2 - nW8;
    N8 - 1.0025*N2;
    N6 - 0.590625*N2;
    0.5743*N2 - N7;
    nCH3OH4 - 0.2*N6 - 0.18*N7;
    ];
   
end