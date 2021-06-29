% GVR 3.20

t = fsolve(@(x)balance(x), [1 1 1 1 1 1 1])
function y=balance(x)

r1 = x(1);
r2 = x(2);
r3 = x(3);
nC3 = x(4);
nD3 = x(5);
nE3 = x(6);
nF3 = x(7);

y = [
    2*r1 + r2 - 160;
    90 - r1 - 2*r3;
    -nC3 - r3 + 2*r2;
    -nD3 - r2 + 2*r1;
    r1 + r2 - nE3;
    nF3 - 2*r3;
    300 - nC3 - nD3 - nE3 - nF3;
    ]

end