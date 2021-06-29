% GVR 4.10

x = fsolve(@(x)balance(x), [1 1 1 1])

function y=balance(x)

N3d = x(1);
nW3 = x(2);
nCH42 = x(3);
nC2H62 = x(4);

N2 = 100;

y = [
    2*1.2*(2*nCH42 + 7/2*nC2H62) - 2*0.084*N3d - 0.012*N3d - 2*0.042*N3d - nW3;
    0.79/0.21*1.2*(2*nCH42 + 7/2*nC2H62) + (N2 - nCH42 - nC2H62) - 0.862*N3d;
    nCH42 + 2*nC2H62 - 0.084*N3d - 0.012*N3d;
    4*nCH42 + 6*nC2H62 - 2*nW3;
    ]

end