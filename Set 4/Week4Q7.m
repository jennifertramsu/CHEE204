% GVR 4.28

M = [
    1 0 1 2 2 1 2 1 0;
    1 0 0 0 0 1 1 2 1;
    0 2 4 6 8 4 6 0 2;
    ];

rref(M)

F = [
    0.48 0.4;
    0.52 0.6;
    ];

b = [
    0.44;
    0.56;
    ];

x = fsolve(@(x)reactor(x), [0 0 0 0])

nC2 = 1/12*(x(1) - x(4) - 9/10e6 * x(1) - 1/10e6 * x(1))

N2 = nC2 + 6.94 + 9/(32*10e6)*x(1) + 1/(14*10e6)*x(1)

function r=reactor(x)

F2 = x(1);
nH24 = x(2);
nH2S4 = x(3);
mH2 = x(4);

nC1 = 4.2;
nCH43 = 0.56;
nH1 = 6.84;
nH23 = 0.44;
nN1 = 0.06;
nS1 = 0.06;

r = [
    nC1 + nCH43 - 1/12*(F2 - mH2 - 9/(10e6) * F2 - 1/(10e6) * F2) - 7/3*nH24;
    nH1 + 2*nH23 + 4*nCH43 - 2*nH24 - 4*7/3*nH24 - 2*nH2S4 - 3/4*nH24 - mH2;
    nN1 - 1/4*nH24 - 1/(14*10e6)*F2;
    nS1 - nH2S4 - 9/(32*10e6)*F2;
];

end