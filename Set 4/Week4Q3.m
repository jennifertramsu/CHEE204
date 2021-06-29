% Reactor + Separator + Splitter
clc
format bank

N1 = 100; % mol/h
nI1 = N1 / 12;
nR1 = nI1 * 11;

nI8 = nI1;
nR8 = 0.5 * nR1;
N8 = nI8 + nR8;

nP6 = 0.25 * 0.8 * nR1;
N6 = nP6;
r1 = nP6;

nW4 = nR1 - nR8 - r1;
r2 = nW4 - r1;
nB4 = r2;

N4 = nW4 + nB4;

wR8 = nR8 / N8;
wI8 = 1 - wR8;

wR2 = wR8;
wR7 = wR2;
wI2 = wI8;
wI7 = wI2;

F1 = [
    0.85 -wR2;
    0.15 -wI2;
    ];

b1 = [nR8 + 2*r1 + r2; nI8]

x1 = F1 \ b1 % N3, N2

N3 = x1(1);
N2 = x1(2);

N7 = 0.15 / wI7 * N3;
N5 = N7 - N6;

% Gotta be honest, I know these values are weird
% They are also weirdly sensitive to rounding
% I'm not really sure why :(
