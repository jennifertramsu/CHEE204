function y = speciesBalance(x)

N1 = x(1);
nO23 = x(2);
percent = x(3);

nCS22 = 35;
nSO22 = 10;
nW2 = 55
nW3 = nW2;
r1 = nCS22;
nCO23 = r1;
nSO23 = nSO22 + 2*r1;
N3 = nSO23/0.02;

y = [
    0.21*N1 - nO23 - 3*r1;
    0.79*N1 - (N3 - nO23 - nCO23 - nW3 - 0.02*N3);
    0.21*N1 - percent*3*nCS22;
]

return