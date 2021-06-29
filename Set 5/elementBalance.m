function y = elementBalance(x) 

N1 = x(1);
nO23 = x(2);
percent = x(3); % Excess is percent - 1

nCS22 = 35;
nSO22 = 10;
nW2 = 55;
nCO23 = nCS22;
nW3 = nW2;
nSO23 = nSO22 + 2*nCS22;
N3 = nSO23/0.02;

y = [
    2*nSO22 + nW2 + 2*0.21*N1 - 2*nCO23 - 2*nSO23 - 2*nO23 - nW3;
    2*0.79*N1 - 2*(N3 - nCO23 - nO23 - 0.02*N3 - nW3);
    0.21*N1 - percent*3*nCS22;
    ]

return