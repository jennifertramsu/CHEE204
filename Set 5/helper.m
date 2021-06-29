function M=helper(x)

methane = x(1);
ethane = x(2);
r1 = x(3);
r2 = x(4);
r3 = x(5);

nCO23 = 84;
nCO3 = 12;
nO23 = 42;
nN23 = 862;

M = [
    1.2*(2*methane + 7/2 * ethane) - r1 - 0.5*r2 - 0.5*r3 - nO23;;
    methane + 3*r1 + 3*r2 - 2*r3;
    ethane - 2*r1 - 2*r2 + r3;
    -nCO23 + r1;
    -nCO3 + r2;
    ]

return