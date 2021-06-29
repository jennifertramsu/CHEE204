%Modifying Page 20

F1 = 100; %kg/s

A = [
    0 0.24 0;
    0 0.01 1;
    1 -0.75 0
    ];

b = [
    0.5 * F1;
    0.5 * F1;
    0
    ];

F = A \ b;

F4 = F(3);

mE1 = 0.5 * F1;
mW1 = 0.5 * F1;

mE3 = 0.01 * F(2);
mW3 = 0.24 * F(2);
mB3 = 0.75 * F(2);

fprintf('The outlet flow rates are F3 = %.2f kg/s and F4(mE4) = %.2f kg/s.\n' , F(2), F(3))
fprintf('The inlet flow rates are F1 = 100 kg/s and F2(mB2) = %.2f kg/s.\n', F(1))
fprintf('The mass flows for stream 1 are mW1 = %.2f kg/s and mE1 = %.2f kg/s.\n', mW1, mE1)
fprintf('The mass flow rates for stream 3 are mW3 = %.2f kg/s, mE3 = %.2f kg/s and mB3 = %.2f kg/s.\n', mW3, mE3, mB3)