F2 = 16.3; %mol / h

A = [
    0.95 -(1 - 0.1);
    0.05 -0.1
    ];

b = [
    0;
    -F2
    ];

F = A \ b

fprintf('The process stream has flow rate %.2f mol/h.\n', F(1))