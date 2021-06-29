%2.16
%Part a)

%DOF = -1, so system is overspecified
% Proof by contradiction
fprintf('Part a.\n')

F4 = 10000; %lbm/h

A_a = [
    0.6 0.2 0.2;
    0.2 0.6 0;
    0.2 0 0.6;
    0 0.2 0.2
    ];

b_a = [
    0.25;
    0.25;
    0.25;
    0.25
    ] * F4;

F_a = A_a \ b_a

fprintf('The mass flow rates of the feed alloy are calculated to be F1 = %.2f lbm/h, F2 = %.2f lbm/h, and F3 = %.2f lbm/h.\n', F_a)
fprintf('However, these values contradict the total mass balance, where their sum, %.2f lbm/h, is not equal to %.2f lbm/h.\n', F_a(1) + F_a(2) + F_a(3), F4)

%Part b)

fprintf('\nPart b.\n')

F5 = 10000; %lbm/h

A_b = [
    0.6 0.2 0.2 0;
    0.2 0.6 0 0.2;
    0.2 0 0.6 0.2;
    0 0.2 0.2 0.6;
    ];

b_b = [
    0.25*F5;
    0.25*F5;
    0.25*F5;
    0.25*F5
    ];

F_b = A_b \ b_b

fprintf('The mass flow rates of the feed alloys are F1 = F2 = F3 = F4 = %.2f lbm/h.\n', F_b(1))

%Part c
fprintf('\nPart c.\n')

A_c = [
    0.6 0.2 0.2 0;
    0.2 0.6 0 -F5;
    0.2 0 0.6 -F5;
    0 0.2 0.2 2*F5;
    ];

b_c = [
    0.4*F5;
    0;
    0;
    (1 - 0.4)*F5
    ];

F_c = A_c \ b_c

xb5 = F_c(4);
xc5 = xb5;
xd5 = 1 - 0.4 - 2*xb5;

fprintf('The mass flow rates of the feed alloys are F1 = %.2f lbm/h, F2 = %.2f lbm/h,\n', F_c(1), F_c(2))
fprintf('and F3 = %.2f lbm/h.\n', F_c(3))

fprintf('The composition of the alloy is 0.40 A, %.2f B, %.2f C, and %.2f D.\n', xb5, xc5, xd5)