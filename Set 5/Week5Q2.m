% GVR 4.14

% Reducing atom matrix

F = [
    1 0 0 0 0 1;
    0 0 2 0 0 0;
    0 2 1 2 0 2;
    2 1 0 0 0 0;
    0 0 0 0 2 0;
];

b = zeros(1, 6)';

rank(F) % Rank is 5

rref(F)



