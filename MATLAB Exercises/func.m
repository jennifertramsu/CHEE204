function y = func(x)

a = 1;
b = 2;
c = 3;
d = 4;

y = a ./(x - 10) + b .* x .^ 2 - 1 ./(c .* log(x)) - d.*x;

return