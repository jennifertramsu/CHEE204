% Question 1
% the purpose of the period preceding *, ^, and / is to indicate an
% element-by-element operation (ie. does not follow rules of linear
% algebra)

% Question 2
% H = CpdT -- evaluate from T  = 80C to T = 320C
T1 = 80 + 273;
T2 = 320 + 273;
H = integral(@(T) heatCapacityWater(T), T1, T2);

sprintf('Change in enthalpy is %.2f J/mol (%.2f J/g)', H, H / 18)

%steam tables

H1 = 2643.8;
H2 = 2703.7;

sprintf('From the steam tables, change in enthalpy is %.2f J/g', H2 - H1)

% Question 3


% Question 4

T = linspace(0, 100, 10);
y = zeros(1, length(T));

for i = 1:length(T)
    y(i) = heatCapacityWater(T(i) + 273);
end

plot(T, y / 18);
xlabel('Temperature (deg C)');
ylabel('Heat Capacity (kJ/(kg K))');
title('Heat Capacity vs Temperature');

temp = 0:10:100;
Hcapacities = zeros(1, length(temp));

for i = 1:11
    Hcapacities(i) = integral(@(T) heatCapacityWater(T), 273, temp(i) + 273);
end

Hcapacities / 18

% Question 5

root = fzero(@(x) func(x), 2)






