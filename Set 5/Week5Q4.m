
% Plotting enthalpy change vs exit temperature

% Range: 25-300

format bank

temp = linspace(25, 300, 50);

%deltaH = enthalpyChange(temp, 0)
deltaH = zeros(1, 50);

for i=1:50
    
   deltaH(i) = enthalpyChange(temp(i), 0);

end

plot(temp, deltaH)
grid on;
title('Change in Enthalpy vs Exit Temperature')
xlabel('Exit Temperature (degC)')
ylabel('Change in Enthalpy (kW)')

T2 = fsolve(@(T)enthalpyChange(T, 0.5), 100)

hold on
plot(T2, 0,'o', 'MarkerSize', 10, 'LineWidth', 3)
text(T2, 0.1, 'T2 = 193.57 C')
text(T2 - 80, -0.1, 'Temperature where heat transfer is 0.5 kW')
hold off

print("Week5Q4 Plot", '-dpdf')
print("Week5Q4 Plot", '-depsc')

function cp=cpH2Og(T) 

% Ideal-gas heat capacity (kJ/(molK)) for H2O with T in deg C 
% Source: Appendix B of Felder & Rousseau 
a=33.46e-3; 
b=0.6880e-5; 
c=0.7604e-8; 
d=-3.593e-12; 
cp=a+b.*T+c.*T.^2+d*T.^3; 

if(any(T>1500)||any(T<0)) 
    fprintf('ERROR: T=%e deg C is out of range\n', T); 
end

end

function y=enthalpyChange(T2, Q) 
T1=100;
Tref=25;
N1=10e3/18/3600; 
y=N1*(integral(@(T)cpH2Og(T),Tref,T1)-integral(@(T)cpH2Og(T),Tref,T2)) + Q; 
end