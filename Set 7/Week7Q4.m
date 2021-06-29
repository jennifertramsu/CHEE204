% GVR 8.14

clc
clear
format bank

% Heat of reaction at 300 K

% part b

Hrxn_ref = -10000; % cal/gmol
Tref = 300; % K

T = linspace(300, 1000);
H = zeros(1, length(T));

for i=1:length(T)
    
    H(i) = Hrxn(T(i));
    
end

plot(T, H)
title("Heat of Reaction (cal/gmol) vs Temperature (K)")
xlabel("Temperature (K)")
ylabel("Heat of Reaction (cal/gmol)")

% part c

T = 500; % K

Hrxn500 = Hrxn_ref + integral(@(T)cpC(T), Tref, T) - 2*(integral(@(T)cpA(T), Tref, 400)+ 928 + 10*(500 - 400)) - integral(@(T)cpB(T), Tref, T)

% Part a

function y=Hrxn(T)

Hrxn_ref = -10000; % cal/gmol
Tref = 300; % K

y = Hrxn_ref - integral(@(T)cpC(T), Tref, T) + 2*integral(@(T)cpA(T), Tref, T) + integral(@(T)cpB(T), Tref, T);

end

function cp=cpA(T) % cal/(gmol K)

cp = 16 - (1.5e3)./T;

end

function cp=cpB(T) % cal/(gmol K)

cp = 11 - (0.5e3)./T;

end

function cp=cpC(T) % cal/(gmol K)

cp = 25 - (1.0e3)./T;

end