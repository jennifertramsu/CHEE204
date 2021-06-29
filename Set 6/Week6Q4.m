% GVR 7.8

% Given 1/3 toluene, 1/3 p-xylene, 1/3 ethylbenzene at 1 atm (101.7 kPa)

% Part a, find dewpoint temperature

    % Solve 1/3 * P/PsatToluene + 1/3*P/psatXylene + 1.3*P/psatEthyl - 1 = 0  
    
    Tdew = fsolve(@(Td)(1/3*101.7/psatToluene(Td) + 1/3*101.7/psatXylene(Td) + 1/3*101.7/psatEthyl(Td) - 1), 100);
    
    fprintf('The dew temperature is %.2f K(or %.2f C).', Tdew, Tdew - 273);
    
% Part b, find bubble temperature

    % Solve 0.6 * PsatPentene/P + 0.4 * psatHeptene/P - 1 = 0  
    
    Tbub = fsolve(@(Td)(1/3/101.7*psatToluene(Td) + 1/3/101.7*psatXylene(Td) + 1/3/101.7*psatEthyl(Td) - 1), 1000);
    
    fprintf('The bubble temperature is %.2f K(or %.2f C).\n', Tbub, Tbub - 273);
    
% Part c, find vapour fraction at 90degC

    % Solve 0.6*(1 - Kpentene)/(Kpentene*X + 1 - X) + 0.4*(1 - Kheptene)/(Kheptene*X + 1 - X)
    
    Ktoluene = psatToluene(400)/101.7;
    Kxylene = psatXylene(400)/101.7;
    Kethyl = psatEthyl(400)/101.7;
    
    vapourFraction = fsolve(@(X)(1/3*(1 - Ktoluene)/(Ktoluene*X + 1 - X) + 1/3*(1 - Kxylene)/(Kxylene*X + 1 - X) + 1/3*(1 - Kethyl)/(Kethyl*X + 1 - X)), 0.5);
    
    Ytoluene = 1/3*Ktoluene / (Ktoluene*vapourFraction + 1 - vapourFraction);
    Yxylene = 1/3*Kxylene / (Kxylene*vapourFraction + 1 - vapourFraction);
    Yethyl = 1/3*Kethyl / (Kethyl*vapourFraction + 1 - vapourFraction);
    
    fprintf("The vapour fraction at 400 K is %.3f.\n", vapourFraction)
    fprintf("The vapour fractions of toluene, p-xylene, and ethylbenzene are %.2f, %.2f, and %.2f.\n", Ytoluene, Yxylene, Yethyl)


    function p=psatToluene(T) 
    
    A = 14.2515;
    B = 3242.38;
    C = -47.1806;
    
    p=exp(A - B/(T + C));  % P in Kpa, T in K
      
    end
           
    function p=psatXylene(T)
    
    A = 14.0891;
    B = 3351.69;
    C = -57.6;
    
    p=exp(A - B/(T + C));
    
    end
    
    function p=psatEthyl(T)
    
    A = 13.9698;
    B = 3257.17;
    C = -61.0096;
    
    p=exp(A - B/(T + C));
    
    end