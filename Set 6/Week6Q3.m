% GVR 7.6

% Given 0.6 1-pentene, 0.4 1-heptene at 2 bar (200 kPa)

% Part a, find dewpoint temperature

    % Solve 0.6 * P/PsatPentene + 0.4 * P/psatHeptene - 1 = 0  
    
    Tdew = fsolve(@(Td)(0.6*200/psatPentene(Td) + 0.4*200/psatHeptene(Td) - 1), 100);
    
    fprintf('The dew temperature is %.2f K(or %.2f C).', Tdew, Tdew - 273);
    
% Part b, find bubble temperature

    % Solve 0.6 * PsatPentene/P + 0.4 * psatHeptene/P - 1 = 0  
    
    Tbub = fsolve(@(Tb)(0.6/200*psatPentene(Tb) + 0.4/200 * psatHeptene(Tb) - 1), 1000);
    
    fprintf('The bubble temperature is %.2f K(or %.2f C).\n', Tbub, Tbub - 273);
    
% Part c, find vapour fraction at 90degC

    % Solve 0.6*(1 - Kpentene)/(Kpentene*X + 1 - X) + 0.4*(1 - Kheptene)/(Kheptene*X + 1 - X)
    
    Kpentene = psatPentene(90 + 273)/200;
    Kheptene = psatHeptene(90 + 273)/200;
    
    vapourFraction = fsolve(@(X)(0.6*(1 - Kpentene)/(Kpentene*X + 1 - X) + 0.4*(1 - Kheptene)/(Kheptene*X + 1 - X)), 0.5);
    
    fprintf("The vapour fraction at 90degC is %.3f.\n", vapourFraction)


    function p=psatPentene(T) 
    
    A = 13.7564;
    B = 2409.11;
    C = -39.4834;
    
    p=exp(A - B/(T + C));  % P in Kpa, T in K
      
    end
           
    function p=psatHeptene(T)
    
    A = 13.8747;
    B = 2895.90;
    C = -53.9388;
    
    p=exp(A - B/(T + C));
    
    end