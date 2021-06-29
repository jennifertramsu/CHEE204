function distillationcolumn
    
    clc;
 
    % N-stage distillation column for M species obeying Raoult's law.
    % R.J.Hill,  McGill University, 2021.
    % This is a teaching tool, provided with no warranty whatsoever.

    close all
    clear all

    %set(0,'defaultaxesfontsize',20);
    %set(0,'defaulttextfontsize',20);
    %set(0,'defaultaxeslinewidth',2);
    %set(0,'defaultlinelinewidth',2);
    
    global Tref;

    M=3; % number of species
    N=10; % number of stages
    p=103.3; % operating pressure kPa (assumed uniform)
    Tf=273.15+25; % feed temperature K
    Tref=Tf;
    Nf=1; % feed flow mol / s ?
    zf=[1/3 1/3 1/3]; % feed mole fractions
    F = 4; % feed tray % VARIABLE
    Xf=0; % vapor fraction of the feed
    
    Qr = 250; 
    Qc = -232.6;
    
    R = abs(Qc/Qr);
    %profit = []
    
    %for i=1:10
    
        %caseid='A';

        X0=ones(3*N+M*N,1); % a very crude initial guess (scale temperatures with Tref to solve)

        X=fsolve(@(X)mydistillationcol(X,M,N,p,Tf,Nf,zf,Xf,F(i),Qr,Qc,Tref),X0);

            % STREAM VARIABLES AT EACH STAGE (1-N)

        [T, Nl, Nv, x, K, y]=getvarsrommaster(X,Tref,p,N,M);

        x_bot=x(1,:); % liquid mole fractions at the bottom
        y_top=y(N,:); % vapor mole fractions at the top
        %Nl(1) + Nv(N) - Nf; % check overall balance with unit molar feed

        stages=linspace(1,N,N);

        % OPTIMIZATION======================================================


            % ADDITIONAL HEAT DUTIES

        Q_condenser_top = fsolve(@(Q) condenseTop(Q, T(N), Nv(N), y_top), 0); % kJ
        Q_cooler_bottom = fsolve(@(Q) coolBottom(Q, T(1), Nl(1), x_bot), 0); % kJ

            % COST

                % FEED COST

                MW1 = 58.08; % g/mol
                MW2 = 78.11; % g/mol
                MW3 = 92.14; % g/mol

                feed_price = 0.6; % USD / kg
                feed_weight = Nf*zf(1)*(MW1 + MW2 + MW3) / 1000; % kg
                feed_cost = feed_price * feed_weight; %USD

                % HEAT EXCHANGER COSTS

                X_fraction = Nv(N)/Nf; 

                % ALL LATENT HEATS

                dHvl1=29087.2/1000; % kJ/mol
                dHvl2=30763.4/1000; % kJ/mol
                dHvl3=33460.6/1000; % kJ/mol

                dHvl_total = zf(1)*dHvl1 + zf(2)*dHvl2 + zf(3)*dHvl3;
                heat_duty_cost = 0.02/3600; % USD / kWh

                % FROM TUTORIAL

                cost_cooler = abs(heat_duty_cost * Q_cooler_bottom); % USD
                cost_condense_vapour = abs(heat_duty_cost * Q_condenser_top); % USD


                %cost_total = (X_fraction*dHvl_total)/(1-R) * (1+R) * heat_duty_cost...
                 %                       + cost_cooler + cost_condense_vapour ...
                  %                      + feed_cost

            % REVENUE

                % RATES

                toluene_rate = 907 / 907.185; % USD / kg
                benzene_rate = 1026 / 907.185; % USD / kg
                acetone_rate = 2000 / 907.185; % USD / kg

                % TOP REVENUE

                if y_top(1) >= 0.99
                    rev = Nv(N)*MW1/1000 * acetone_rate

                elseif y_top(2) >= 0.99
                    rev = Nv(N)*MW2/1000 * benzene_rate

                elseif y_top(3) >= 0.99
                    rev = Nv(N)*MW3/1000 * toluene_rate       
                
                elseif x_bot(1) >= 0.99
                    rev = Nl(1)*MW1/1000 * acetone_rate

                elseif x_bot(2) >= 0.99
                    rev = Nl(1)*MW2/1000 * benzene_rate

                elseif x_bot(3) >= 0.99
                    rev = Nl(1)*MW3/1000 * toluene_rate

                else
                    rev = 0 % NOT PURE ENOUGH
                end
               

            % PROFIT

                profit(i) = rev - cost_total

                %{
            % FIGURE

                % LIQUID FRACTIONS VS STAGE
                        figure(4)
                        plot(stages,x(:,1),'-o'); hold on
                        plot(stages,x(:,2),'-o')
                        plot(stages,x(:,3),'-o')
                        ylabel('x_i')
                        xlabel('n')
                        xlim([1,N])
                        legend('Acetone', 'Benzene', 'Toluene')
                        filename=['distillationcolumn4',caseid];
                        print('-depsc',char(filename));

                % VAPOUR FRACTIONS VS STAGE
                        figure(5)
                        plot(stages,y(:,1),'-o'); hold on
                        plot(stages,y(:,2),'-o')
                        plot(stages,y(:,3),'-o')
                        ylabel('y_i')
                        xlabel('n')
                        xlim([1,N])
                        legend('Acetone', 'Benzene', 'Toluene')
                        filename=['distillationcolumn5',caseid];
                        print('-depsc',char(filename));

        %}
                %plot(F, profit)

    end
    
    plot(F, profit)
    grid on
    
    
end

% commented out figures

%{

    % TEMPERATURE VS STAGE
figure(1)
plot(stages,T-273.15,'-o')
ylabel('T (\circC)')
xlabel('n')
xlim([1,N])
filename=['distillationcolumn1',caseid];
print('-depsc',char(filename));

    % Nl VS STAGE
figure(2)
plot(stages,Nl,'-o')
ylabel('N_l')
xlabel('n')
xlim([1,N])
filename=['distillationcolumn2',caseid];
print('-depsc',char(filename));

    % Nv VS STAGE
figure(3)
plot(stages,Nv,'-o')
ylabel('N_v')
xlabel('n')
xlim([1,N])
filename=['distillationcolumn3',caseid];
print('-depsc',char(filename));

    % VAPOUR FRACTIONS VS LIQUID FRACTIONS
figure(6)
plot(x(:,1),y(:,1),'-o'); hold on
plot(x(:,2),y(:,2),'-o')
plot(x(:,3),y(:,3),'-o')
ylabel('y_i')
xlabel('x_i')
filename=['distillationcolumn6',caseid];
print('-depsc',char(filename));
return
%}

% OPTIMIZATION CODE ==================================

    % top and bottom heat exchangers

function f=condenseTop(Q, T, Nv, y) % kJ

y1 = y(1);
y2 = y(2);
y3 = y(3);

% dHv(T, i) - vapour enthalpies, liquid reference

f = Q + Nv*(y1*dHv(T, 1) + y2*dHv(T, 2) + y3*dHv(T, 3)); 

end

function f=coolBottom(Q, T, Nl, x)

x1 = x(1);
x2 = x(2);
x3 = x(3);

% dHl(T, i) = liquid enthalpies, liquid reference

f = Q + Nl*(x1*dHl(T, 1) + x2*dHl(T, 2) + x3*dHl(T, 3)); 

end % kJ


%=====================================================



% SOLVES THE DISTILLATION COLUMN

function f=mydistillationcol(X,M,N,p,Tf,Nf,zf,Xf,F,Qr,Qc,Tref)
    % X = vector of independent variables X=[Ts||Nls|Nvs|x1,x2,...]
    % N = number of stages (inc reboiler and condenser)
    % M = number of species
    % p = pressure kPa
    % Tf = feed temperature K
    % Nf = feed flow mol/unit time
    % zf = vector of feed mole fractions
    % Xf = feed vapor fraction (assumed vapor and liquid have the same composition---not at equilibrium!)
    % F = feed stage
    % Qr = reboiler heat duty (>0)
    % Qc = condenser heat duty (<0)
    % Tref = reference temperature (to scale temperatures and as reference state for enthalpies)

    % convert master X to stream variables for each stage...
    [T, Nl, Nv, x, K, y]=getvarsrommaster(X,Tref,p,N,M);

    %============================

    % energy balances (skip condenser and reboiler)    
    for n=2:N-1
        f(n)=0;        
        for i=1:M
            f(n)=f(n)+Nl(n+1)*x(n+1,i)*dHl(T(n+1),i)+...
                      Nv(n-1)*y(n-1,i)*dHv(T(n-1),i)-...
                      Nl(n)*x(n,i)*dHl(T(n),i)-...
                      Nv(n)*y(n,i)*dHv(T(n),i);
        end
    end

    % energy balances reboiler
    for n=1
      f(n)=Qr;        
      for i=1:M
            f(n)=f(n)+...
                      Nl(n+1)*x(n+1,i)*dHl(T(n+1),i)-...
                      Nl(n)*x(n,i)*dHl(T(n),i)-...
                      Nv(n)*y(n,i)*dHv(T(n),i);
        end
    end

    % energy balances condenser    
    for n=N
        f(n)=Qc;        
        for i=1:M
            f(n)=f(n)+...
                      Nv(n-1)*y(n-1,i)*dHv(T(n-1),i)-...
                      Nl(n)*x(n,i)*dHl(T(n),i)-...
                      Nv(n)*y(n,i)*dHv(T(n),i);
        end
    end

    %  energy balances feed (add enthalpy of feed stream)    
    for n=F
        for i=1:M
            f(n)=f(n)+Nf*zf(i)*Xf*dHv(Tf,i)+...
                      Nf*zf(i)*(1-Xf)*dHl(Tf,i);
        end
    end

    %============================

    % species balances
    for n=2:N-1
        for i=1:M
        e=N+(n-1)*M+i;
            f(e)=Nl(n+1)*x(n+1,i)+...
                 Nv(n-1)*y(n-1,i)-...
                 Nl(n)*x(n,i)-...
                 Nv(n)*y(n,i);
        end
    end

    % species balances reboiler
    for n=1
        for i=1:M
            e=N+(n-1)*M+i;
            f(e)=Nl(n+1)*x(n+1,i)-Nl(n)*x(n,i)-Nv(n)*y(n,i);
        end
    end

    % species balances condenser
    for n=N
        for i=1:M
            e=N+(n-1)*M+i;
            f(e)=Nv(n-1)*y(n-1,i)-Nl(n)*x(n,i)-Nv(n)*y(n,i);
        end
    end

    % species balances feed
    for n=F
        for i=1:M
            e=N+(n-1)*M+i;
            f(e)=f(e)+Nf*zf(i);
        end
    end
    
 %============================

    % sum x_i = 1 for each stage
    for n=1:N
        e=N+N*M+n;
        f(e)=1;    
        for i=1:M
            f(e)=f(e)-x(n,i);
        end
    end

    % sum y_i = 1 for each stage
    for n=1:N
        e=2*N+N*M+n;
        f(e)=1;    
        for i=1:M
            f(e)=f(e)-y(n,i);
        end
    end

%============================

end


% GIVEN


function [T, Nl, Nv, x, K, y]=getvarsrommaster(X,Tref,p,N,M)
    for n=1:N
        T(n)=X(n)*Tref;
        Nl(n)=X(N+n);
        Nv(n)=X(2*N+n);
        for i=1:M
            x(n,i)=X(3*N+(n-1)*M+i);
            K(n,i)=Keqm(T(n),p,i);
            y(n,i)=x(n,i)*K(n,i);
        end
    end
end



% RETURNS PSATi/P

function k=Keqm(T,p,i)
    if i==1
        k=psat1(T)/p;
    end
    if i==2
        k=psat2(T)/p;
    end
    if i==3
        k=psat3(T)/p;
    end
end




%=====================================================




% CHARACTERISTIC SPECIES DATA



% ACETONE PSAT
function p=psat1(T)
% acetone
ACOEFF=14.7171;
BCOEFF=2975.95;
CCOEFF=-34.5228;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
end

% BENZENE PSAT
function p=psat2(T)
% benzene
ACOEFF=14.1603;
BCOEFF=2948.78;
CCOEFF=-44.5633;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
end

% TOLUENE PSAT
function p=psat3(T)
% toluene
ACOEFF=14.2515;
BCOEFF=3242.38;
CCOEFF=-47.1806;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
end




% LIQUID ENTHALPIES (kJ/mol) FROM TREF TO T (K), LIQUID REFERENCE

function h=dHl(T,i)
% liquid enthalpies kJ/mol (T in K)
global Tref;
if i==1
    h=integral(@(T)cpacetonel(T),Tref,T);
end
if i==2
    h=integral(@(T)cpbenzenel(T),Tref,T);
end
if i==3
    h=integral(@(T)cptoluenel(T),Tref,T);
end
end

% VAPOUR ENTHALPIES (kJ/mol) FROM TREF TO T (K), LIQUID REFERENCE
function h=dHv(T,i)
% vapor enthalpies kJ/mol (T in K)
global Tref;
Tb1=329.281; % acetone
Tb2=353.261; % benzene
Tb3=383.786; % toluene
dHvl1=29087.2/1000; % kJ/mol
dHvl2=30763.4/1000; % kJ/mol
dHvl3=33460.6/1000; % kJ/mol
if i==1
    h=integral(@(T)cpacetonel(T),Tref,Tb1)+dHvl1+integral(@(T)cpacetonev(T),Tb1,T);
end
if i==2
    h=integral(@(T)cpbenzenel(T),Tref,Tb2)+dHvl2+integral(@(T)cpbenzenev(T),Tb2,T);
end
if i==3
    h=integral(@(T)cptoluenel(T),Tref,Tb3)+dHvl3+integral(@(T)cptoluenev(T),Tb3,T);
end
end






% ALL HEAT CAPACITIES

function cp=cpbenzenel(T)
% liquid heat capacity (kJ/molK) for benzene with T in deg C
T=T-273.15;
a=126.5e-3;
b=23.4e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
end

function cp=cpbenzenev(T)
% vapor heat capacity (kJ/molK) for benzene with T in deg C
T=T-273.15;
a=74.06e-3;
b=32.95e-5;
c=-25.20e-8;
d=77.57e-12;
cp=a+b*T+c*T.^2+d*T.^3;
end

function cp=cptoluenel(T)
% liquid heat capacity (kJ/molK) for toluene with T in deg C
T=T-273.15;
a=148.8e-3;
b=32.4e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
end

function cp=cptoluenev(T)
% vapor heat capacity (kJ/molK) for toluene with T in deg C
T=T-273.15;
a=94.18e-3;
b=38.00e-5;
c=-27.86e-8;
d=80.33e-12;
cp=a+b*T+c*T.^2+d*T.^3;
end

function cp=cpacetonel(T)
% liquid heat capacity (kJ/molK) for acetone with T in deg C
T=T-273.15;
a=123.0e-3;
b=18.6e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
end

function cp=cpacetonev(T)
% vapor heat capacity (kJ/molK) for acetone with T in deg C
T=T-273.15;
a=71.96e-3;
b=20.10e-5;
c=-12.78e-8;
d=34.76e-12;
cp=a+b*T+c*T.^2+d*T.^3;
end




