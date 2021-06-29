function distillation_breakeven
    % N-stage distillation column for M species obeying Raoult's law.
    % R.J.Hill,  McGill University, 2021.
    % This is a teaching tool, provided with no warranty whatsoever.

    clc 
    
    close all
    clear all
    format long

    set(0,'defaultaxesfontsize',20);
    set(0,'defaulttextfontsize',20);
    set(0,'defaultaxeslinewidth',2);
    set(0,'defaultlinelinewidth',2);

    global Tref;
    caseid='A';
    M=3; % number of species
    N=10; % number of stages
    p=101.325; % operating pressure kPa (assumed uniform)
    Tf=273.15+25; % feed temperature K
    Tref=Tf; 
    Nf=1; % feed flow mol
    zf=[1/3 1/3 1/3]; % feed mole fractions
    Mw1 = 58.08/1000 ; Mw2 = 78.11/1000; Mw3 = 92.14/1000; % kg/mol
    F=4; % feed tray
    Xf=0; % vapor fraction of the feed
    X0=ones(3*N+M*N,1); % a very crude initial guess (scale tempertaures with Tref to solve)


    liquid_price = 0.6; % $/kg
    cost_heating = 0.02; % $/kWh
    cost_of_feed = (Nf*(zf(1)*Mw1 + zf(2)*Mw2 + zf(3)*Mw3)) * liquid_price; % The cost of the feed stream 
    price_per_kJ = cost_heating/3600; % $/kJ 

    Qr = 250; % At R = 0.87, and F = 5
    Qc = -232.6;
    R = abs(Qc/Qr);
    profit_on_sold_liquid = [];
    y_acetone = [];
    % R = [0.5:0.05:0.99]
    % R = [0.84:0.1:0.85]; R = flip(R);

    for i = 1:10
            profit_on_sold_liquid
            %Qc = -Qr*R(i);

            % delta_hvl = zf(1) * 29087.2/1000 + zf(2)* 30763.4/1000 +zf(3)* 33460.6/1000
            % abs_qr_minus_qc = Qr - abs(Qc)


            X=fsolve(@(X)mydistillationcol(X,M,N,p, Tf,Nf,zf,Xf,F,Qr,Qc,Tref),X0);

            [T, Nl, Nv, x, K, y]=getvarsrommaster(X,Tref,p,N,M);

            x_bot = x(1,:); % liquid mole fractions at the bottom
            y_top=y(N,:); % vapor mole fractions at the top
            % Nl(1)+Nv(N)-1; % check overall balance with unit molar feed

            X_vapor_fraction = Nv(N)/Nf;
            deltaHvl = zf(1) * 29087.2/1000 + zf(2)* 30763.4/1000 +zf(3)* 33460.6/1000;

            fprintf('y(acetone) is %f and x(acetone) is %f\n', y_top(1), x_bot(1))
            fprintf('y(benzene) is %f and x(benzene) is %f\n', y_top(2), x_bot(2))
            fprintf('y(toluene) is %f and x(toluene) is %f\n', y_top(3), x_bot(3))

            cost_acetone = linspace(2000, 3000, 10)
            cost_benzene = 1026;
            cost_toluene = 907;

            if y_top(1) >= 0.99
                total_revenue_function = Nv(N)*y_top(1)*Mw1*(cost_acetone(i)/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            elseif y_top(2) >= 0.99
                total_revenue_function = Nv(N)*y_top(2)*Mw2*(cost_benzene/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            elseif y_top(3) >= 0.99
                total_revenue_function = Nv(N)*y_top(3)*Mw3*(cost_toluene/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            elseif x_bot(1) >= 0.99
                total_revenue_function = Nl(1)*x_bot(1)*Mw1*(cost_acetone(i)/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            elseif x_bot(2) >= 0.99
                total_revenue_function = Nl(1)*x_bot(1)*Mw2*(cost_benzene/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            elseif x_bot(3) >= 0.99
                total_revenue_function = Nl(1)*x_bot(1)*Mw3*(cost_toluene/907.185)
                % This calculates the total revenues from selling the top and bottom stream
            else
                profit_on_sold_liquid(i) = NaN; 
                y_acetone(i) = NaN;  % If the seperation of acetone isn't above 99%, we discard these points.
                continue
            end

            Q_condense_final = fsolve(@(Q) condensevapor(Q, Nv(N), T(end), y_top), 10); % Solve for the heat required for condensing vapor
            Q_cooler_final = fsolve(@(Q) coolerliquid(Q, Nl(1), T(1), x_bot),10); % Solve for the heat required to cool the liquid

            heating_costs = (Nf*(Qr + abs(Qc)) + abs(Q_condense_final + Q_cooler_final)) * price_per_kJ;
            
            total_cost_function = heating_costs + cost_of_feed;
            
            profit_on_sold_liquid(i) = total_revenue_function - total_cost_function % Simply determines the difference in revenue and cost (i.e. Profit)        

    end

    breakeven = total_cost_function / (Nv(N)*y_top(1)*Mw1) * 907.185
    
    figure(1)
    plot(cost_acetone, profit_on_sold_liquid)
    grid on
    
    hold on
    title("Profit vs Acetone Market Rate")
    xlabel("Cost of Acetone (USD/ton)")
    ylabel("Profit (USD)")
    plot(breakeven, 0, 'o', 'MarkerSize', 10, 'LineWidth', 3)
    % t = sprintf("Breakeven cost is %.2f USD/ton", breakeven)
    % text(breakeven, 0, t)
    hold off
    
    % plot(R, profit_on_sold_liquid)
    % table(profit_on_sold_liquid', y_acetone')
    
    figure(2)

return


% stages=linspace(1,N,N);
% 
% 
% figure(1)
% plot(stages,T-273.15,'-o')
% ylabel('T (\circC)')
% xlabel('n')
% xlim([1,N])
% filename=['distillationcolumn1',caseid];
% print('-depsc',char(filename));
% 
% figure(2)
% plot(stages,Nl,'-o')
% ylabel('N_l')
% xlabel('n')
% xlim([1,N])
% filename=['distillationcolumn2',caseid];
% print('-depsc',char(filename));
% 
% figure(3)
% plot(stages,Nv,'-o')
% ylabel('N_v')
% xlabel('n')
% xlim([1,N])
% filename=['distillationcolumn3',caseid];
% print('-depsc',char(filename));
% 
% figure(4)
% plot(stages,x(:,1),'-o'); hold on
% plot(stages,x(:,2),'-o')
% plot(stages,x(:,3),'-o')
% ylabel('x_i')
% xlabel('n')
% xlim([1,N])
% filename=['distillationcolumn4',caseid];
% print('-depsc',char(filename));
% 
% figure(5)
% plot(stages,y(:,1),'-o'); hold on
% plot(stages,y(:,2),'-o')
% plot(stages,y(:,3),'-o')
% ylabel('y_i')
% xlabel('n')
% xlim([1,N])
% filename=['distillationcolumn5',caseid];
% print('-depsc',char(filename));
% 
% figure(6)
% plot(x(:,1),y(:,1),'-o'); hold on
% plot(x(:,2),y(:,2),'-o')
% plot(x(:,3),y(:,3),'-o')
% ylabel('y_i')
% xlabel('x_i')
% filename=['distillationcolumn6',caseid];
% print('-depsc',char(filename));
% return
% 
function f = condensevapor(Q, N1, T, y)

dHvl1=29087.2/1000; % kJ/mol
dHvl2=30763.4/1000; % kJ/mol
dHvl3=33460.6/1000; % kJ/mol

yacetone = y(1);
ybenzene = y(2);
ytoluene = y(3);

f = Q + N1*(yacetone*dHv(T, 1) + ybenzene*dHv(T, 2) + ytoluene*dHv(T, 3) - (yacetone*dHvl1 + ybenzene*dHvl2 + ytoluene*dHvl3));
return


function f = coolerliquid(Q, Nl,T, x)

xacetone = x(1);
xbenzene = x(2);
xtoluene = x(3);

f = Q + Nl*(xacetone*dHl(T,1) + xbenzene*dHl(T,2) + xtoluene*dHl(T,3));
return


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

return


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
return


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
return


%==================================================
function p=psat1(T)
% acetone
ACOEFF=14.7171;
BCOEFF=2975.95;
CCOEFF=-34.5228;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
return
%==================================================
function p=psat2(T)
% benzene
ACOEFF=14.1603;
BCOEFF=2948.78;
CCOEFF=-44.5633;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
return
%==================================================
function p=psat3(T)
% toluene
ACOEFF=14.2515;
BCOEFF=3242.38;
CCOEFF=-47.1806;
p=exp(ACOEFF-BCOEFF./(T+CCOEFF));
return
%==================================================

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
return


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
return

function cp=cpbenzenel(T)
% liquid heat capacity (kJ/molK) for benzene with T in deg C
T=T-273.15;
a=126.5e-3;
b=23.4e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
return

function cp=cpbenzenev(T)
% vapor heat capacity (kJ/molK) for benzene with T in deg C
T=T-273.15;
a=74.06e-3;
b=32.95e-5;
c=-25.20e-8;
d=77.57e-12;
cp=a+b*T+c*T.^2+d*T.^3;
return

function cp=cptoluenel(T)
% liquid heat capacity (kJ/molK) for toluene with T in deg C
T=T-273.15;
a=148.8e-3;
b=32.4e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
return

function cp=cptoluenev(T)
% vapor heat capacity (kJ/molK) for toluene with T in deg C
T=T-273.15;
a=94.18e-3;
b=38.00e-5;
c=-27.86e-8;
d=80.33e-12;
cp=a+b*T+c*T.^2+d*T.^3;
return

function cp=cpacetonel(T)
% liquid heat capacity (kJ/molK) for acetone with T in deg C
T=T-273.15;
a=123.0e-3;
b=18.6e-5;
c=0.0;
d=0.0;
cp=a+b*T+c*T.^2+d*T.^3;
return

function cp=cpacetonev(T)
% vapor heat capacity (kJ/molK) for acetone with T in deg C
T=T-273.15;
a=71.96e-3;
b=20.10e-5;
c=-12.78e-8;
d=34.76e-12;
cp=a+b*T+c*T.^2+d*T.^3;
return




