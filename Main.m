clear
close all
%%
% Input parameters for fit
k1= 70; km1= 10; k2= 100; km2= 20; k3= 50; % Chemical rates
N = 20; % Number of motors
f0 = 3; % Force of individual motors

% The two lines below are the only thing for the fit
k = [k1,km1,k2,km2,k3]; % Stack chemical rates for input function
Iformula4 = Analytical_Curve_Fit(k,f0,N);
% Note that Iformula4 is a function, we do not need to specify a grid. We
% can then evaluate it on the grid given by the data.

timeit(@() Analytical_Curve_Fit(k,f0,N)) % Check how long it takes: for me it was ~10^(-5)

%% This is just to compare with the exact integral
[yMF,zMF,SigmaFluct] = get_MeanField_Fluct(k);
f_grid = -5:0.5:30; f_grid = f_grid.';
% Exact integral 
q1 = 0.1: 0.01:15; q2 = 0.1: 0.01:15;
pdf = Exact_Integral(q1,q2,yMF,zMF,SigmaFluct,N,f0,f_grid);

figure(1)
plot(f_grid,pdf,'DisplayName','Full integral','LineWidth',1.2)
hold on 
plot(f_grid,Iformula4(f_grid),'-','DisplayName','4th order','LineWidth',1.2)
xlabel("$f$",'Interpreter','latex','FontSize',20)
ylabel("pdf",'Interpreter','latex','FontSize',20)
legend()