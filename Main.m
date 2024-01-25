clear
close all
% Input parameters 
yMF = 0.28; zMF = 0.40; % Mean Field concentrations
SigmaFluct = [0.2 -0.11
    -0.11 0.24]; % Covariance of 
N = 20; % Number of motors
f0 = 3; % Force of individual motor

%% Evaluate integral numerically
% Grid over q space
q1 = 0.1: 0.01:15; q2 = 0.1: 0.01:15;
% Grid over the f-axis
f_grid = -5:0.5:30;
% Exact integral 
[f_grid,pdf] = Exact_Integral(q1,q2,yMF,zMF,SigmaFluct,N,f0,f_grid);
% Approximate formulas to different orders
[gForceMF,Iformula2,Iformula4] = Approximate_Formula(yMF,zMF,SigmaFluct,N,f0); 


%%
figure(1)
plot(f_grid,pdf,'DisplayName','Full integral','LineWidth',1.2)
hold on 
plot(f_grid,gForceMF(f_grid),'--','DisplayName','Mean Field','LineWidth',0.5)
plot(f_grid,Iformula2(f_grid),'--','DisplayName','2nd order','LineWidth',0.5)
plot(f_grid,Iformula4(f_grid),'-','DisplayName','4th order','LineWidth',1.2)
xlabel("$f$",'Interpreter','latex','FontSize',20)
ylabel("pdf",'Interpreter','latex','FontSize',20)
legend()







