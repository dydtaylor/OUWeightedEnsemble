function xout = oneDOU(x0, paramsWE,paramsDE,paramsModel)

%% Simulates a 1D Ornstein-Ulenbeck process. 
% Parameters
tauSlow = paramsModel.tauSlow; % s

sigmax = paramsModel.sigmax; %nm
ntMax = paramsWE.tau;
dt = paramsDE.dt;

%% time loop

% initial condition
x = x0;

% Euler-Maruyama
for nt=1:ntMax
    
    x = x - 1/(tauSlow)*x*dt + (sigmax/(sqrt(tauSlow)))*sqrt(2*dt).*randn(size(x));
    
end % finished time loop

xout = x;
end