%% Simulates a 2D Ornstein-Ulenbeck process. This is equivalent to two independent OU processes x(t) and y(t) that are linearly combined to z(t).

% Parameters
tauFast =1; % s
tauSlow = 1; % s

sigmax = 1; %nm
sigmay = 1; %nm

%% time loop

% timestep and max run time
dt = 0.005*tauFast;
ntMax = 1000*ceil(5e0*tauSlow/dt);

% variable for data storage
xyArray = zeros(ntMax,2);

% initial condition
xy = [0;0];

% Euler-Maruyama
for nt=1:ntMax
    
    xy = xy - 1./[tauSlow; tauFast].*xy*dt + ([sigmax; sigmay]./(sqrt([tauSlow; tauFast])))*sqrt(2*dt).*randn(2,1);
    
    xyArray(nt,:) = xy';
    
end % finished time loop


%% plot and analyze

% variables to specify autocorrelation function's "horizontal axis" (i.e.,
% which deltat's to find alpha at). 

ideltaTauMax  = floor(1e-5/dt);
ideltaTauSkip = 5e1;

ideltaTauArray = 1:ideltaTauSkip:ideltaTauMax;
deltaTauArray = ideltaTauArray*dt;

%% -------- process x -------- 
figure(21); clf; 

% plot time series
subplot(3,1,1); hold on; box on;
title('Process x(t) - the slow one');
variableToPlot = xyArray(1:100:end,1);
plot(dt*(1:(numel(variableToPlot))),variableToPlot);

% plot stationary distribution
subplot(3,1,2); hold on; box on;
hist(xyArray(:,1));

std(xyArray(:,1))


% %% autocorrelation
% alpha_x   = autocorr(xyArray(:,1),  ideltaTauArray);
%%

subplot(3,1,3); hold on; box on;
plot(deltaTauArray,alpha_x);
plot(deltaTauArray,0*deltaTauArray,'-k');

%% -------- process y --------
figure(22); clf; 

% plot time series
subplot(3,1,1); hold on; box on;
title('Process y(t) - the fast one');
variableToPlot = xyArray(1:100:end,2);
plot(dt*(1:(numel(variableToPlot))),variableToPlot);

% plot stationary distribution
subplot(3,1,2); hold on; box on;
hist(xyArray(:,2));

std(xyArray(:,2))

% %% autocorrelation
% alpha_y   = autocorr(xyArray(:,2),  ideltaTauArray);

%%
subplot(3,1,3); hold on; box on;
plot(deltaTauArray,alpha_y);
plot(deltaTauArray,0*deltaTauArray,'-k');

%% -------- transformed process --------

c = 0.5; % mixture variable

z = c*xyArray(:,1) + (1-c)*xyArray(:,2);

figure(23); clf; 

% plot time series
subplot(3,1,1); hold on; box on;
title('Transformed process z(t) - the composite one')
variableToPlot = z(1:100:end);
plot(dt*(1:(numel(variableToPlot))),variableToPlot);

% plot stationary distribution
subplot(3,1,2); hold on; box on;
hist(z);

%% autocorrelation
alpha_z   = autocorrelation(z,  ideltaTauArray);

%% 
subplot(3,1,3); hold on; box on;
plot(deltaTauArray,alpha_z, '-b');

% two-component autocorrelation ansatz

c2 = 0.98;

alpha_twoComponent = (1-c2)*exp(-deltaTauArray/tauFast) + c2*exp(-deltaTauArray/tauSlow);

plot(deltaTauArray,alpha_twoComponent, '-r');
%plot(deltaTauArray,0*deltaTauArray,'-k');


legend('From OU simulation', 'Ansatz')

plot(deltaTauArray,0*deltaTauArray,'-k');

