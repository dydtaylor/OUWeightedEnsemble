%First, create some empty vectors then run the 1d OU and WE processes

% --- put dynamics steps into for-loops


% WE params
paramsWE.tau = 10; %number of time steps to run OU before WE
paramsWE.repsPerBin = 50; %Amount of replicas for each bin in WE
paramsWE.tauMax = 100;

% Dynamics engine parameters
paramsDE.dt = .005; %time step for OU
% Model parameters
paramsModel.tauSlow = 1;
paramsModel.sigmax = 1;

% for order parameter definition and binning
paramsBinning=cell{nBins,2};
for j = 1:nBins
    paramsBinning{1,j}
end    

% -- put parameters into struct

% --- separate m-file for generating inital xout and weights.
% initial distribution and weights of trajectories
xout = initialDistribution(paramsModel,1000); %data containing the replicas, runs 
weights = ones(numel(xout),1)./numel(xout); %data containing the weights of each replica


for nWE = 1:paramsWE.tauMax
    [xout, weights] = WEshell(repsPerBin,xout,weights,paramsBinning);
    xout = oneDOU(xout,tau,paramsDE,paramsModel);
end
fprintf('The sum of all the weights is %f \n',sum(weights))


save('WERun.mat', 'xout', 'weights');