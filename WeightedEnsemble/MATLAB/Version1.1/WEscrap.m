%First, create some empty vectors then run the 1d OU and WE processes

% --- put dynamics steps into for-loops
function WEscrap(zstar,paramsWE)
% WE params
% paramsWE.tau = 10; %number of time steps to run OU before WE
% paramsWE.repsPerBin = 32; %Amount of replicas for each bin in WE
% paramsWE.tauMax = 10000;
% paramsWE.nBins = 302;
% paramsWE.binDefs = cell(paramsWE.nBins,2);
% paramsWE.fluxBin = 302;

% Dynamics engine parameters
paramsDE.dt = 1e-9; %time step for OU
% Model parameters
paramsModel.tauSlow = 1.05 * 10^-6;
paramsModel.sigmax = 3.1385;

% --- separate m-file for generating inital xout and weights.
% initial distribution and weights of trajectories
[xout, weights] = initialDistribution(paramsModel,10); %data containing the replicas, runs 
wtSum = zeros(paramsWE.tauMax,1);
fluxAtTauStep = zeros(ceil(paramsWE.tauMax/2),1);

for nWE = 1:250
    fprintf('%f \n', nWE)
    [xout, weights] = WESplitMerge(xout,weights,paramsWE);
    xout = oneDOU(xout,paramsWE,paramsDE,paramsModel);
    if abs(sum(weights) - 1 )> 10^-4
        error('Weights no longer normalized')
    end
end

for nWE = 250: paramsWE.tauMax
    fprintf('%f \n',nWE)
    [xout, weights, fluxAtTauStep(nWE)] = WEcomputeFluxes(xout,weights,paramsWE);
    [xout, weights] = WESplitMerge(xout,weights,paramsWE);
    xout = oneDOU(xout,paramsWE,paramsDE,paramsModel);
    if abs(sum(weights) - 1 )> 10^-4
        error('Weights no longer normalized')
    end
end

fprintf('The sum of all the weights is %f \n',sum(weights))
if paramsWE.fluxBin >= 0
meanFlux=mean(fluxAtTauStep(end/4:end)/(paramsWE.tau * paramsDE.dt));
MFPT = 1/meanFlux;
analyticSoln = paramsModel.tauSlow*pi*erfi(paramsWE.binDefs{paramsWE.fluxBin,1}(1)/sqrt(2) / paramsModel.sigmax);
err = MFPT/analyticSoln;
end
save("WERunZ" + zstar +".mat", 'xout', 'weights','fluxAtTauStep','paramsWE','paramsModel','paramsDE');
end