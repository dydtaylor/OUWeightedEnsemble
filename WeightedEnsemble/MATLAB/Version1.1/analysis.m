%% Post-analysis
analyticSolnVector = [];
MFPTvector = zeros(5,10);
fluxVector = NaN(5,10);
endpts = [10,15,20,25];
for j = 1:4
    for k = 1:10
    filename = "WERunZ" + endpts(j) + "." + k + ".mat";
    load(filename)
        if paramsWE.fluxBin >= 0
        meanFlux=mean(fluxAtTauStep(end/4:end)/(paramsWE.tau * paramsDE.dt));
        fluxVector(j,k) = meanFlux;
        MFPTvector(j,k) = 1/meanFlux;
        end
    end
        analyticSolnVector = [analyticSolnVector, paramsModel.tauSlow*pi*erfi(paramsWE.binDefs{paramsWE.fluxBin,1}(1)/sqrt(2) / paramsModel.sigmax)];
end

for k = 1:5
    filename = "WERunZ30."+k+".mat";
    load(filename)
        if paramsWE.fluxBin >= 0
        meanFlux=mean(fluxAtTauStep(end/4:end)/(paramsWE.tau * paramsDE.dt));
        fluxVector(5,k) = meanFlux;
        MFPTvector(5,k) = 1/meanFlux;
        end
end

analyticSolnVector = [analyticSolnVector,paramsModel.tauSlow*pi*erfi(paramsWE.binDefs{paramsWE.fluxBin,1}(1)/sqrt(2) / paramsModel.sigmax)];

MFPTvector(5,6:10) = NaN;
MFPTavgs = nanmean(MFPTvector.');

errorVect = nanstd(MFPTvector.');
errorVect(1:4) = errorVect(1:4)/sqrt(10);
errorVect(5) = errorVect(5) / sqrt(5);


figure()
errorbar(10:5:30,MFPTavgs,errorVect,'o')
hold on
plot(10:5:30,analyticSolnVector)
set(gca,'Yscale','log')
%Add in the equilibrium distribution
% z = linspace(min(xout),max(xout));
% pz =  exp(-z .^2 / (paramsModel.tauSlow * paramsModel.sigmax ^2))/sqrt(paramsModel.sigmax ^2 * pi * paramsModel.tauSlow);
% plot(z,pz,'LineWidth', 2)