%% Post-analysis
analyticSolnVector = [];
MFPTvector = NaN(5,10);
fluxvectorC = NaN(5,10);
endpts = [10,15,20,25, 30];
for j = 1:5
    for k = 1:10
    filename = "Fluxes/Data/fluxOutZ" + endpts(j) + "Run" + k + ".txt";
    if(isfile(filename))
        data=load(filename);
        if paramsWE.fluxBin >= 0
        meanFlux=mean((data)/(paramsDE.dt));
        fluxvectorC(j,k) = meanFlux;
        MFPTvector(j,k) = 1/meanFlux;
        end
    end
    end
        analyticSolnVector = [analyticSolnVector, paramsModel.tauSlow*pi*erfi(endpts(j)/sqrt(2) / paramsModel.sigmax)];
end

MFPTavgs = nanmean(MFPTvector.');

errorVect = nanstd(MFPTvector.');
errorVect(1:4) = errorVect(1:4)/sqrt(10);
errorVect(5) = errorVect(5) / sqrt(5);


figure()
errorbar(10:5:30,MFPTavgs,errorVect,'o')
hold on
plot(10:5:30,analyticSolnVector)
set(gca,'Yscale','log')
title('OU MFPTs')
xlabel('Distance from 0')
ylabel('MFPT')
legend('Measured MFPT','Analytic Solution')
%Add in the equilibrium distribution
% z = linspace(min(xout),max(xout));
% pz =  exp(-z .^2 / (paramsModel.tauSlow * paramsModel.sigmax ^2))/sqrt(paramsModel.sigmax ^2 * pi * paramsModel.tauSlow);
% plot(z,pz,'LineWidth', 2)