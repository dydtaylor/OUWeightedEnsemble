load('fluxes.mat')
load('WERun.mat')

[fRead, xRead] = ecdf(fluxes(end/4:end));
[fRob, xRob] = ecdf(fluxAtTauStep(end/4:end));
figure()
semilogx(xRob,fRob);
hold on
semilogy(xRead,fRead);

meanFlux=mean(fluxAtTauStep(end/4:end)/(paramsWE.tau * paramsDE.dt));
robMFPT = 1/meanFlux
estFlux = mean(fluxes(end/4:end)/(5E-9));
readMFPT = 1/estFlux
analyticSoln = paramsModel.tauSlow*pi*erfi(paramsWE.binDefs{paramsWE.fluxBin,1}(1)/sqrt(2) / paramsModel.sigmax)

figure()
semilogy(fluxAtTauStep)
hold on
semilogy(fluxes)