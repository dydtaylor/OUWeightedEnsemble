runningRead = zeros(9000,1);
runningRob = runningRead;

for j = 1:9000
    runningRead(j) = mean(fluxes(j:j+1000));
    runningRob(j) = mean(fluxAtTauStep(j:j+1000));
end
% figure()
% plot(runningRead)
% hold on
% plot(runningRob)
% set(gca,'YScale','log')
% hline = refline(0,(1/analyticSoln) * paramsWE.tau * paramsDE.dt);
% hline.LineWidth = 1.5;
% hline.Color = 'g';
% legend('Read','Rob','Analytic')
% robline = refline(0,meanFlux);
% robline.LineWidth = 1.5;
% robline.Color = 'r';
% readline = refline(0,

[robf1, robx1] = ecdf(fluxAtTauStep(3333:6666));
[robf2, robx2] = ecdf(fluxAtTauStep(6666:end));
[readF1, readX1] = ecdf(fluxes(3333:6666));
[readF2, readX2] = ecdf(fluxes(6666:end));
figure()
semilogx(readX1,readF1)
hold on
semilogx(readX2,readF2)
figure()
semilogx(robx1,robf1)
hold on
semilogx(robx2,robf2)

figure()
semilogx(robx2,robf2)
hold on
semilogx(readX2,readF2)

figure()
semilogy(fluxes)
hold on
% semilogy(fluxAtTauStep)
hline = refline(0,(1/analyticSoln) * paramsWE.tau * paramsDE.dt);
hline.LineWidth = 1.5;
hline.Color = 'g';