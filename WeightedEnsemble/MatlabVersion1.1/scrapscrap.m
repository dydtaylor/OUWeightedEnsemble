figure()
hold on
for j = 1:10
filename = "WERun"+j+".mat";
load(filename)
[f,x] = ecdf(fluxAtTauStep(3*end/4 : end));
semilogx(x,f)
end
set(gca,'Xscale','log')
title('Rob OU Last Quarter CDFs')
figure()
hold on
for j = 1:10
    filename = "ou_z25_c0_morebins/ou_tau1E-8_" + j + "/output/WE_flux.mat";
    load(filename)
    [f,x] = ecdf(flux_est(75000:100000));
    semilogx(x,f)
end
set(gca,'Xscale','log')
title('Brian OU Last Quarter CDFs')