%% Post-analysis

load('WERun.mat');

%Create a weighted histogram with associated bins
histostep = .5;
histobins = -6:histostep:6;

wHBin = [];

%Sum the total WE weights contained in each histogram bin
for n = histobins
    xInBin = xout > n & xout < n+histostep;
    wInBin = sum(weights(xInBin));
    wHBin = [wHBin, wInBin];
end
figure();
%Make the histogram. Add histostep/2 to keep graph centered, scale by
%histostep to keep area of graph = 1
bar(histobins+histostep/2, wHBin/histostep);
hold on
%Add in the equilibrium distribution
z = linspace(min(xout),max(xout));
pz = exp(-z.^2 / 2) /sqrt(2*pi);
plot(z,pz,'LineWidth', 2)