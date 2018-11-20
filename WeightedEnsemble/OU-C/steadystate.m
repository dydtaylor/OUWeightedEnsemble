%Create a weighted histogram with associated bins
histostep = .5;
histobins = -12:histostep:12;
sigmaX = 3.1385;

wHBin = [];

data = load('SteadyState/simOut7.txt');

%Sum the total WE weights contained in each histogram bin
for n = histobins
    xInBin = data(:,1) > n & data(:,1) < n+histostep;
    wInBin = sum(data((xInBin),2));
    wHBin = [wHBin, wInBin];
end
figure();
%Make the histogram. Add histostep/2 to keep graph centered, scale by
%histostep to keep area of graph = 1
bar(histobins+histostep/2, wHBin/histostep);
hold on
%Add in the equilibrium distribution
z = linspace(min(data(:,1)),max(data(:,1)));
pz = exp(-z.^2 / (2*sigmaX^2)) /sqrt(2*pi*sigmaX^2);
plot(z,pz,'LineWidth', 2)