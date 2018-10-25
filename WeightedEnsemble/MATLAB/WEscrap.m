%First, create some empty vectors then run the 1d OU and WE processes

x0 = zeros(100,1);
dt = .005; %time step for OU
tau = 10; %number of time steps to run OU before WE
paramsModel.tauSlow = 1;
paramsModel.sigmax = 1;
xout = oneDOU(x0,dt,tau); %data containing the replicas, runs 
weights = ones(numel(xout),1)./numel(xout); %data containing the weights of each replica
binWidth = .5; %binWidth for WE step
repsPerBin = 100; %Amount of replicas for each bin in WE

for nWE = 1:100
    [xout, weights] = WEshell(repsPerBin,xout,weights,binWidth);
    xout = oneDOU(xout,dt,tau, paramsModel);
end
fprintf('The sum of all the weights is %f \n',sum(weights))

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