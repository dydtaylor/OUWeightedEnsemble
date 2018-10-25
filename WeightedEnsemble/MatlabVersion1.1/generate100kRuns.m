paramsWE.tau = 10; %number of time steps to run OU before WE
paramsWE.repsPerBin = 32; %Amount of replicas for each bin in WE
paramsWE.tauMax = 100000;
compTime = zeros(10,5);
for Zstar = [10,15,20,25,30]
    paramsWE.nBins = 2+10*Zstar;
    paramsWE.fluxBin = 2+10*Zstar;
    paramsWE.binDefs = cell(paramsWE.nBins,2);
    paramsWE.binDefs{1,1} = [-Inf,0];
    nc = 2;
    for j = 0:0.1:Zstar
        paramsWE.binDefs{nc,1} = [j,j+.1;];
        nc = nc+1;
    end
    paramsWE.binDefs{paramsWE.nBins,1} = [Zstar, Inf];
    for runNo = 1:10
        filename = Zstar + "." + runNo;
        tic
        WEscrap(filename,paramsWE)
        compTime(runNo,Zstar) = toc;
    end
end