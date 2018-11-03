function [repOut, weightOut, fluxOut] = WEcomputeFluxes(repIn, weightIn, paramsWE)

fluxBin = paramsWE.fluxBin;
binParams = paramsWE.binDefs;
binIDs = findBin(repIn,binParams); %Gives the identity of each bin 
repOut = repIn;
weightOut = weightIn;

%Measure the flux out
if fluxBin >= 0
    leaveIndices = find(binIDs == fluxBin);
    fluxOut = sum(weightOut(leaveIndices));
    weightOut(leaveIndices) = 0;
    repOut(leaveIndices) = NaN;
    weightOut = weightOut/sum(weightOut);
    weightOut(leaveIndices) = NaN;
else
    fluxOut = NaN;
end

filledIndex = ~isnan(repOut);
repOut = repOut(filledIndex);
weightOut = weightOut(filledIndex);

end