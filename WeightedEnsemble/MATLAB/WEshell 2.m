function [newX, newWeights, fluxOut] = WEshell(repPerBin, x, weights,binParams,fluxBin) 
%Basic shell for weighted ensemble functionality
%Assumes that the data is in the current kernel, obviously will need to load data when adapting to C
%Things that will necessarily change in the next iteration
%1. Data location / access / storage
%2. Bin definitions
%x is the variable that will store the data

 %Specify how many bins exactly that is
binIDs = findBin(x,binParams); %Gives the identity of each bin 
nBins = numel(unique(binIDs))-1;

%Measure the flux out
if fluxBin >= 0
    leaveIndices = find(binIDs == fluxBin);
    fluxOut = sum(newWeights(leaveIndices));
    weights(leaveIndices) = 0;
    x(leaveIndices) = NaN;
    weights = weights/sum(weights);
    weights(leaveIndices) = NaN;
end


% move preallocation out of WE loop
newWeights = NaN((nBins+1) * repPerBin, 1); %Preallocate vectors with enough room for copies
newWeights(1:length(weights)) = weights;
newX = NaN((nBins+1) * repPerBin, 1,'double');
newX(1:length(x)) = x;
newIDs = NaN((nBins+1) * repPerBin, 1)
newIDs(1:length(x)) = binIDs;

for n = 1:nBins;
    xInBin = newX(newIDs == n); %Boolean true / false for whether a specific point in the x vector is in the bin
    nInBin = numel(newX(newIDs == n)); %For each bin, find the number of reps in the bin
    
%The result of this allows us to create a reduced vector with only the x components
%in the bin through newX(xInBin).

%In the following comments: "reduced vector" refers to the vector
%newX(xInBin). Extended vector refers to the vector newX.
    
%Below through the following process for each bin. 
%Firstly, checks to make sure that there are replicas in the bin
%Second, while there are more replicas in the bin than necessary, we check
%to see how many replicas inside the bin share the smallest weight. If it's more than one,
%we combine the first two listed in the vector together. Otherwise, we
%combine the smallest with the first found second smallest.
    if nInBin>0 && nInBin > repPerBin 
        while nInBin > repPerBin
            if length(find(newWeights(xInBin)== min(newWeights(xInBin)))) > 1
                mergeNumber = find(newWeights(xInBin) == min(newWeights(xInBin)),2,'first');
                %mergeNumber: Gives the index of the replicas to be combined in the REDUCED vector
                mergeIndex = find(xInBin);
                %Gives the indices inside the extended vector of all of the
                %reduced vector components. E.g. if element 1,5,10 of
                %extended vector are inside the reduced vector, mergeIndex
                %becomes the vector [1,5,10]
                mergeIndex = mergeIndex(mergeNumber);
                %Gives the indices of the replicas to be combined inside
                %the extended vector
            else
                %The following is a similar process, but specifies a
                %"second smallest" weight for when there's only one
                %smallest
                mergeNum1 = find(newWeights(xInBin) == min(newWeights(xInBin)),1,'first');
                mergeNum2 = find(newWeights(xInBin) == min(setdiff(newWeights(xInBin),min(newWeights(xInBin)))),1,'first');
                mergeIndex = find(xInBin);
                mergeIndex = mergeIndex([mergeNum1,mergeNum2]);
            end
            %Now that we know the locations of the replicas to be combined
            %in the overall vector, we combine their weights and replace
            %the extras with NaNs.
            newWeights(mergeIndex(1)) = sum(newWeights(mergeIndex));
            newWeights(mergeIndex(2)) = NaN;
            newX(mergeIndex(1)) = mean(newX(mergeIndex)); %This can be commented to cause one to absorb the other, currently combines them
            newX(mergeIndex(2)) = NaN;
            %Following 2 lines updates the counts for number of entries in
            %the bin and the boolean true false that helps us locate
            %replicas
            binIDs = findBin(x,binParams)
            xInBin = newX(newIDs == n);
            nInBin = numel(newX(newIDs == n));
        end
%This following section does something similar, except it creates copies of
%the larger replicas. splitNumber gives the index of the splittee in the
%reduced vector, splitIndex becomes the location in the extended vector,
%then we place the copy in the first NaN contained in the vector.
    else if nInBin > 0 && nInBin < repPerBin
        while nInBin < repPerBin
            splitNumber = find(newWeights(xInBin) == max(newWeights(xInBin)),1,'first');
            splitIndex = find(xInBin);
            splitIndex = splitIndex(splitNumber);
            newWeights(splitIndex(1)) = newWeights(splitIndex(1))/2;
            newWeights(find(isnan(newWeights),1,'first')) = newWeights(splitIndex);
            newX(find(isnan(newX),1,'first')) = newX(splitIndex);
            binIDs = findBin(x,binParams)
            xInBin = newX(newIDs == n);
            nInBin = numel(newX(newIDs == n));
        end
    end
    end
end

if fluxTF
   
    
    
end


%Remove all NaN entries from the output
    filledIndex = ~isnan(newX);
    newX = newX(filledIndex);
    newWeights = newWeights(filledIndex);
end
        