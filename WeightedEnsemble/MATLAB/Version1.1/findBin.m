function iBin = findBin(x,paramsBinning)
%Question of how we define the parameters. We can set up boundaries or just
%have each component fully define the compartment. I'm going to write it
%such that each element of paramsBinning fully defines one bin.

%For a 1-d parameter space, we'll assume that paramsBinning has the
%following format:
%Cell with size nBins x 2
%First cell: Array with lower and upper bounds of intervals
%Second cell: Array with individual points

iBin = NaN(length(x),1);

for j = 1:length(paramsBinning)
    xInBin = zeros(length(x),1);
    lowerBounds = paramsBinning{j,1}(1:2:end);
    upperBounds = paramsBinning{j,1}(2:2:end);
    points = paramsBinning{j,2};

    for i = 1:length(lowerBounds)
        xInInt = x>lowerBounds(i) & x < upperBounds(i);
        xInBin = xInBin + xInInt;
    end
    for i = 1:length(points)
        xAtPoint = x == points(i);
        xInBin = xInBin + xAtPoint;
    end
    xInBin = logical(xInBin);
    iBin(xInBin) = j;
end

iBin(isnan(iBin)) = 0;
end