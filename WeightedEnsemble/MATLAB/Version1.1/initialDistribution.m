function [x0, weights] = initialDistribution(paramsModel,n)

%     u = rand(n,1);
%     sigma = paramsModel.sigmax;
%     xo = sigma * sqrt(2) * erfinv(2.*u - 1);
    x0 = zeros(1,n);
    weights = ones(numel(x0),1)./numel(x0);
end