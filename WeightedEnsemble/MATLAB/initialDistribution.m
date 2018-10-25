function xo = initialDistribution(paramsModel,n)

    u = rand(n,1);
    sigma = paramsModel.sigmax;
    xo = sigma * sqrt(2) * erfinv(2.*u - 1);
    
end