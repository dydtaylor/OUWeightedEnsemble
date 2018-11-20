figure()
hold on
endpts = [10,15,20,25,30];
    for k = 1:10
       
        filename = "Fluxes/Data/fluxOutZ" + endpts(3) + "Run"+k + ".txt";
    if isfile(filename)
        data=load(filename);
        [f,x] = ecdf(data(floor(end/3) : end));
        semilogx(x,f)
    end
    end
set(gca,'Xscale','log')
title('Rob OU Last Quarter CDFs')