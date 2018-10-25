x0 = zeros(100000,1);
xout = oneDOU(x0,.005,100000);
figure()
histogram(xout,50,'Normalization', 'pdf');
hold on
z = linspace(min(xout),max(xout));
pz = exp(-z.^2 /2) / sqrt(2*pi);
plot(z,pz,'LineWidth',2)