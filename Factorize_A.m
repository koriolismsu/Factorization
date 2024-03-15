%%%% let us try to factorize matrix A

%% Introduce a loop
k = 1+0.01i;
step = 0.01;
xi_plus_loop = [0:step:1.5,1.5 + 1i*(0:step:0.1),1.5 + 0.1i - (0:step:1.5), 1i*(0.1:-step:0)];

gamma_0 = sqrt(k);

figure;
plot(xi_plus_loop,'*')