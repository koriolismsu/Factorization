%% Introduce a loop
k = 1+0.01i;
step1 = 0.01;
step2 = 0.0001;
xi_plus_loop = [0:step1:1.5,1.5 + 1i*(0:step2:0.1),1.5 + 0.1i - (0:step1:1.5), 1i*(0.1:-step2:0)];

xi_minus_loop =[0:-step1:-1.5,-1.5 - 1i*(0:step2:0.1), -1.5-0.1i + (0:step1:1.5), -1i*(0.1:-step2:0)];

xi_whole_loop = [-1.5:step1:1.5,1.5 + 1i*(0:step2:0.1),1.5 + 0.1i - (0:step1:3), -1.5 + 1i*(0.1:-step2:0)] - 0.04i;

xi_cur_loop = [xi_minus_loop];

gamma_0 = sqrt(k^2  -  xi_cur_loop(1).^2);
root3_1 = (xi_cur_loop(1) + 1i*gamma_0)^(1/3);
root3_2 = (xi_cur_loop(1) - 1i*gamma_0)^(1/3);

% 
 figure;
 plot(xi_cur_loop,'*')
 hold all
 plot(k,'*')
 hold all
 plot(-k,'*')
  hold all
 plot(0,'*')

%%%%  evaluate roots we need
square_cont = root_cont(1/2,gamma_0,k^2 - xi_cur_loop.^2); 
%%%% whats below lives on the triple covering of a tube
covering_function = root_cont(1/3,root3_1, xi_cur_loop + 1i*square_cont)+...
    root_cont(1/3,root3_2, xi_cur_loop - 1i*square_cont);

figure;
plot(covering_function,'*')


covering_function(end)/covering_function(1)

