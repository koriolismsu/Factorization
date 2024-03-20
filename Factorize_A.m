clear all
%%%% let us try to factorize matrix A = Q_plus Q_minus 

%% Introduce a loop
k = 1+0.01i;
step1 = 0.01;
step2 = 0.0001;
xi_plus_loop = [-0.5:step1:1.5,1.5 + 1i*(0:step2:0.1),1.5 + 0.1i - (0:step1:2),-0.5+ 1i*(0.1:-step2:0)];

xi_minus_loop = xi_plus_loop - 0.1i - 1;

xi_whole_loop = [-1.5:step1:1.5,1.5 + 1i*(0:step2:0.1),1.5 + 0.1i - (0:step1:3), -1.5 + 1i*(0.1:-step2:0)] - 0.04i;

xi_cur_loop = xi_minus_loop;

gamma_0 = sqrt(k^2  -  xi_cur_loop(1).^2);
root3_1a = (xi_cur_loop(1) + 1i*gamma_0)^(2/3);
root3_1b = (xi_cur_loop(1) - 1i*gamma_0)^(2/3);
root3_2a = (xi_cur_loop(1) + 1i*gamma_0)^(1/3);
root3_2b = (xi_cur_loop(1) - 1i*gamma_0)^(1/3);

root3_1 = root3_2a+root3_2b;
root3_2 = root3_2a+root3_2b;
% 
 figure;
 plot(xi_cur_loop)
 hold all
 plot(k,'*')
 hold all
 plot(-k,'*')

%%%%  evaluate roots we need
square_cont = root_cont(1/2,gamma_0,k^2 - xi_cur_loop.^2); 
%%%% whats below lives on the triple covering of a tube
covering_function_1a = root_cont(2/3,root3_1a, xi_cur_loop + 1i*square_cont);
covering_function_2a = root_cont(1/3,root3_2a, xi_cur_loop + 1i*square_cont);

covering_function_1b = root_cont(2/3,root3_1b, xi_cur_loop - 1i*square_cont);
covering_function_2b = root_cont(1/3,root3_2b, xi_cur_loop - 1i*square_cont);


figure;
plot(covering_function_1a,'*')

% covering_function_1a(end)/root3_1b


%%%% Build Q_plus matrix, factor of A analytic in upper half_plane
%% Check that it satisfies bypass problem

M_mat = @(x) [0,1,1;-sqrt(3)/2/x,-1/2,1;sqrt(3)/2/x,-1/2,1].';
Eig_mat = @(x1,x2) diag([1,x1,x2]);

H1_plus_start = M_mat(gamma_0)*Eig_mat(root3_2a,root3_1a)*inv(M_mat(gamma_0));


H1_plus_end = M_mat(square_cont(end))*Eig_mat(covering_function_2a(end),covering_function_1a(end))*inv(M_mat(square_cont(end)));


H2_plus_start = M_mat(-gamma_0)*Eig_mat(root3_2b,root3_1b)*inv(M_mat(-gamma_0));


H2_plus_end = M_mat(-square_cont(end))*Eig_mat(covering_function_2b(end),covering_function_1b(end))*inv(M_mat(-square_cont(end)));

%bypass Matrix % P_matrix(gamma_0) = A_mat(square_cont(end))*inv(A_mat(gamma_0))

P_matrix = @(x) [-1/2,1/2i/x,-1/2i/x;-1i*x/2,1/2,1/2;1i*x,1,0];

%%%%% check that H_plus = P_matrix*H_plus after bypass whole_loop


H2_plus_end - P_matrix(gamma_0)*H2_plus_start


%%%% build Qplus 

Qplus_end = H1_plus_end + H2_plus_end;
Qplus_start = H1_plus_start + H2_plus_start;

Qplus_end- Qplus_start 

Qplus_end - P_matrix(gamma_0)*Qplus_start


%%% check that Qminus does not have branching in lower half-plane

A_mat = @(x) [-1/2i/x,1/2i/x,1/2;1/2,1/2,1i/2*x; 0,1,-1i*x];

Qminus_start = inv(Qplus_start)*A_mat(gamma_0);
Qminus_end =  inv(Qplus_end)*A_mat(square_cont(end));

Qminus_end - Qminus_start % this is A = QplusQminus factorization


%%% factorization is not commutative

Qminus_start*Qplus_start - Qplus_start*Qminus_start


%% let us build A = QminusQplus factorization

%%%% P_matrix(gamma_0) - inv(A_mat(gamma_0))*A_mat(square_cont(end))

P_matrix = @(x)[0, 1, -1i*x; 1/2, 1/2, (1i*x/2); (1/(2i*x)), -1/(2i*x), -(1/2)];

M_mat = @(x) [1,1,0;2/sqrt(3)*x,-1/sqrt(3)*x,1;-2/sqrt(3)*x,1/sqrt(3)*x,1].';


H1_plus_start = M_mat(gamma_0)*Eig_mat(root3_2a,root3_1a)*inv(M_mat(gamma_0));


H1_plus_end = M_mat(square_cont(end))*Eig_mat(covering_function_2a(end),covering_function_1a(end))*inv(M_mat(square_cont(end)));


H2_plus_start = M_mat(-gamma_0)*Eig_mat(root3_2b,root3_1b)*inv(M_mat(-gamma_0));


H2_plus_end = M_mat(-square_cont(end))*Eig_mat(covering_function_2b(end),covering_function_1b(end))*inv(M_mat(-square_cont(end)));

%%%% build Qplus 

Qplus_end = H1_plus_end + H2_plus_end;
Qplus_start = H1_plus_start + H2_plus_start;

Qplus_end - Qplus_start 

Qplus_end - Qplus_start*P_matrix(gamma_0)

Qminus_start = A_mat(gamma_0)*inv(Qplus_start);
Qminus_end =  A_mat(square_cont(end))*inv(Qplus_end);

Qminus_end - Qminus_start % this is A = QminusQplus factorization
