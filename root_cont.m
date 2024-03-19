function result = root_cont(mu,xi_start,xi)
         Nit = 20;
         N = length(xi);
         % define the first value like this and then continue along the boundary
         xi_prev =  xi_start;  
         result = zeros(1,N);
         for n = 1:N
            xi_prev = mynthroot(xi_prev,xi(n),mu,Nit);
            result(n) = xi_prev; 
         end
         
end




function outvar = mynthroot(x_cur,val,mu,Nit)
            for n_int=1:Nit
                x_new = x_cur - (x_cur^(1/mu) - val)/((1/mu)*x_cur^(1/mu-1));
                x_cur = x_new;
            end
            outvar = x_cur;
end