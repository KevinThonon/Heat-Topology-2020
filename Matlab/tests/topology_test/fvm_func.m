function [solution, LL, RL] = fvm_func(kmat, N, q)
    %q = 2/(0.01*0.01*0.001);
    h = 0.01/N;
    dp = N+1;

        
    [LL,RL] = FVM_mat(q,h,dp,kmat);
    solution = LL\RL;
    
    
end

