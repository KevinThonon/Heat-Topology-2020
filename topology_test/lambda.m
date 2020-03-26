function [l] = lambda(T,K, N)
    ne = N*N;
    
    dgdu = 2*transpose(T)/ne;
    
    l = transpose(K)\(-transpose(dgdu));
end

