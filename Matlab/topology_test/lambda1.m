function [lambda] = lambda1(T, K, dp)
lambda = transpose(K)\((-2/dp^2)*T);
end

