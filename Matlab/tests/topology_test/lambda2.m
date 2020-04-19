function [lambda] = lambda2(K, dp)
w = -1*ones(dp^2,1);
lambda = transpose(K)\w;
end

