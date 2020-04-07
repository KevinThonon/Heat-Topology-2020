function [mat] = x4func(N,h)
dp = N+1;
k = linspace(0,0.01,dp);
K = repmat(k,dp,1);
Solmat = -12*(K.^2);
for i = 2:(dp-1)
    for j = 2:(dp-1)
        Solmat(i,j) = -Solmat(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    Solmat(1,i) = -Solmat(1,i)*h^2/2;
    Solmat(dp,i) = -Solmat(dp,i)*h^2/2;
end

Solmat(1:dp,1) = zeros(1,dp);
Solmat(1:dp,dp) = ones(1,dp)*0.01^4;

mat = Solmat(:);
end

