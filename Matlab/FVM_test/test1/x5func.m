function [mat,der2] = x4func(dp,h)
k = linspace(0,0.01,dp);
K = repmat(k,dp,1);
%gewone oplossing

mat = K.^5;
mat = mat(:);

der2 = -20*K.^3;
for i = 2:(dp-1)
    for j = 2:(dp-1)
        der2(i,j) = -der2(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    der2(1,i) = -der2(1,i)*h^2/2;
    der2(dp,i) = -der2(dp,i)*h^2/2;
end

der2(1:dp,1) = zeros(1,dp);
der2(1:dp,dp) = ones(1,dp)*0.01^5;
der2 = der2(:);
end
