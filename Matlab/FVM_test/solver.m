function [solmat,sol] = solver(dp,LL,RL)
%SOLVER Summary of this function goes here
%   Detailed explanation goes here

sol = sparse(LL)\RL;

solmat = zeros(dp,dp);

for i = 1:dp
    solmat(:,i)=sol(i*dp-dp+1:i*dp);
end

end

