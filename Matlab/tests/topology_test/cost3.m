function [cost] = cost3(T,dp)
A = (1/dp)^2;
w = A*ones(dp^2,1);
w(0) = A/4;
w(N) = A/4;
w(dp^2) = A/4;
w(dp^2-dp+1) = A/4;

for i = 2:1:(dp-1)
    w(i) = A/2;
    w(dp^2-dp + i) = A/2;
    w(i*dp) = A/2;
    w(i*dp-dp+1) = A/2;
end

cost = dot(w,T);

end

