function [c] = costfunc(temp, dp)

    c = transpose(temp)*temp/(dp*dp);
end

