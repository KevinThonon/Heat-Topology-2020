
% inputs: N, pctmetal, 0.4, dcda_mat
function [xnew] = OC(N, x, volfrac, dc)
%OC Summary of this function goes here
%   Detailed explanation goes here
l1 = 0;
l2 = 100000;
move = 0.2;
while (l2-l1)>1e-4
    lmid = 0.5*(l2+l1);
    xnew = max(0.001, max(x-move,min(1.,min(x+move,x.*sqrt(abs(-dc./lmid))))));
    if sum(sum(xnew))-volfrac*N*N>0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end

