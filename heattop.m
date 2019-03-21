function [u] = heattop(N,pctmetal,pen)
K = zeros(N^2);
h = 0.01/N; %?
f = ones(N^2,1)*2*h^2 *10^7; %lengte
%pctmetal(1:N,1:N) = 0;
k = (65-0.2)*pctmetal.^pen + 0.2;

Dir = [round(.3*N):round(.7*N)+1 ...
    pos(1,N)-1+(round(.3*N):round(.7*N)+1)];

for i = 1:N
    for j = 1:N
        if ~ismember(pos(i,j), Dir)
            xin=0; xis=0; xiw=0; xie=0;
            if i > 1
                xin = meank(k(pos(i,j)),k(pos(i-1,j)));
                K(pos(i,j),pos(i-1,j)) = -xin;
            end
            if i < N
                xis = meank(k(pos(i,j)),k(pos(i+1,j)));
                K(pos(i,j),pos(i+1,j)) = -xis;
            end
            if j > 1
                xiw = meank(k(pos(i,j)),k(pos(i,j-1)));
                K(pos(i,j),pos(i,j-1)) = -xiw;
            end
            if j < N
                xie = meank(k(pos(i,j)),k(pos(i,j+1)));
                K(pos(i,j),pos(i,j+1)) = -xie;
            end
            K(pos(i,j),pos(i,j)) = xin+xis+xiw+xie;
        end
    end
end

%Dirichlet
for ij = Dir
    K(ij,ij) = 1;
    f(ij) = 293;
end

%u
K = sparse(K);
u = K\f;




function m = meank(k1,k2)
    m = (k1+k2)/2; %arithmetic
    %m = 2/(1/k1+1/k2); %harmonic
end


function ij = pos(i,j)
    ij = i+N*(j-1);
end

end