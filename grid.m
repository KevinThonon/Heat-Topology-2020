function [] = grid(N)
K = zeros(N^2);
h = 0.01/N; 
f = ones(N^2,1)*2*h^2 *10^6; %fout om verschil te tonen
pctmetal = [.4*ones(N,.4*N) zeros(N,.6*N)];
k = (65-0.2)*pctmetal + 0.2;

Dir = [round(.3*N):round(.7*N) ...
    pos(1,N)+(round(.3*N):round(.7*N))];

for i = 1:N
    for j = 1:N
        if ~ismember(pos(i,j), Dir)
            [xin,xis,xiw,xie] = xi(i,j);
            if i > 1
                K(pos(i,j),pos(i-1,j)) = -xin;
            end
            if i < N
                K(pos(i,j),pos(i+1,j)) = -xis;
            end
            if j > 1
                K(pos(i,j),pos(i,j-1)) = -xiw;
            end
            if j < N
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

%plot
field = zeros(N);
for j = 1:N
    field(:,j) = u(pos(1:N,j));
end
figure
heatmap(field)
figure
heatmap(k)




function [xin,xis,xiw,xie] = xi(i,j)
    xin=0; xis=0; xiw=0; xie=0;
    if i > 1
        xin = meank(k(i,j),k(i-1,j));
    end
    if i < N
        xis = meank(k(i,j),k(i+1,j));
    end
    if j > 1
        xiw = meank(k(i,j),k(i,j-1));
    end
    if j < N
        xie = meank(k(i,j),k(i,j+1));
    end
end


function m = meank(k1,k2)
    m = (k1+k2)/2; %arithmetic
    %m = 2/(1/k1+1/k2); %harmonic
end


function ij = pos(i,j)
    ij = i+N*(j-1);
end

end
