function [] = grid(N)
K = zeros(N^2);
f = zeros(N^2,1);
h = 0.01/N;
kmatrix = zeros(N);

for i = 1:N
    for j = 1:N
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
        
        kmatrix(i,j) = k(i,j);
    end
end
figure
heatmap(kmatrix)
figure
f = f-ones(N^2,1)*2*h^2;
K = sparse(K);
u = K\f;

field = zeros(N);
for j = 1:N
    field(:,j) = u(pos(1:N,j));
end
heatmap(field)




function y = k(i,j)
    %pctmetal = 0.4;
    pctmetal = 0.1;
    if j<=0.2*N || j > 0.8*N
        pctmetal = 0.9;
    end  
    y = (65-0.2)*pctmetal^3 + 0.2;
    
    if i<1 || j<1 || i>N || j>N
        y = 0;
        disp('error')
    end
end


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