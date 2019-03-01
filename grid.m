function [] = grid(N)
K = zeros(N^2);
h = 0.01/N;
kmatrix = zeros(N);

for i = 1:N
    for j = 1:N
        esum = 0;
        if j+1 <= N
            en = condsquare(i,j+1);
            esum = esum + en;
            K(pos(i,j),pos(i,j+1)) = -en;
        end
        if i+1 <= N
            ee = condsquare(i+1,j);
            esum = esum + ee;
            K(pos(i,j),pos(i+1,j)) = -ee;
        end
        if j-1 >= 1
            es = condsquare(i,j-1);
            esum = esum + es;
            K(pos(i,j),pos(i,j-1)) = -es;
        end
        if i-1 >= 1
            ew = condsquare(i-1,j);
            esum = esum + ew;
            K(pos(i,j),pos(i-1,j)) = -ew;
        end
        K(pos(i,j),pos(i,j)) = esum;
        
        kmatrix(i,j) = k(i,j);
    end
end
figure
heatmap(kmatrix)
figure
f = -ones(N^2,1)*0.2*h^2;
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
end
end

function kresult = condsquare(i,j)

k1 = k(i+0.5,j+0.5);
k2 = k(i+0.5,j-0.5);
k3 = k(i-0.5,j-0.5);
k4 = k(i-0.5,j+0.5);

kresult = (k1+k2+k3+k4)/((k1>0)+(k2>0)+(k3>0)+(k4>0)); %arithmetic
%kresult = 4/(1/k1+1/k2+1/k3+1/k4); %harmonic
end


function ij = pos(i,j)
ij = i+N*(j-1);
end

end