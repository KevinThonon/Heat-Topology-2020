function [der] = dcda(lambda, T, pctmetal, N)
dcdam = zeros(N,N);
dcdk = zeros(N*N,1);
penal = 3.0;

%dcda element per element opvullen

for i = 1:1:N
    for j = 1:1:N
        dKdk_u = zeros((N+1)*(N+1),1);
        dKdk_u(i + (j-1)*(N+1)) = -1.0*T(i + (j-1)*(N+1)) + mean(T(i+1 + (j-1)*(N+1)), T(i + (j)*(N+1)));
        dKdk_u(i+1 + (j-1)*(N+1)) = -1.0*T(i+1 + (j-1)*(N+1)) + mean(T(i + (j-1)*(N+1)), T(i+1 + (j)*(N+1)));
        dKdk_u(i + (j)*(N+1)) = -1.0*T(i + (j)*(N+1)) + mean(T(i + (j-1)*(N+1)), T(i+1 + (j)*(N+1)));
        dKdk_u(i+1 + (j)*(N+1)) = -1.0*T(i+1 + (j)*(N+1)) + mean(T(i+1 + (j-1)*(N+1)), T(i + (j)*(N+1)));
        
        dcdk(i+(j-1)*N) = transpose(lambda)*dKdk_u;
        dcdam(i,j) = penal*(65.0-0.2)*pow(pctmetal(i,j),penal-1)*dcdk(i + (j-1)*N);
        
        
    end
end

der = dcdam;


end

function [sol] = pow(a,b)
sol = a^b;
end
