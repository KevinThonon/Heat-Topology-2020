function [K] = createK(pctmetal, N)
%dichtheid k in elk element met simp methode
penal = 3;
K = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        K(i,j) = (65-0.2)*(pctmetal(i,j)^penal) + 0.2;
    end
end

end