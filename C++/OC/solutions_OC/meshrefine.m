function [pct] = meshrefine(pctmetal, times)
ny_nx = times*size(pctmetal);
N = ny_nx(1);
pct=zeros(N,N);
for i = 1:N/times
    for j = 1:N/times
        pct(times*i-times+1:times*i,times*j-times+1:times*j) = pctmetal(i,j);
    end
end


end

