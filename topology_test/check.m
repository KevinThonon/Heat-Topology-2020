function [dcn] = check(N, rmin, x, dc)
% x = designvariable(matrix)
% dc = dcda
% --> modify the element sensitivities adhv de formule gegeven in de paper
% gevolg is dat er geen checkerboard is

dcn = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        sum = 0.0;
        for k = max(i-round(rmin),1):1:min(i+round(rmin),N)
            for l = max(j-round(rmin),1):1:min(j+round(rmin),N)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcn(j,i) = dcn(j,i) + max(0, fac)*x(l,k)*dc(l,k);
            end
        end
        dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
    end
end
end 

%Emma typte hier! Philippe vindt Emma leuk!

