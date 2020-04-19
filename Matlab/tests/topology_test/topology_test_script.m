clear all
close all

% gegevens
N = 10;
q = 2/(0.01*0.01*0.001);
rmin = 2;

dp = N+1;

% percentage metaal
pctmetal = 0.4*ones(N,N);

loop = 0;
change = 1.;
while change > 0.01
    loop = loop +1
    pctmetal_old = pctmetal;


    [T, K, f] = fvm_func(pctmetal, N, q);

%     figure()
%     spy(K)
% 
%     Tmat = zeros(dp,dp);
%     for i = 1:dp
%          Tmat(:,i)=T(i*dp-dp+1:i*dp);
%     end
% 
%     figure()
%     surf(Tmat)

    cost = costfunc(T, dp);

    lambda_vec = lambda(T, K, dp);
    dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
    dcda_mat = check(N, rmin, pctmetal, dcda_mat1);

    pctmetal = OC(N, pctmetal, 0.4, dcda_mat);

    change = max(max(abs(pctmetal-pctmetal_old)));
end

[T, K, f] = fvm_func(pctmetal, N, q);

Tmat = zeros(dp,dp);
for i = 1:dp
     Tmat(:,i)=T(i*dp-dp+1:i*dp);
end

figure()
surf(Tmat)

% N = 100 -> 179 loops nodig





