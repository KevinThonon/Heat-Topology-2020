clear all
close all

% gegevens
N = 10;
q = 2/(0.01*0.01*0.001);
rmin = 2;

dp = N+1;

% percentage metaal
pctmetal = 0.3*ones(N,N);

loop1 = 0;
change = 1.;
while change > 0.01
    loop1 = loop1 +1
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

    %cost = cost1(T, dp);

    lambda_vec = lambda1(T, K, dp);
    dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
    dcda_mat = check(N, rmin, pctmetal, dcda_mat1);

    pctmetal = OC(N, pctmetal, 0.4, dcda_mat);

    change = max(max(abs(pctmetal-pctmetal_old)))
end


figure()
surface(pctmetal)

times = 4;
N = times*N;
dp = N+1;
pctmetal = meshrefine(pctmetal,times);

loop2 = loop1;
change = 1.;
while change > 0.02
    loop2 = loop2 +1
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

    %cost = cost2(T, dp);

    lambda_vec = lambda1(T, K, dp);
    dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
    dcda_mat = check(N, rmin, pctmetal, dcda_mat1);

    pctmetal = OC(N, pctmetal, 0.4, dcda_mat);

    change = max(max(abs(pctmetal-pctmetal_old)))
end



% times = 2;
% N = times*N;
% dp = N+1;
% pctmetal = meshrefine(pctmetal,times);
% 
% loop = 0;
% change = 1.;
% while change > 0.01
%     loop = loop +1
%     pctmetal_old = pctmetal;
%     [T, K, f] = fvm_func(pctmetal, N, q);
% 
% %     figure()
% %     spy(K)
% % 
% %     Tmat = zeros(dp,dp);
% %     for i = 1:dp
% %          Tmat(:,i)=T(i*dp-dp+1:i*dp);
% %     end
% % 
% %     figure()
% %     surf(Tmat)
% 
%     %cost = cost2(T, dp);
% 
%     lambda_vec = lambda1(T, K, dp);
%     dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
%     dcda_mat = check(N, rmin, pctmetal, dcda_mat1);
% 
%     pctmetal = OC(N, pctmetal, 0.4, dcda_mat);
% 
%     change = max(max(abs(pctmetal-pctmetal_old)))
% end

figure()
surface(pctmetal)

times = 2;
N = times*N;
dp = N+1;
pctmetal = meshrefine(pctmetal,times);

loop3 = loop2;
change = 1.;
while change > 0.03
    loop3 = loop3 +1
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

    %cost = cost2(T, dp);

    lambda_vec = lambda1(T, K, dp);
    dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
    dcda_mat = check(N, rmin, pctmetal, dcda_mat1);

    pctmetal = OC(N, pctmetal, 0.4, dcda_mat);

    change = max(max(abs(pctmetal-pctmetal_old)))
end

figure()
surface(pctmetal)

times = 10/8;
N = times*N;
dp = N+1;
pctmetal = meshrefine(pctmetal,times);

loop4 = loop3;
change = 1.;
while change > 0.01
    loop4 = loop4 +1
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

    %cost = cost2(T, dp);

    lambda_vec = lambda1(T, K, dp);
    dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
    dcda_mat = check(N, rmin, pctmetal, dcda_mat1);

    pctmetal = OC(N, pctmetal, 0.4, dcda_mat);

    change = max(max(abs(pctmetal-pctmetal_old)))
end

figure()
surface(pctmetal)

%%
heatsink = pctmetal*60 + (1-pctmetal)*0.2;

[T, K, f] = fvm_func(heatsink, N, q);

Tmat = zeros(dp,dp);
for i = 1:dp
     Tmat(:,i)=T(i*dp-dp+1:i*dp);
end

figure()
surface(dcda_mat1)


figure()
surface(Tmat)

figure()
surface(heatsink)

figure()
surface(pctmetal)

% N = 100 -> 179 loops nodig





