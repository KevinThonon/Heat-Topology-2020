clear all
close all
%% kwadratische oplossing
% y = -(x-0.01)*x
% invullen in diffvgl geeft oplossing D2T/Dx2 = 2 --> Q = 2

% ana sol

%num sol
N = 10:10:200;
resnorm = zeros(1,20);

for i = 1:20
    n = N(i);
    dp = n+1;
    
    k = linspace(0,0.01,dp);
    K = repmat(k,dp,1);
    Solmat = -(K-0.01).*K;
    Solmat2 = Solmat(:);
   
    fvmsol = fvmfunc(n,1);
    fvmsol = fvmsol(:);
    
    resnorm(i) = (norm(Solmat2-fvmsol))/norm(fvmsol);
end

plot(N, resnorm)


%% 4 oplossing
%y = x^4
%D2T/Dx2 = 12x^2

N = 10:10:200;
resnorm = zeros(1,20);

for i = 1:20
    n = N(i);
    dp = n+1;
    k = linspace(0,0.01,dp);
    K = repmat(k,dp,1);
    Solmat = K.^4;
    Solmat2 = Solmat(:);

    fvmsol = fvmfunc(n,3);
    fvmsol = fvmsol(:);
    resnorm(i) = norm(Solmat2-fvmsol);
end

plot(N, resnorm)
loglog(N, resnorm)

%% x^3

N = 10:10:200;
resnorm = zeros(1,20);

for i = 1:20
    n = N(i);
    dp = n+1;
    k = linspace(0,0.01,dp);
    K = repmat(k,dp,1);
    Solmat = K.^3;
    Solmat2 = Solmat(:);

    fvmsol = fvmfunc(n,2);
    fvmsol = fvmsol(:);
    resnorm(i) = norm(Solmat2-fvmsol)/norm(fvmsol);
end

plot(N, resnorm)

%% x^5

N = 10:10:200;
resnorm = zeros(1,20);

for i = 1:20
    n = N(i);
    dp = n+1;
    k = linspace(0,0.01,dp);
    K = repmat(k,dp,1);
    Solmat = K.^5;
    Solmat2 = Solmat(:);

    fvmsol = fvmfunc(n,4);
    fvmsol = fvmsol(:);
    resnorm(i) = norm(Solmat2-fvmsol);
end

plot(N, resnorm)


%%


fvmsol = fvmfunc(10,4);

dp = 10+1;
    k = linspace(0,0.01,dp);
    K = repmat(k,dp,1);
    Solmat = K.^5;
    Solmat2 = Solmat(:);
    
X = linspace(0,0.01,dp);
Y = linspace(0,0.01,dp);
figure()
surf(X,Y,Solmat)



