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

loglog(N, resnorm)

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

loglog(N, resnorm)


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



function [solmat] = fvmfunc(N,i)
h = 0.01/N; 
dp = N+1;
kmat = ones(N,N); 
[LL,RL] = FVM_mat(i,h,dp,kmat);
[solmat,~] = solver(dp,LL,RL);
end


function [ll,RL] = FVM_mat(j,h,dp,k)
%LL*T = RL
% T is de kolomvector met de te berekenen variabelen:
% [HP1,N,N,D,D,D,N,N,HP2,N,X,X,X,X,X,X,X,N,N,X,X,X,X,X,X,X,N,...,HP3,N,N,D,D,D,N,N,HP4]

if (j==1)
    q=2;
    qhp = -q*h^2/4;
    qn = -q*h^2/2;
    qq = -q*h^2;
    %Rechterlid
    col11 = qn*ones(dp,1);
    col11(1,1) = qhp;
    col11(dp,1) = qhp;
    col12 = col11;
    counter2 = 0;
    %DirichletRVW in RL
    for i=1:dp
            col11(i) = 0;
            counter2 = counter2+1;
    end
    for i=1:dp    
            col12(i) = 0;
    end
    col2 = qq*ones(dp,1);
    col2(1,1) = qn;
    col2(dp,1) = qn;

    RL = repmat(col2,dp,1);
    RL(1:dp,1) = col11;
    RL((dp*(dp-1)+1):end,1)=col12;
    %col1 hardcoden in LL
elseif (j==2)
    N = dp-1;
    RL = x3func(N,h);
elseif (j==3)
    N = dp-1;
    RL = x4func(N,h);
elseif (j==4)
    N = dp-1;
    RL = x5func(N,h);
end

%ll = zeros(dp^2,dp^2);
ll = sparse(dp^2,dp^2);
N = dp-1;

%dirichlet
for i = (0.3*N+1):1:(0.7*N+1)
    ll(i,i) = 1;
    ll(dp^2-i+1, dp^2-i+1) = 1;
end

%hoekpunten
%linksboven
ll(1,1) = -mean(k(1,1), k(1,1));
ll(1,2) = 0.5*k(1,1);
ll(1,dp+1) = 0.5*k(1,1);
%linksonder
ll(dp,dp) = -mean(k(N,1), k(N,1));
ll(dp,dp-1) = 0.5*k(N,1);
ll(dp,2*dp) = 0.5*k(N,1);
%rechtsboven
ll((dp-1)*dp+1,(dp-1)*dp+1) = -mean(k(1,N), k(1,N));
ll((dp-1)*dp+1,(dp-1)*dp+2) = 0.5*k(1,N);
ll((dp-1)*dp+1,dp*(dp-2)+1) = 0.5*k(1,N);
%rechtsonder
ll(dp^2,dp^2) = -mean(k(N,N), k(N,N));
ll(dp^2,(dp-1)*dp) = 0.5*k(N,N);
ll(dp^2,dp^2-1) = 0.5*k(N,N);

%Neumann links&rechts
i = 1;
while i<0.3*N
    ll(i+1,i+1) = -0.5*k(i,1) -0.5*k(i+1,1) - mean(k(i,1),k(i+1,1));
    ll(i+1,i) = 0.5*k(i,1);
    ll(i+1,i+2) = 0.5*k(i+1,1);
    ll(i+1,i+N+2) = mean(k(i,1), k(i+1,1));
    
    ll(N-i+1,N-i+1) = -0.5*k(N-i+1,1) -0.5*k(N-i,1) - mean(k(N-i+1,1), k(N-i,1));
    ll(N-i+1,N-i) = 0.5*k(N-i,1);
    ll(N-i+1,N-i+2) = 0.5*k(N-i+1,1);
    ll(N-i+1,N-i+N+2) = mean(k(N-i+1,1), k(N-i,1));
    
    i = i+1;
end


i = 1;
while i<0.3*N
    ll(N*(N+1)+i+1,N*(N+1)+i+1) = -0.5*k(i,N) -0.5*k(i+1,N) - mean(k(i,N), k(i+1,N));
    ll(N*(N+1)+i+1,N*(N+1)+i) = 0.5*k(i,N);
    ll(N*(N+1)+i+1,N*(N+1)+i+2) = 0.5*k(i+1,N);
    ll(N*(N+1)+i+1,N*(N+1)+i-(N+1)+1) = mean(k(i,N), k(i+1,N));
    
    ll(N*(N+2)-i+1,N*(N+2)-i+1) = -0.5*k(N-i+1,N) -0.5*k(N-1-i+1,N) - mean(k(N-i+1,N), k(N-1-i+1,N));
    ll(N*(N+2)-i+1,N*(N+2)-i) = 0.5*k(N-1-i+1,N);
    ll(N*(N+2)-i+1,N*(N+2)-i+2) = 0.5*k(N-i+1,N);
    ll(N*(N+2)-i+1,N*(N+2)-i-(N+1)+1) = mean(k(N-i+1,N), k(N-1-i+1,N));
    
    i = i+1;
end

%i = 2;
%while i<(N+1)
for i = 2:1:N
    
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+1) = -0.5*k(1,i-1) -0.5*k(1,i) - mean(k(1,i-1), k(1,i));
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+2) = mean(k(1,i-1), k(1,i));
    ll((i-1)*(N+1)+1,(i-2)*(N+1)+1) = 0.5*k(1,i-1);
    ll((i-1)*(N+1)+1,(i)*(N+1)+1) = 0.5*k(1,i);
    
    ll(i*(N+1),i*(N+1)) = -0.5*k(N,i-1) -0.5*k(N,i) - mean(k(N,i-1), k(N,i));
    ll(i*(N+1),i*(N+1)-1) = mean(k(N,i-1), k(N,i));
    ll(i*(N+1),(i-1)*(N+1)) = 0.5*k(N,i-1);
    ll(i*(N+1),(i+1)*(N+1)) = 0.5*k(N,i);
    
end

for j = 2:1:N
    for i = 2:1:N
        
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i) = -(mean(k(i-1,j-1), k(i-1,j)) + mean(k(i,j-1), k(i,j)) + mean(k(i-1,j-1), k(i,j-1)) + mean(k(i-1,j), k(i,j)));
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i-1) = mean(k(i-1,j-1), k(i-1,j));
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i+1) = mean(k(i,j-1),k(i,j));
        ll((j-1)*(N+1)+i,(j-2)*(N+1)+i) = mean(k(i-1,j-1),k(i,j-1));
        ll((j-1)*(N+1)+i,j*(N+1)+i) = mean(k(i-1,j),k(i,j));
    end
end

end





