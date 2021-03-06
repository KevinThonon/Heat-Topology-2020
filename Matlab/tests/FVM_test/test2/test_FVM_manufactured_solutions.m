clear all
close all

% N = 500;
% h = 0.01/N;
% dp = N+1;
% q = 2/(0.01*0.01*0.001);
% pctmetal = 0.3*ones(N,N);
% k = createK(pctmetal,N,3);
% [sol, LL, RL] = fvm_func(k, N, q);
% 
% spy(LL)

%


%%
N = 100;
X = 0:0.01/N:0.01;
y = X;
dp = N+1;
q = 2/(0.01*0.01*0.001);
kmat = 0.2*ones(N,N);
kmat(30:70,:) = 100000;
%kmat = ones(N,N);

[solution, LL, RL] = fvm_func(kmat, N, q);
Tmat = zeros(dp,dp);
for i = 1:1:dp
    Tmat(:,i) = solution((i-1)*dp+1:i*dp);
end

figure()
surf(X,X,Tmat)
xlabel('x')
ylabel('y')

figure()
surf(kmat)
%%
N = 10;
X = 0:0.01/N:0.01;
y = X;
dp = N+1;
q = 0;
kmat = 50*ones(N,N);

[solution, LL, RL] = fvm_func(kmat, N, q);
Tmat = zeros(dp,dp);
for i = 1:1:dp
    Tmat(:,i) = solution((i-1)*dp+1:i*dp);
end

figure()
surf(X,X,Tmat)
xlabel('x')
ylabel('y')


%%
N = 100;
X = 0:0.01/N:0.01;
y = X;
dp = N+1;
q = 2/(0.01*0.01*0.001);
kmat = 0.2*ones(N,N);
kmat(:,1:10) = 10000000;
kmat(:,91:100) = 10000000;

[solution, LL, RL] = fvm_func(kmat, N, q);
Tmat = zeros(dp,dp);
for i = 1:1:dp
    Tmat(:,i) = solution((i-1)*dp+1:i*dp);
end

figure()
surf(X,X,Tmat)
xlabel('x')
ylabel('y')






%%

function [K] = createK(pctmetal, N, penal)
%dichtheid k in elk element met simp methode
K = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        K(i,j) = 64.8*(pctmetal(i,j)^penal) + 0.2;
    end
end

end

function [solution, LL, RL] = fvm_func(kmat, N, q)
h = 0.01/N;
dp = N+1;

[LL,RL] = FVM_mat(q,h,dp,kmat);
solution = LL\RL;
end


function [ll,RL] = FVM_mat(q,h,dp,k)
%LL*T = RL
% T is de kolomvector met de te berekenen variabelen:
% [HP1,N,N,D,D,D,N,N,HP2,N,X,X,X,X,X,X,X,N,N,X,X,X,X,X,X,X,N,...,HP3,N,N,D,D,D,N,N,HP4]


qhp = -q*h^2/4;
qn = -q*h^2/2;
qq = -q*h^2;

%Rechterlid
col1 = qn*ones(dp,1);
col1(1,1) = qhp;
col1(dp,1) = qhp;

counter2 = 0;
%DirichletRVW in RL
for i=1:dp
    if(h*(i-1)>=0.003)&&(h*(i-1)<=0.007)
        col1(i) = 293;
        counter2 = counter2+1;
    end
end

col2 = qq*ones(dp,1);
col2(1,1) = qn;
col2(dp,1) = qn;

RL = repmat(col2,dp,1);
RL(1:dp,1) = col1;
RL((dp*(dp-1)+1):end,1)=col1;



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


function k = mean(a, b)
%k = (a+b)/2;           % arithmetic mean
k = 2.0*(a*b)/(a+b);    % harmonic mean
end