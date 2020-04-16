clear all
close all


%% kwadratische functie
% y = -(x-0.01)*x
% method of manufactured solutions: D2T/Dx2 = 2 --> Q = 2
j = 2;


N = 10:10:200;
resnorm = zeros(1,20);

for i=1:1:20
    n = N(i);
    h = 0.01/n;
    dp=n+1;
    [Solmat,~] = x2func(dp,h);
    Solmat2 = zeros(dp,dp);
%     for i = 1:1:dp
%         Solmat2(:,i) = Solmat((i-1)*dp+1:i*dp);
%     end

    [fvmsol,~,~] = fvm_func(j,n);
%     Tmat = zeros(dp,dp);
%     for i = 1:1:dp
%         Tmat(:,i) = fvmsol((i-1)*dp+1:i*dp);
%     end
    
    resnorm(i) = norm(Solmat-fvmsol);
end

loglog(N,resnorm)

%% Kubische functie
% X^3
j = 3;
N = 10:10:200;
resnorm = zeros(1,20);

for i=1:1:20
    n = N(i);
    h=0.01/n;
    dp=n+1;
    [Solmat,~] = x3func(dp,h);
    Solmat2 = zeros(dp,dp);
%     for i = 1:1:dp
%         Solmat2(:,i) = Solmat((i-1)*dp+1:i*dp);
%     end

    [fvmsol,~,~] = fvm_func(j,n);
%     Tmat = zeros(dp,dp);
%     for i = 1:1:dp
%         Tmat(:,i) = fvmsol((i-1)*dp+1:i*dp);
%     end
    
    resnorm(i) = norm(Solmat-fvmsol);
end
figure()
loglog(N,resnorm)
%% X^4
j = 4;
N = 50:50:1000;
resnorm = zeros(1,20);

for i=1:1:20
    i
    i*50
    
    n = N(i);
    h=0.01/n;
    dp=n+1;
    [Solmat,~] = x4func(dp,h);
    Solmat2 = zeros(dp,dp);
%     for i = 1:1:dp
%         Solmat2(:,i) = Solmat((i-1)*dp+1:i*dp);
%     end

    [fvmsol,~,~] = fvm_func(j,n);
%     Tmat = zeros(dp,dp);
%     for i = 1:1:dp
%         Tmat(:,i) = fvmsol((i-1)*dp+1:i*dp);
%     end
    
    resnorm(i) = norm(Solmat-fvmsol);
end
figure()
loglog(N,resnorm)


%% X^5
j = 5;
N = 10:10:200;
resnorm = zeros(1,20);

for i=1:1:20
    n = N(i);
    h=0.01/n;
    dp=n+1;
    [Solmat,~] = x5func(dp,h);
    Solmat2 = zeros(dp,dp);
%     for i = 1:1:dp
%         Solmat2(:,i) = Solmat((i-1)*dp+1:i*dp);
%     end

    [fvmsol,~,~] = fvm_func(j,n);
%     Tmat = zeros(dp,dp);
%     for i = 1:1:dp
%         Tmat(:,i) = fvmsol((i-1)*dp+1:i*dp);
%     end
    
    resnorm(i) = norm(Solmat-fvmsol);
end
figure()
loglog(N,resnorm)



%%
function [solution, LL, RL] = fvm_func(j, N)
h = 0.01/N;
dp = N+1;

[LL,RL] = FVM_mat(j,h,dp);
solution = LL\RL;
end


function [K] = createK(N, penal)
%dichtheid k in elk element met simp methode
pctmetal = 0.3*ones(N,N);
K = zeros(N,N);
for i = 1:1:N
    for j = 1:1:N
        K(i,j) = 64.8*(pctmetal(i,j)^penal) + 0.2;
    end
end

end

function k = mean(a, b)
k = (a+b)/2;           % arithmetic mean
%k = 2.0*(a*b)/(a+b);    % harmonic mean
end


function [ll,RL] = FVM_mat(j,h,dp)
%LL*T = RL
% T is de kolomvector met de te berekenen variabelen:
% [HP1,N,N,D,D,D,N,N,HP2,N,X,X,X,X,X,X,X,N,N,X,X,X,X,X,X,X,N,...,HP3,N,N,D,D,D,N,N,HP4]
N=dp-1;

    
if (j==1)
    k=createK(dp-1,3);
    q = 2/(0.01*0.01*0.001);
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
    
elseif (j==2)
    k=ones(N,N);
    [~,RL] = x2func(dp,h);
elseif (j==3)
    k=ones(N,N);
    [~,RL] = x3func(dp,h);
elseif (j==4)
    k=ones(N,N);
    [~,RL] = x4func(dp,h);
    
elseif (j==5)
    k=ones(N,N);
    [~,RL] = x5func(dp,h);
end
   




%ll = zeros(dp^2,dp^2);
ll = sparse(dp^2,dp^2);
N = dp-1;

%dirichlet
for i = 1:1:dp%1:1:(0.7*N+1)
    ll(i,i) = 1;
    ll(dp^2-i+1, dp^2-i+1) = 1;
end

% %hoekpunten
% %linksboven
% ll(1,1) = -mean(k(1,1), k(1,1));
% ll(1,2) = 0.5*k(1,1);
% ll(1,dp+1) = 0.5*k(1,1);
% %linksonder
% ll(dp,dp) = -mean(k(N,1), k(N,1));
% ll(dp,dp-1) = 0.5*k(N,1);
% ll(dp,2*dp) = 0.5*k(N,1);
% %rechtsboven
% ll((dp-1)*dp+1,(dp-1)*dp+1) = -mean(k(1,N), k(1,N));
% ll((dp-1)*dp+1,(dp-1)*dp+2) = 0.5*k(1,N);
% ll((dp-1)*dp+1,dp*(dp-2)+1) = 0.5*k(1,N);
% %rechtsonder
% ll(dp^2,dp^2) = -mean(k(N,N), k(N,N));
% ll(dp^2,(dp-1)*dp) = 0.5*k(N,N);
% ll(dp^2,dp^2-1) = 0.5*k(N,N);

%Neumann links&rechts
% i = 1;
% while i<0.3*N
%     ll(i+1,i+1) = -0.5*k(i,1) -0.5*k(i+1,1) - mean(k(i,1),k(i+1,1));
%     ll(i+1,i) = 0.5*k(i,1);
%     ll(i+1,i+2) = 0.5*k(i+1,1);
%     ll(i+1,i+N+2) = mean(k(i,1), k(i+1,1));
%     
%     ll(N-i+1,N-i+1) = -0.5*k(N-i+1,1) -0.5*k(N-i,1) - mean(k(N-i+1,1), k(N-i,1));
%     ll(N-i+1,N-i) = 0.5*k(N-i,1);
%     ll(N-i+1,N-i+2) = 0.5*k(N-i+1,1);
%     ll(N-i+1,N-i+N+2) = mean(k(N-i+1,1), k(N-i,1));
%     
%     i = i+1;
% end
% 
% 
% i = 1;
% while i<0.3*N
%     ll(N*(N+1)+i+1,N*(N+1)+i+1) = -0.5*k(i,N) -0.5*k(i+1,N) - mean(k(i,N), k(i+1,N));
%     ll(N*(N+1)+i+1,N*(N+1)+i) = 0.5*k(i,N);
%     ll(N*(N+1)+i+1,N*(N+1)+i+2) = 0.5*k(i+1,N);
%     ll(N*(N+1)+i+1,N*(N+1)+i-(N+1)+1) = mean(k(i,N), k(i+1,N));
%     
%     ll(N*(N+2)-i+1,N*(N+2)-i+1) = -0.5*k(N-i+1,N) -0.5*k(N-1-i+1,N) - mean(k(N-i+1,N), k(N-1-i+1,N));
%     ll(N*(N+2)-i+1,N*(N+2)-i) = 0.5*k(N-1-i+1,N);
%     ll(N*(N+2)-i+1,N*(N+2)-i+2) = 0.5*k(N-i+1,N);
%     ll(N*(N+2)-i+1,N*(N+2)-i-(N+1)+1) = mean(k(N-i+1,N), k(N-1-i+1,N));
%     
%     i = i+1;
% end

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


function [mat,der2] = x2func(dp,h)
k = linspace(0,0.01,dp);
K = repmat(k,dp,1);

%gewone oplossing
mat = -(K-0.01).*K;
mat = mat(:);

%2e afgeleide
der2 = 2*ones(dp,dp);

for i = 2:(dp-1)
    for j = 2:(dp-1)
        der2(i,j) = -der2(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    der2(1,i) = -der2(1,i)*h^2/2;
    der2(dp,i) = -der2(dp,i)*h^2/2;
end

der2(1:dp,1) = zeros(dp,1);
der2(1:dp,dp) = zeros(dp,1);

der2 = der2(:);

end


function [mat,der2] = x3func(dp,h)

k = linspace(0,0.01,dp);
K = repmat(k,dp,1);
%gewone oplossing

mat = K.^3;
mat = mat(:);

der2 = -6*K;
for i = 2:(dp-1)
    for j = 2:(dp-1)
        der2(i,j) = -der2(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    der2(1,i) = -der2(1,i)*h^2/2;
    der2(dp,i) = -der2(dp,i)*h^2/2;
end

der2(1:dp,1) = zeros(1,dp);
der2(1:dp,dp) = ones(1,dp)*0.01^3;
der2 = der2(:);

end



function [mat,der2] = x4func(dp,h)
k = linspace(0,0.01,dp);
K = repmat(k,dp,1);
%gewone oplossing

mat = K.^4;
mat = mat(:);

der2 = -12*K.^2;
for i = 2:(dp-1)
    for j = 2:(dp-1)
        der2(i,j) = -der2(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    der2(1,i) = -der2(1,i)*h^2/2;
    der2(dp,i) = -der2(dp,i)*h^2/2;
end

der2(1:dp,1) = zeros(1,dp);
der2(1:dp,dp) = ones(1,dp)*0.01^4;
der2 = der2(:);
end

function [mat,der2] = x5func(dp,h)
k = linspace(0,0.01,dp);
K = repmat(k,dp,1);
%gewone oplossing

mat = K.^5;
mat = mat(:);

der2 = -20*K.^3;
for i = 2:(dp-1)
    for j = 2:(dp-1)
        der2(i,j) = -der2(i,j)*h^2;
    end
end

for i = 2:(dp-1)
    der2(1,i) = -der2(1,i)*h^2/2;
    der2(dp,i) = -der2(dp,i)*h^2/2;
end

der2(1:dp,1) = zeros(1,dp);
der2(1:dp,dp) = ones(1,dp)*0.01^5;
der2 = der2(:);
end








