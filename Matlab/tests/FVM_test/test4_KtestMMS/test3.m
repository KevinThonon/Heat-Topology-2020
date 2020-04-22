clear all
close all

N=10;
dp = N+1;
j=1;
h = 0.01/10;
dp =N+1;

%k maken


surf(bigKmat(meshgrid(1:1:10),10))

%% k=100/0.01*X

n1 = 10:10:90;
resnorm1 = zeros(9,1);
j=1;

for i =1:1:9
    N = n1(i)
    dp = N+1;
    x = (0.01/N)/2:(0.01/N):0.01;


    k = meshgrid(x)*100/0.01;
    [sol, ~,q] = x4func_1(k,N);
    size(sol);

    figure()
    surf(sol)

    [solution, LL, RL] = fvm_func(j, N, k);

    % solmat2 = zeros(dp,dp);
    % for i = 1:1:dp
    %     solmat2(:,i) = solution((i-1)*dp+1:i*dp);
    % end
    % 
    % figure()
    % surf(solmat2)

    resnorm1(i) = norm(sol(:)-solution)/norm(sol(:));
    
end

n2 = 100:50:550;
resnorm2 = zeros(1,10);
for i =1:1:10
    N = n2(i)
    dp = N+1;
    x = (0.01/N)/2:(0.01/N):0.01;


    k = meshgrid(x)*100/0.01;
    [sol, ~,q] = x4func_1(k,N);
    size(sol);

    figure()
    surf(sol)

    [solution, LL, RL] = fvm_func(j, N, k);

    % solmat2 = zeros(dp,dp);
    % for i = 1:1:dp
    %     solmat2(:,i) = solution((i-1)*dp+1:i*dp);
    % end
    % 
    % figure()
    % surf(solmat2)

    resnorm2(i) = norm(sol(:)-solution)/norm(sol(:));
    
end

resnorm = [resnorm1, resnorm2];
n = [n1, n2];

figure()
loglog(n,resnorm)

figure()
surf(x,x,k)











%% functions 


function bigK = bigKmat(k,N)
bigK = zeros(N+1,N+1);
dp = N+1;
%boven en onderrand; links & rechts; internal
for i = 2:1:N
    bigK(1,i) = (k(1,i-1)+k(1,i))/2;
    bigK(dp,i) = (k(N,i-1)+k(N,i))/2;
    bigK(i,1) = (k(i-1,1)+k(i,1))/2;
    bigK(i,dp) = (k(i-1,N)+k(i,N))/2;
    
    for j = 2:1:N
        bigK(i,j) = (k(i-1,j-1)+k(i-1,j)+k(i,j-1)+k(i,j))/4;
    end
end

bigK(1,1) = k(1,1);
bigK(1,dp) = k(1,N);
bigK(dp,1) = k(N,1);
bigK(dp,dp) = k(N,N);

end

function [solution, LL, RL] = fvm_func(j, N,k)
h = 0.01/N;
dp = N+1;

[LL,RL] = FVM_mat(j,dp,k);
solution = LL\RL;
end


function k = mean(a, b)
k = (a+b)/2;           % arithmetic mean
%k = 2.0*(a*b)/(a+b);    % harmonic mean
end


function [ll,RL] = FVM_mat(j,dp,k)
%LL*T = RL
% T is de kolomvector met de te berekenen variabelen:
% [HP1,N,N,D,D,D,N,N,HP2,N,X,X,X,X,X,X,X,N,N,X,X,X,X,X,X,X,N,...,HP3,N,N,D,D,D,N,N,HP4]
N=dp-1;

    
if (j==1)
    [~,RL,~] = x4func_1(k,N);
elseif (j==2)
    [~,RL,~] = x4func_2(k,N);
elseif(j==3)
    [~,RL,~] = x4func_11(k,N);
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




function [sol, RL,q] = x4func_1(k,N)

h = 0.01/N;
x = 0:0.01/N:0.01;
y = 0:0.01/N:0.01;

[X,~] = meshgrid(x,y);
sol = X.^4;

%sol diff eq

K = bigKmat(k,N);

%[Kx,Ky] = gradient(K);
Kx = 100/0.01;
q=Kx*(4*X.^3) +K.*(12*X.^2);


for i = 2:(N)
    for j = 2:(N)
        q(i,j) = q(i,j)*h^2;
    end
end
 
for i = 2:(N)
    q(1,i) = q(1,i)*h^2/2;
    q((N+1),i) = q((N+1),i)*h^2/2;
end


q(1:(N+1),1) = sol(1:(N+1),1);
q(1:(N+1),(N+1)) = sol(1:(N+1),(N+1));


RL = q(:);


end

function [sol, RL,q] = x4func_11(k,N)

h = 0.01/N;
x = 0:0.01/N:0.01;
y = 0:0.01/N:0.01;

[X,~] = meshgrid(x,y);
sol = X.^4;

%sol diff eq

K = bigKmat(k,N);

%[Kx,Ky] = gradient(K);
Ky = 100/0.01;
q=Ky*(4*X.^3) +K.*(12*X.^2);


for i = 2:(N)
    for j = 2:(N)
        q(i,j) = q(i,j)*h^2;
    end
end
 
for i = 2:(N)
    q(1,i) = q(1,i)*h^2/2;
    q((N+1),i) = q((N+1),i)*h^2/2;
end


q(1:(N+1),1) = sol(1:(N+1),1);
q(1:(N+1),(N+1)) = sol(1:(N+1),(N+1));

RL = q(:);


end


function [sol, RL,q] = x4func_2(k,N)
%-10x(x-0.01)


h = 0.01/N;
x = 0:0.01/N:0.01;
y = 0:0.01/N:0.01;

[X,~] = meshgrid(x,y);
sol = X.^4;

%sol diff eq

K = bigKmat(k,N);

%[Kx,Ky] = gradient(K);
Kx = -20*X+0.1;
q=Kx*(4*X.^3) +K.*(12*X.^2);


for i = 2:(N)
    for j = 2:(N)
        q(i,j) = q(i,j)*h^2;
    end
end
 
for i = 2:(N)
    q(1,i) = q(1,i)*h^2/2;
    q((N+1),i) = q((N+1),i)*h^2/2;
end


q(1:(N+1),1) = sol(1:(N+1),1);
q(1:(N+1),(N+1)) = sol(1:(N+1),(N+1));


RL = q(:);


end






