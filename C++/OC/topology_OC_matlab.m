clear all
close all

% gegevens
N = 30;
q = 2/(0.01*0.01*0.001);
rmin = 2;

dp = N+1;

% percentage metaal
pctmetal = 0.3*ones(N,N);
penal = 3.0;

loop1 = 0;
change = 1.0;
while loop1 < 30 %change > 0.01
    loop1 = loop1 + 1
        
    pctmetal_old = pctmetal;
    
    [k] = createK(pctmetal, N, penal);
    [T, K, f] = fvm_func(k, N, q);
    
    lambda_vec = lambda1(T, K, N);
    dcda_mat1 = dcda_fd(T, pctmetal, N, penal);
    %dcda_mat = check(N, rmin, pctmetal, dcda_mat1);
    
    pctmetal = OC(N, pctmetal, 0.4, dcda_mat1);
    
    change = max(max(abs(pctmetal-pctmetal_old)));
end


figure()
pctmetal = meshrefine(pctmetal,2);
surface(pctmetal)


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
    ll(i+1,i+1+1) = 0.5*k(i+1,1);
    ll(i+1,i+N+1+1) = mean(k(i,1), k(i+1,1));
    
    ll(N-i+1,N-i+1) = -0.5*k(N-i+1,1) -0.5*k(N-i,1) - mean(k(N-i+1,1), k(N-i,1));
    ll(N-i+1,N-i-1+1) = 0.5*k(N-i,1);
    ll(N-i+1,N-i+1+1) = 0.5*k(N-i+1,1);
    ll(N-i+1,N-i+N+1+1) = mean(k(N-i+1,1), k(N-i,1));
    
    i = i+1;
end


i = 1;
while i<0.3*N
    ll(N*(N+1)+i+1,N*(N+1)+i+1) = -0.5*k(i-1+1,N-1+1) -0.5*k(i+1,N-1+1) - mean(k(i-1+1,N-1+1), k(i+1,N-1+1));
    ll(N*(N+1)+i+1,N*(N+1)+i-1+1) = 0.5*k(i-1+1,N-1+1);
    ll(N*(N+1)+i+1,N*(N+1)+i+1+1) = 0.5*k(i+1,N-1+1);
    ll(N*(N+1)+i+1,N*(N+1)+i-(N+1)+1) = mean(k(i-1+1,N-1+1), k(i+1,N-1+1));
    
    ll(N*(N+2)-i+1,N*(N+2)-i+1) = -0.5*k(N-i+1,N-1+1) -0.5*k(N-1-i+1,N-1+1) - mean(k(N-i+1,N-1+1), k(N-1-i+1,N-1+1));
    ll(N*(N+2)-i+1,N*(N+2)-i-1+1) = 0.5*k(N-1-i+1,N-1+1);
    ll(N*(N+2)-i+1,N*(N+2)-i+1+1) = 0.5*k(N-i+1,N-1+1);
    ll(N*(N+2)-i+1,N*(N+2)-i-(N+1)+1) = mean(k(N-i+1,N-1+1), k(N-1-i+1,N-1+1));
    
    i = i+1;
end

%i = 2;
%while i<(N+1)
for i = 2:1:N
    
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+1) = -0.5*k(1,i-1) -0.5*k(1,i) - mean(k(1,i-1), k(1,i));
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+1+1) = mean(k(1,i-1), k(1,i));
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




% %col1 hardcoden in LL
%
% %dirichletRVW in LL
% LL = zeros(dp^2,dp^2);
% dirll = eye(counter2);
% LL(dp*dp-dp+1:dp*dp,dp*dp-dp+1:dp*dp)=dirll;
% LL(1:dp,1:dp) = dirll;
%
% %hoekpunten
% %linskboven (1,1), LL(1,:)
% LL(1,1) = -0.5*(bigkmat(2,1)+bigkmat(1,2));
% LL(1,2) = 0.5* bigkmat(2,1);%south
% LL(1,1+dp) = 0.5*bigkmat(1,2);%east
% %linksonder (dp,1), LL(dp,:)
% LL(dp,dp) = -0.5*(bigkmat((2*dp-1)-1,1)+bigkmat(2*dp-1,2));
% LL(dp,dp-1) = 0.5*bigkmat((2*dp-1)-1,1);%north
% LL(dp,dp+dp) = 0.5*bigkmat(2*dp-1,2);%east
% %rechtsboven (1,dp),
% LL(dp*dp-dp+1,dp*dp-dp+1) = -0.5*(bigkmat(1,dp*2-1-1)+bigkmat(2,dp*2-1));
% LL(dp*(dp-1)+1,dp*(dp-1)+1+1)=0.5*bigkmat(2,dp*2-1);%south
% LL(dp*(dp-1)+1,dp*dp-dp+1-dp)= 0.5*bigkmat(1,dp*2-1-1);%west
% %rechtsonder (dp,dp), LL(dp*dp,:)
% LL(dp*dp,dp*dp) = -0.5*(bigkmat(dp*2-1,dp*2-1-1)+bigkmat(dp*2-1-1,dp*2-1));
% LL(dp*dp,dp*dp-1) = 0.5*bigkmat(dp*2-1-1,dp*2-1);%north
% LL(dp*dp,dp*dp-dp)=0.5*bigkmat(dp*2-1,dp*2-1-1);%west
%
% %Neumann links & rechts
% for j = 2:(dp-1)
%     if(h*(j-1)<0.003)||(h*(j-1)>0.007)
%         %LL(1:dp,1:2*dp)
%         LL(j,j) = -(bigkmat(j*2-1,2)+0.5*bigkmat(j*2-1-1,1)+0.5*bigkmat(j*2-1+1,1));%east,north,south
%         LL(j,j+dp) = bigkmat(j*2-1,2); %east
%         LL(j,j-1) = 0.5*bigkmat(j*2-1-1,1);%north
%         LL(j,j+1) = 0.5*bigkmat(j*2-1+1,1);%south
%
%         %LL(dp*(dp-1)+1:dp*dp,dp*(dp-2)+1:dp*dp)
%         LL(dp*dp-dp+j,dp*dp-dp+j) = -(bigkmat(j*2-1,dp*2-1-1)+0.5*bigkmat(j*2-1-1,2*dp-1)+0.5*bigkmat(j*2-1+1,2*dp-1));%west,north,south
%         LL(dp*dp-dp+j,dp*dp-dp+j-dp) = bigkmat(j*2-1,dp*2-1-1); %west
%         LL(dp*dp-dp+j,dp*dp-dp+j-1) = 0.5*bigkmat(j*2-1-1,2*dp-1);%north
%         LL(dp*dp-dp+j,dp*dp-dp+j+1) = 0.5*bigkmat(j*2-1+1,2*dp-1);%south
%
%     end
% end
%
% %Neumann boven & onder
%
% for j = 1:(dp-2)
%     k = j+1;
%     %boven
%     LL(j*dp+1,j*dp+1) = -(bigkmat(2,2*k-1)+0.5*bigkmat(1,2*k-1-1)+0.5*bigkmat(1,2*k-1+1));%south,west,east
%     LL(j*dp+1,j*dp+1+1) =bigkmat(2,2*k-1);%south
%     LL(j*dp+1,j*dp+1-dp) = 0.5*bigkmat(1,2*k-1-1);%west
%     LL(j*dp+1,j*dp+1+dp) = 0.5*bigkmat(1,2*k-1+1);%east
%
%     %onder
%     LL(j*dp+dp,j*dp+dp)= -(bigkmat(2*dp-1-1,2*k-1)+0.5*bigkmat(2*dp-1,2*k-1-1)+0.5*bigkmat(2*dp-1,2*k-1+1));%north,west,east
%     LL(j*dp+dp,j*dp+dp-1) = bigkmat(2*dp-1-1,2*k-1);%north
%     LL(j*dp+dp,j*dp+dp-dp) = 0.5*bigkmat(2*dp-1,2*k-1-1);%west
%     LL(j*dp+dp,j*dp+dp+dp) = 0.5*bigkmat(2*dp-1,2*k-1+1);%east
% end
%
% %gewonekolommen
%
% for i = 2:(dp-1)
%     for j = 2:(dp-1)
%         k = (i-1)*dp+j;
%         i2 = 2*i-1;
%         j2 = 2*j-1;
%         LL(k,k) = -(bigkmat(i2+1,j2)+bigkmat(i2-1,j2)+bigkmat(i2,j2+1)+bigkmat(i2,j2-1));%east,west,south,north
%         LL(k,k+dp)=bigkmat(i2+1,j2);%east
%         LL(k,k-dp)=bigkmat(i2-1,j2);%west
%         LL(k,k+1)=bigkmat(i2,j2+1);%south
%         LL(k,k-1)=bigkmat(i2,j2-1);%north
%     end
% end
%
% Z = LL-ll
% normZ = norm(Z)
%
% ll=LL;

end


function k = mean(a, b)
k = (a+b)/2;           % arithmetic mean
%k = 2.0*(a*b)/(a+b);    % harmonic mean
end



function [lambda] = lambda1(T, K, N)
lambda = transpose(K)\((-2.0/N^2)*T);
end


function [der] = dcda(lambda, T, pctmetal, N)
dcdam = zeros(N,N);
dcdk = zeros(N*N,1);
penal = 3.0;

%dcda element per element opvullen

for i = 1:1:N
    for j = 1:1:N
        dKdk_u = zeros((N+1)*(N+1),1);
        
        dKdk_u(i + (j-1)*(N+1)) = -1.0*T(i + (j-1)*(N+1)) + 0.5*T(i+1 + (j-1)*(N+1)) + 0.5*T(i + (j)*(N+1));
        dKdk_u(i+1 + (j-1)*(N+1)) = -1.0*T(i+1 + (j-1)*(N+1)) + 0.5*T(i + (j-1)*(N+1)) + 0.5*T(i+1 + (j)*(N+1));
        dKdk_u(i + (j)*(N+1)) = -1.0*T(i + (j)*(N+1)) + 0.5*T(i + (j-1)*(N+1)) + 0.5*T(i+1 + (j)*(N+1));
        dKdk_u(i+1 + (j)*(N+1)) = -1.0*T(i+1 + (j)*(N+1)) + 0.5*T(i+1 + (j-1)*(N+1)) + 0.5*T(i + (j)*(N+1));
        
        %dKdk_u(i + (j-1)*(N+1)) = -1.0*T(i + (j-1)*(N+1)) + mean(T(i+1 + (j-1)*(N+1)), T(i + (j)*(N+1)));
        %dKdk_u(i+1 + (j-1)*(N+1)) = -1.0*T(i+1 + (j-1)*(N+1)) + mean(T(i + (j-1)*(N+1)), T(i+1 + (j)*(N+1)));
        %dKdk_u(i + (j)*(N+1)) = -1.0*T(i + (j)*(N+1)) + mean(T(i + (j-1)*(N+1)), T(i+1 + (j)*(N+1)));
        %dKdk_u(i+1 + (j)*(N+1)) = -1.0*T(i+1 + (j)*(N+1)) + mean(T(i+1 + (j-1)*(N+1)), T(i + (j)*(N+1)));
        
        dcdk(i+(j-1)*N) = transpose(lambda)*dKdk_u;
        dcdam(i,j) = penal*(65.0-0.2)*pow(pctmetal(i,j),penal-1)*dcdk(i + (j-1)*N);
        
        
    end
end

der = dcdam;


end

function [sol] = pow(a,b)
sol = a^b;
end


function [obj_f] = objective_function1(T, N)
obj_f = dot(T,T)/(N*N);
end


function [dcda_f] = dcda_fd(T, pctmetal, N, penal)
dcda_f = zeros(N,N);
q = 2/(0.01*0.01*0.001);

cost_original = objective_function1(T, N);

for i = 1:1:N
    for j = 1:1:N
        pctmetal(i,j) = pctmetal(i,j) - 0.01;
        
        [k] = createK(pctmetal, N, penal);
        [u] = fvm_func(k, N, q);
        
        
        cost_update = objective_function1(u, N);
        
        dcda_f(i,j) = (cost_original - cost_update)/0.01;
        pctmetal(i,j) = pctmetal(i,j) + 0.01;
    end
end

end


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


% inputs: N, pctmetal, 0.4, dcda_mat
function [xnew] = OC(N, x, volfrac, dc)
%OC Summary of this function goes here
%   Detailed explanation goes here
l1 = 0;
l2 = 100000;
move = 0.2;
while (l2-l1)>1e-4
    lmid = 0.5*(l2+l1);
    xnew = max(0.001, max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
    if sum(sum(xnew))-volfrac*N*N>0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end





% times = 4;
% N = times*N;
% dp = N+1;
% pctmetal = meshrefine(pctmetal,times);
%
% loop2 = loop1;
% change = 1.;
% while loop2<200%change > 0.04
%     loop2 = loop2 +1
%     pctmetal_old = pctmetal;
%     [k] = createK(pctmetal, N);
%     [T, K, f] = fvm_func(k, N, q);
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
% %     surface(Tmat)
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
%
%
%
% % times = 2;
% % N = times*N;
% % dp = N+1;
% % pctmetal = meshrefine(pctmetal,times);
% %
% % loop = 0;
% % change = 1.;
% % while change > 0.01
% %     loop = loop +1
% %     pctmetal_old = pctmetal;
% %     [T, K, f] = fvm_func(pctmetal, N, q);
% %
% % %     figure()
% % %     spy(K)
% % %
% % %     Tmat = zeros(dp,dp);
% % %     for i = 1:dp
% % %          Tmat(:,i)=T(i*dp-dp+1:i*dp);
% % %     end
% % %
% % %     figure()
% % %     surface(Tmat)
% %
% %     %cost = cost2(T, dp);
% %
% %     lambda_vec = lambda1(T, K, dp);
% %     dcda_mat1 = dcda(lambda_vec, T, pctmetal, N);
% %     dcda_mat = check(N, rmin, pctmetal, dcda_mat1);
% %
% %     pctmetal = OC(N, pctmetal, 0.4, dcda_mat);
% %
% %     change = max(max(abs(pctmetal-pctmetal_old)))
% % end
%
% figure()
% surface(pctmetal)
%
% times = 2;
% N = times*N;
% dp = N+1;
% pctmetal = meshrefine(pctmetal,times);
%
% loop3 = loop2;
% change = 1.;
% while loop3<300%change > 0.03
%     loop3 = loop3 +1
%     pctmetal_old = pctmetal;
%     [k] = createK(pctmetal, N);
%     [T, K, f] = fvm_func(k, N, q);
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
% %     surface(Tmat)
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
%
% figure()
% surface(pctmetal)
%
% times = 10/8;
% N = times*N;
% dp = N+1;
% pctmetal = meshrefine(pctmetal,times);
%
% loop4 = loop3;
% change = 1.;
% while change > 0.01
%     loop4 = loop4 +1
%     pctmetal_old = pctmetal;
%     [k] = createK(pctmetal, N);
%     [T, K, f] = fvm_func(k, N, q);
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
% %     surface(Tmat)
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
%
% figure()
% surface(pctmetal)
%
% %%
% heatsink = pctmetal*60 + (1-pctmetal)*0.2;
%
% [T, K, f] = fvm_func(heatsink, N, q);
%
% Tmat = zeros(dp,dp);
% for i = 1:dp
%      Tmat(:,i)=T(i*dp-dp+1:i*dp);
% end
%
% figure()
% surface(dcda_mat1, 'FaceColor', 'interp')
%
%
% figure()
% surface(Tmat)
%
% figure()
% surface(heatsink)
%
% figure()
% surface(pctmetal)
%
% % N = 100 -> 179 loops nodig




