clear all
close all


% gegevens
N = 50;
q = 2/(0.01*0.01*0.001);
rmin = 2;

dp = N+1;

% percentage metaal
pctmetal = 0.3*ones(N,N);
penal = 3.0;
dcda_mat = zeros(N,N);

loop1 = 0;
change = 1.0;
while loop1 <= 50 %change > 0.01
    
        
    pctmetal_old = pctmetal;
    
    [k] = createK(pctmetal, N, penal);
    [T, K, f] = fvm_func(k, N, q);
    
    lambda_vec = lambda1(T, K, N);
    %dcda_f = dcda_fd(T, pctmetal, N, penal);
    %dcda_a = dcda_arit(lambda_vec, T, pctmetal, N, penal);
    dcda_h = dcda_harm(lambda_vec, T, pctmetal, k, N, penal);
    
    for i = 1:1:N
           %dcda_mat(:,i) = dcda_a((i-1)*N+1:i*N);
           dcda_mat(:,i) = dcda_h((i-1)*N+1:i*N);
    end
    
    
    %dcda_check = check(N, rmin, pctmetal, dcda_f);
    dcda_check = check(N, rmin, pctmetal, dcda_mat);
   
    
    pctmetal = OC(N, pctmetal, 0.4, dcda_check);
    pctmetal_iter{loop1+1} = pctmetal;
    
    change = max(max(abs(pctmetal-pctmetal_old)));
    
    % Load C++ solution
    
    iter = int2str(loop1);
    
    m = strcat('metal_',iter);
    metal = strcat(m,'.txt');
    fileID = fopen(metal,'r');
    formatSpec = '%f';
    metal = fscanf(fileID,formatSpec);
    
    g = strcat('gradient_',iter);
    gradient = strcat(g,'.txt');
    fileID = fopen(gradient,'r');
    formatSpec = '%f';
    gradient = fscanf(fileID,formatSpec);

    metal_mat = zeros(N,N);
    gradient_mat = zeros(N,N);
    
    dp = N;
    for i = 1:dp
        metal_mat(:,i)=metal(i*dp-dp+1:i*dp); 
        gradient_mat(:,i)=gradient(i*dp-dp+1:i*dp);
    end
    
    difference_metal{loop1+1} = metal_mat-pctmetal;
    difference_gradient{loop1+1} = gradient_mat-dcda_check;
 
    loop1 = loop1 + 1;
end

for i = 1:1:loop1-1
    hfig = figure;
    pos = get(hfig,'position');
    set(hfig,'position',pos.*[.5 1 2 1]);
    subplot(1,2,1)
    surface(difference_gradient{i},'FaceColor','interp')
    colorbar
    subplot(1,2,2)
    surface(pctmetal_iter{i},'FaceColor','interp')
    colorbar
    pause(0.5)
end





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



function [lambda] = lambda1(T, K, N)
lambda = transpose(K)\((-2.0/N^2)*T);
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



function [dcda_m] = dcda_arit(lambda, T, a, N, penal)
dcda_m = sparse(N*N,1);
dcdk = sparse(N*N,1);

% element per element opvullen
for i = 2:1:N-1
    for j = 1:1:N-2
        dKdk_u = sparse((N+1)*(N+1),1);
        dKdk_u(j*(N+1)+i) = -1.0*T(j*(N+1)+i) + 0.5*T(j*(N+1)+i+1) + 0.5*T((j+1)*(N+1)+i);
        dKdk_u(j*(N+1)+i+1) = -1.0*T(j*(N+1)+i+1) + 0.5*T(j*(N+1)+i) + 0.5*T((j+1)*(N+1)+i+1);
        dKdk_u((j+1)*(N+1)+i) = -1.0*T((j+1)*(N+1)+i) + 0.5*T(j*(N+1)+i) + 0.5*T((j+1)*(N+1)+i+1);
        dKdk_u((j+1)*(N+1)+i+1) = -1.0*T((j+1)*(N+1)+i+1) + 0.5*T(j*(N+1)+i+1) + 0.5*T((j+1)*(N+1)+i);
        
        dcdk(i + (j)*N) = dot(lambda, dKdk_u);
        dcda_m(i + (j)*N) = penal*(65.0-0.2)*pow(a(i + N*(j)),penal-1)*dcdk(i + (j)*N);
    end
end

%Linksboven
dKdk_ulb = sparse((N+1)*(N+1),1);
dKdk_ulb(1) = -1.0*T(1) + 0.5*T(2) + 0.5*T(N+2);
dKdk_ulb(2) = -1.0*T(2) + 0.5*T(1) + 0.5*T(N+3);
dKdk_ulb(N+2) = -1.0*T(N+2) + 0.5*T(1) + 0.5*T(N+3);
dKdk_ulb(N+3) = -1.0*T(N+3) + 0.5*T(2) + 0.5*T(N+2);

dcdk(1) = dot(lambda, dKdk_ulb);
dcda_m(1) = penal*(65.0-0.2)*pow(a(1),penal-1)*dcdk(1);

%Linksonder
dKdk_ulo = sparse((N+1)*(N+1),1);
dKdk_ulo(N+1) = -1.0*T(N+1) + 0.5*T(N) + 0.5*T(2*N+2);
dKdk_ulo(N) = -1.0*T(N) + 0.5*T(N+1) + 0.5*T(2*N+1);
dKdk_ulo(2*N+2) = -1.0*T(2*N+2) + 0.5*T(N+1) + 0.5*T(2*N+1);
dKdk_ulo(2*N+1) = -1.0*T(2*N+1) + 0.5*T(N) + 0.5*T(2*N+2);

dcdk(N) = dot(lambda, dKdk_ulo);
dcda_m(N) = penal*(65.0-0.2)*pow(a(N),penal-1)*dcdk(N);

%Rechtsboven
dKdk_urb = sparse((N+1)*(N+1),1);
dKdk_urb(N*N+N+1) = -1.0*T(N*N+N+1) + 0.5*T(N*N+N+2) + 0.5*T(N*N);
dKdk_urb(N*N) = -1.0*T(N*N) + 0.5*T(N*N+N+1) + 0.5*T(N*N+1);
dKdk_urb(N*N+N+2) = -1.0*T(N*N+N+2) + 0.5*T(N*N+N+1) + 0.5*T(N*N+1);
dKdk_urb(N*N+1) = -1.0*T(N*N+1) + 0.5*T(N*N) + 0.5*T(N*N+N+2);

dcdk((N-1)*N+1) = dot(lambda, dKdk_urb);
dcda_m((N-1)*N+1) = penal*(65.0-0.2)*pow(a((N-1)*N+1),penal-1)*dcdk((N-1)*N+1);

%Rechtsonder
dKdk_uro = sparse((N+1)*(N+1),1);
dKdk_uro(N*N + 2*N+1) = -1.0*T(N*N + 2*N + 1) + 0.5*T(N*N + 2*N ) + 0.5*T(N*N + N);
dKdk_uro(N*N + 2*N) = -1.0*T(N*N + 2*N) + 0.5*T(N*N + 2*N+1) + 0.5*T(N*N + N - 1);
dKdk_uro(N*N + N) = -1.0*T(N*N + N ) + 0.5*T(N*N + 2*N+1) + 0.5*T(N*N + N - 1);
dKdk_uro(N*N + N - 1) = -1.0*T(N*N + N - 1) + 0.5*T(N*N + 2*N) + 0.5*T(N*N + N);

dcdk(N*N) = dot(lambda, dKdk_uro);
dcda_m(N*N) = penal*(65.0-0.2)*pow(a(N*N),penal-1)*dcdk(N*N);


% Boven en onderrand

for j = 1:1:N-2
    dKdk_ub = sparse((N+1)*(N+1),1);
    dKdk_ub(j*(N+1)+1) = -1*T(j*(N+1)+1) + 0.5*T(j*(N+1)+2) + 0.5*T((j+1)*(N+1)+1);
    dKdk_ub(j*(N+1)+2) = -1*T(j*(N+1)+2) + 0.5*T(j*(N+1)+1) + 0.5*T((j+1)*(N+1)+2);
    dKdk_ub((j+1)*(N+1)+1) = -1*T((j+1)*(N+1)+1) + 0.5*T((j+1)*(N+1)+2) + 0.5*T(j*(N+1)+1);
    dKdk_ub((j+1)*(N+1)+2) = -1*T((j+1)*(N+1)+2) + 0.5*T((j+1)*(N+1)+1) + 0.5*T(j*(N+1)+2);
    
    dcdk(j*N+1) = dot(lambda, dKdk_ub);
    dcda_m(j*N+1) = penal*(65.0-0.2)*pow(a(j*N+1),penal-1)*dcdk(j*N+1);
    
    dKdk_uo = sparse((N+1)*(N+1),1);
    dKdk_uo((j+1)*(N+1)-1) = -1*T((j+1)*(N+1)-1) + 0.5*T((j+1)*(N+1)) + 0.5*T((j+2)*(N+1)-1);
    dKdk_uo((j+1)*(N+1)) = -1*T((j+1)*(N+1)) + 0.5*T((j+1)*(N+1)-1) + 0.5*T((j+2)*(N+1));
    dKdk_uo((j+2)*(N+1)-1) = -1*T((j+2)*(N+1)-1) + 0.5*T((j+1)*(N+1)-1) + 0.5*T((j+2)*(N+1));
    dKdk_uo((j+2)*(N+1)) = -1*T((j+2)*(N+1)) + 0.5*T((j+1)*(N+1)) + 0.5*T((j+2)*(N+1)-1);
    
    dcdk((j+1)*N) = dot(lambda, dKdk_ub);
    dcda_m((j+1)*N) = penal*(65.0-0.2)*pow(a((j+1)*N),penal-1)*dcdk((j+1)*N);
end

%Neumann links
for i = 2:1:(0.3*N-1)
    dKdk_unlb = sparse((N+1)*(N+1),1);
    dKdk_unlb(i) = -1*T(i) + 0.5*T(i+1) + 0.5*T(i + (N+1));
    dKdk_unlb(i+1) = -1*T(i+1) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
    dKdk_unlb(i + (N+1)) = -1*T(i + (N+1)) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
    dKdk_unlb(i+1 + (N+1)) = -1*T(i+1 + (N+1)) + 0.5*T(i+1) + 0.5*T(i + (N+1));
    
    dcdk(i) = dot(lambda, dKdk_unlb);
    dcda_m(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
    
end


for i = (0.7*N+2):1:(N-1)
    dKdk_unlo = sparse((N+1)*(N+1),1);
    dKdk_unlo(i) = -1*T(i) + 0.5*T(i+1) + 0.5*T(i + (N+1));
    dKdk_unlo(i+1) = -1*T(i+1) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
    dKdk_unlo(i + (N+1)) = -1*T(i + (N+1)) + 0.5*T(i) + 0.5*T(i+1 + (N+1));
    dKdk_unlo(i+1 + (N+1)) = -1*T(i+1 + (N+1)) + 0.5*T(i+1) + 0.5*T(i + (N+1));
    
    dcdk(i) = dot(lambda, dKdk_unlo);
	dcda_m(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
    
end



%Neumann rechts
for i = 2:1:(0.3*N-1)
    dKdk_unrb = sparse((N+1)*(N+1),1);
    dKdk_unrb(i + (N-1)*(N+1)) = -1.0*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));
    dKdk_unrb(i+1 + (N-1)*(N+1)) = -1.0*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N)*(N+1));
    dKdk_unrb(i + (N)*(N+1)) = -1.0*T(i + (N)*(N+1)) + 0.5*T(i+1 + (N)*(N+1)) + 0.5*T(i + (N-1)*(N+1));
    dKdk_unrb(i+1 + (N)*(N+1)) = -1.0*T(i+1 + (N)*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));
   
    dcdk(N*(N-1)+i) = dot(lambda, dKdk_unrb);
    dcda_m(N*(N-1)+i) = penal*(65.0-0.2)*pow(a(N*(N-1)+i),penal-1)*dcdk(N*(N-1)+i);
end

for i = (0.7*N+2):1:N-1
    dKdk_unlo = sparse((N+1)*(N+1),1);
    dKdk_unlo(i + (N-1)*(N+1)) = -1.0*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));
    dKdk_unlo(i+1 + (N-1)*(N+1)) = -1.0*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N)*(N+1));
    dKdk_unlo(i + (N)*(N+1)) = -1.0*T(i + (N)*(N+1)) + 0.5*T(i+1 + (N)*(N+1)) + 0.5*T(i + (N-1)*(N+1));
    dKdk_unlo(i+1 + (N)*(N+1)) = -1.0*T(i+1 + (N)*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));
   
    dcdk(N*(N-1)+i) = dot(lambda, dKdk_unlo);
	dcda_m(N*(N-1)+i) = penal*(65.0-0.2)*pow(a(N*(N-1)+i),penal-1)*dcdk(N*(N-1)+i);
end


%Vakje met 1 hoekpunt van 293K, linkerkant, boven dirichlett boundary condition
dKdk_uhlb = sparse((N+1)*(N+1),1);
dKdk_uhlb(0.3*N) = -1*T(0.3*N) + 0.5*T(0.3*N+1) + 0.5*T(0.3*N +(N+1));
%(0.3*N+1)
dKdk_uhlb(0.3*N +(N+1)) = -1*T(0.3*N +(N+1)) + 0.5*T(0.3*N) + 0.5*T(0.3*N +(N+1)+ 1);
dKdk_uhlb(0.3*N +(N+1)+ 1) = -1*T(0.3*N +(N+1)+ 1) + 0.5*T(0.3*N+1) + 0.5*T(0.3*N +(N+1));

dcdk(0.3*N) = dot(lambda, dKdk_uhlb);
dcda_m(0.3*N) = penal*(65.0-0.2)*pow(a(0.3*N),penal-1)*dcdk(0.3*N);

%Vakje met 1 hoekpunt van 293K, linkerkant, onder dirichlett boundary condition
dKdk_uhlo = sparse((N+1)*(N+1),1);
dKdk_uhlo(0.7*N+2) = -1*T(0.7*N+2) + 0.5*T(0.7*N +(N+1)+2) + 0.5*T(0.7*N+1);
%(0.7*N+1)
dKdk_uhlo(0.7*N+1+(N+1)) = -1*T(0.7*N+1+(N+1)) + 0.5*T(0.7*N +(N+1)+2) + 0.5*T(0.7*N+1);
dKdk_uhlo(0.7*N +(N+1)+2) = -1*T(0.7*N +(N+1)+2) + 0.5*T(0.7*N+1+(N+1)) + 0.5*T(0.7*N+2);

dcdk(0.7*N+1) = dot(lambda, dKdk_uhlo);
dcda_m(0.7*N+1) = penal*(65.0-0.2)*pow(a(0.7*N+1),penal-1)*dcdk(0.7*N+1);

%Vakjes direct naast de linkse dirichlett boundary condition. Deze hebben 2 hoekpunten van 293K
for i = 0.3*N+1:1:0.7*N
    dKdk_udl = sparse((N+1)*(N+1),1);
    dKdk_udl(i+(N+1)) = -1*T(i+(N+1))+0.5*T(i+(N+1)+1)+0.5*T(i);
    dKdk_udl(i+(N+1)+1) = -1*T(i+(N+1)+1)+0.5*T(i+(N+1))+0.5*T(i+1);
    %(i)
    %(i+1)
    dcdk(i) = dot(lambda, dKdk_udl);
    dcda_m(i) =  penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
    
end

%Vakje met 1 hoekpunt van 293K, rechterkant, boven dirichlett boundary condition
dKdk_uhrb = sparse((N+1)*(N+1),1);
dKdk_uhrb((N+1)*N + 0.3*N) = -1*T((N+1)*N + 0.3*N) + 0.5*T((N+1)*N + 0.3*N+1) + 0.5*T((N+1)*(N-1) + 0.3*N);
%(N+1)*N + 0.3*N+1
dKdk_uhrb((N+1)*(N-1) + 0.3*N) = -1*T((N+1)*(N-1) + 0.3*N) + 0.5*T((N+1)*N + 0.3*N) + 0.5*T((N+1)*(N-1) + 0.3*N+1);
dKdk_uhrb((N+1)*(N-1) + 0.3*N+1) = -1*T((N+1)*(N-1) + 0.3*N+1) + 0.5*T((N+1)*N + 0.3*N+1) + 0.5*T((N+1)*(N-1) + 0.3*N);

dcdk((N-1)*N + 0.3*N) = dot(lambda, dKdk_uhrb);
dcda_m((N-1)*N + 0.3*N) =  penal*(65.0-0.2)*pow(a((N-1)*N + 0.3*N),penal-1)*dcdk((N-1)*N + 0.3*N);

%Vakje met 1 hoekpunt van 293K, rechterkant, onder dirichlett boundary condition
dKdk_uhro = sparse((N+1)*(N+1),1);
%(N*(N+1) + 0.7*N)
dKdk_uhro(0.7*N + (N-1)*(N+1)+1) = -1*T(0.7*N + (N-1)*(N+1)+1) + 0.5*T(0.7*N+1 + (N-1)*(N+1)+1) + 0.5*T(0.7*N + (N)*(N+1));
dKdk_uhro(0.7*N+1 + (N-1)*(N+1)+1) = -1*T(0.7*N+1 + (N-1)*(N+1)+1) + 0.5*T(0.7*N + (N-1)*(N+1)+1) + 0.5*T(0.7*N+1 + (N)*(N+1)+1);
dKdk_uhro(0.7*N+1 + (N)*(N+1)+1) = -1*T(0.7*N+1 + (N)*(N+1)+1) + 0.5*T(0.7*N+1 + (N-1)*(N+1)+1) + 0.5*T(0.7*N + (N)*(N+1));
%0.7*N + (N)*(N+1)
dcdk(N*(N-1) + 0.7*N +1) = dot(lambda, dKdk_uhro);
dcda_m(N*(N-1) + 0.7*N +1) = penal*(65.0-0.2)*pow(a(N*(N-1) + 0.7*N +1),penal-1)*dcdk(N*(N-1) + 0.7*N +1);



%Vakjes direct naast de rechtse dirichlett boundary condition. Deze hebben 2 hoekpunten van 293K

for i = 0.3*N+1:1:0.7*N
    dKdk_udr = sparse((N+1)*(N+1),1);
    dKdk_udr(i + (N-1)*(N+1)) = -1*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N-1)*(N+1))  + 0.5*T(i + (N)*(N+1));
    dKdk_udr(i+1 + (N-1)*(N+1)) = -1*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N)*(N+1));
    %(i + (N)*(N+1))
    %(i+1 + (N)*(N+1))
    
    dcdk(i + (N-1)*N) = dot(lambda, dKdk_udr);
    dcda_m(i + (N-1)*N) = penal*(65.0-0.2)*pow(a(i + (N-1)*N),penal-1)*dcdk(i + (N-1)*N);
end




end


function [dcda_m] = dcda_harm(lambda, T, a, k, N, penal)
dcda_m = sparse(N*N,1);
dcdk = sparse(N*N,1);
% element per element opvullen
for i = 2:1:N-1
    for j = 2:1:N-1
        dKdk_u = sparse((N+1)*(N+1),1);
        dKdk_u(i + (j-1)*(N+1)) = (-2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2) + pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2)))*T(i + (j-1)*(N+1)) + (2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2)))*T(i+1 + (j-1)*(N+1)) + (2.0*pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2))*T(i + (j)*(N+1));
        dKdk_u(i+1 + (j-1)*(N+1)) = (-2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2) + pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2)))*T(i+1 + (j-1)*(N+1)) + (2.0*(pow(k(i,j-1),2)/pow((k(i,j-1) + k(i,j)),2)))*T(i + (j-1)*(N+1)) + (2.0*pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2))*T(i+1 + (j)*(N+1));
        dKdk_u(i + (j)*(N+1)) = (-2.0*(pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2) + pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i + (j)*(N+1)) + (2.0*(pow(k(i-1,j),2)/pow((k(i-1,j) + k(i,j)),2)))*T(i + (j-1)*(N+1)) + (2.0*pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2))*T(i+1 + (j)*(N+1));
        dKdk_u(i+1 + (j)*(N+1)) = (-2.0*(pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2) + pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i+1 + (j)*(N+1)) + (2.0*(pow(k(i+1,j),2)/pow((k(i+1,j) + k(i,j)),2)))*T(i+1 + (j-1)*(N+1)) + (2.0*(pow(k(i,j+1),2)/pow((k(i,j+1) + k(i,j)),2)))*T(i + (j)*(N+1));
        			
        dcdk(i + (j-1)*N) = dot(lambda, dKdk_u);
        dcda_m(i + (j-1)*N) = penal*(65.0-0.2)*pow(a(i + N*(j-1)),penal-1)*dcdk(i + (j-1)*N);
    end
end

%Linksboven
dKdk_ulb = sparse((N+1)*(N+1),1);
dKdk_ulb(1) = -1.0*T(1) + 0.5*T(2) + 0.5*T(N+2);
dKdk_ulb(2) = (-0.5 - 2.0*(pow(k(2,1),2)/pow((k(2,1) + k(1,1)),2)))*T(2) + 0.5*T(1) + (2.0*pow(k(2,1),2)/pow((k(2,1) + k(1,1)),2))*T(N+3);
dKdk_ulb(N+2) = (-0.5 - 2.0*(pow(k(1,2),2)/pow((k(1,2) + k(1,1)),2)))*T(N+2) + 0.5*T(1) + (2.0*pow(k(1,2),2)/pow((k(1,2) + k(1,1)),2))*T(N+3);
dKdk_ulb(N+3) = (-2.0*(pow(k(2,1),2)/pow((k(2,1) + k(1,1)),2) + pow(k(1,2),2)/pow((k(1,2) + k(1,1)),2)))*T(N+3) + (2.0*(pow(k(2,1),2)/pow((k(2,1) + k(1,1)),2)))*T(2) + (2.0*(pow(k(1,2),2)/pow((k(1,2) + k(1,1)),2)))*T(N+2);

dcdk(1) = dot(lambda, dKdk_ulb);
dcda_m(1) = penal*(65.0-0.2)*pow(a(1),penal-1)*dcdk(1);

%Linksonder
dKdk_ulo = sparse((N+1)*(N+1),1);
dKdk_ulo(N) = (-0.5 - 2.0*(pow(k(N-1,1),2)/pow((k(N-1,1) + k(N,1)),2)))*T(N) + 0.5*T(N+1) + (2.0*pow(k(N-1,1),2)/pow((k(N-1,1) + k(N,1)),2))*T(2*N+1);
dKdk_ulo(N+1) = -1.0*T(N+1) + 0.5*T(N) + 0.5*T(2*N+2);
dKdk_ulo(2*N+1) = (-2.0*(pow(k(N-1,1),2)/pow((k(N-1,1) + k(N,1)),2) + pow(k(N,2),2)/pow((k(N,2) + k(N,1)),2)))*T(2*N+1) + (2.0*(pow(k(N-1,1),2)/pow((k(N-1,1) + k(N,1)),2)))*T(N) + (2.0*pow(k(N,2),2)/pow((k(N,2) + k(N,1)),2))*T(2*N+2);
dKdk_ulo(2*N+2) = (-0.5 - 2.0*(pow(k(N,2),2)/pow((k(N,2) + k(N,1)),2)))*T(2*N+2) + 0.5*T(N+1) + (2.0*(pow(k(N,2),2)/pow((k(N,2) + k(N,1)),2)))*T(2*N+1);

dcdk(N) = dot(lambda, dKdk_ulo);
dcda_m(N) = penal*(65.0-0.2)*pow(a(N),penal-1)*dcdk(N);

%Rechtsboven
dKdk_urb = sparse((N+1)*(N+1),1);
dKdk_urb(N*N) = (-2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(1,N)),2)) - 0.5)*T(N*N) + (2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(1,N)),2)))*T(N*N+1) + 0.5*T(N*(N+1)+1);
dKdk_urb(N*N+1) = (-2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(1,N)),2) + pow(k(2,N),2)/pow((k(2,N) + k(1,N)),2)))*T(N*N+1) + (2.0*(pow(k(1,N-1),2)/pow((k(1,N-1) + k(1,N)),2)))*T(N*N) + (2.0*pow(k(2,N),2)/pow((k(2,N) + k(1,N)),2))*T(N*(N+1) +2);
dKdk_urb(N*(N+1)+1) = -1.0*T(N*(N+1)+1) + 0.5*T(N*N) + 0.5*T(N*(N+1) +2);
dKdk_urb(N*(N+1) +2) = (-2.0*(pow(k(2,N),2)/pow((k(2,N) + k(1,N)),2)) - 0.5)*T(N*(N+1)+1 + 1) + (2.0*(pow(k(2,N),2)/pow((k(2,N) + k(1,N)),2)))*T(N*N+1) + 0.5*T(N*(N+1)+1);
	
dcdk((N-1)*N+1) = dot(lambda, dKdk_urb);
dcda_m((N-1)*N+1) = penal*(65.0-0.2)*pow(a((N-1)*N+1),penal-1)*dcdk((N-1)*N+1);

%Rechtsonder
dKdk_uro = sparse((N+1)*(N+1),1);
dKdk_uro(N*N + 2*N+1) = -1.0*T(N*N + 2*N+1) + 0.5*T(N*N + N) + 0.5*T(N*N + 2*N);
dKdk_uro(N*N + 2*N) = (-2.0*(pow(k(N-1,N),2)/pow((k(N-1,N) + k(N,N)),2)) - 0.5)*T(N*N + 2*N) + (2.0*(pow(k(N-1,N),2)/pow((k(N-1,N) + k(N,N)),2)))*T(N*N + N - 1) + 0.5*T(N*N + 2*N+1);
dKdk_uro(N*N + N) = (-2.0*(pow(k(N,N-1),2)/pow((k(N,N-1) + k(N,N)),2)) - 0.5)*T(N*N + N) + (2.0*(pow(k(N,N-1),2)/pow((k(N,N-1) + k(N,N)),2)))*T(N*N + N - 1) + 0.5*T(N*N + 2*N +1);
dKdk_uro(N*N + N - 1) = (-2.0*(pow(k(N,N-1),2)/pow((k(N,N-1) + k(N,N)),2) + pow(k(N-1,N),2)/pow((k(N-1,N) + k(N,N)),2)))*T(N*N + N - 1) + (2.0*(pow(k(N,N-1),2)/pow((k(N,N-1) + k(N,N)),2)))*T(N*N + N) + (2.0*pow(k(N-1,N),2)/pow((k(N-1,N) + k(N,N)),2))*T(N*N + 2*N);

dcdk(N*N) = dot(lambda, dKdk_uro);
dcda_m(N*N) = penal*(65.0-0.2)*pow(a(N*N),penal-1)*dcdk(N*N);


% Boven en onderrand

for j = 1:1:N-2
    dKdk_ub = sparse((N+1)*(N+1),1);
    dKdk_ub(j*(N+1)+1) = (-2.0*(pow(k(1,j),2)/pow((k(1,j) + k(1,j+1)),2)) - 0.5)*T(j*(N+1)+1) + (2.0*(pow(k(1,j),2)/pow((k(1,j) + k(1,j+1)),2)))*T(2 + j*(N+1)) + 0.5*T((j+1)*(N+1)+1);
    dKdk_ub(j*(N+1)+2) = (-2.0*(pow(k(1,j),2)/pow((k(1,j) + k(1,j+1)),2) + pow(k(2,j+1),2)/pow((k(2,j+1) + k(1,j+1)),2)))*T(2 + j*(N+1)) + (2.0*(pow(k(1,j),2)/pow((k(1,j) + k(1,j+1)),2)))*T(j*(N+1)+1) + (2.0*pow(k(2,j+1),2)/pow((k(2,j+1) + k(1,j+1)),2))*T(2 + (j+1)*(N+1));
    dKdk_ub((j+1)*(N+1)+1) = (-0.5 - 2.0*(pow(k(1,j+2),2)/pow((k(1,j+2) + k(1,j+1)),2)))*T((j+1)*(N+1)+1) + 0.5*T(j*(N+1)+1) + (2.0*pow(k(1,j+2),2)/pow((k(1,j+2) + k(1,j+1)),2))*T(2 + (j+1)*(N+1));
    dKdk_ub((j+1)*(N+1)+2) = (-2.0*(pow(k(2,j+1),2)/pow((k(2,j+1) + k(1,j+1)),2) + pow(k(1,j+2),2)/pow((k(1,j+2) + k(1,j+1)),2)))*T(2 + (j+1)*(N+1)) + (2.0*(pow(k(2,j+1),2)/pow((k(2,j+1) + k(1,j+1)),2)))*T(2 + j*(N+1)) + (2.0*(pow(k(1,j+2),2)/pow((k(1,j+2) + k(1,j+1)),2)))*T((j+1)*(N+1)+1);
    
    dcdk(j*N+1) = dot(lambda, dKdk_ub);
    dcda_m(j*N+1) = penal*(65.0-0.2)*pow(a(j*N+1),penal-1)*dcdk(j*N+1);
    
    dKdk_uo = sparse((N+1)*(N+1),1);
    dKdk_uo((j+1)*(N+1)-1) = (-2.0*(pow(k(N,j),2)/pow((k(N,j) + k(N,j+1)),2) + pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N,j+1)),2)))*T(N + j*(N+1)) + (2.0*(pow(k(N,j),2)/pow((k(N,j) + k(N,j+1)),2)))*T(N +1+ j*(N+1)) + (2.0*pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N,j+1)),2))*T(N + (j+1)*(N+1));
    dKdk_uo((j+1)*(N+1)) = (-2.0*(pow(k(N,j),2)/pow((k(N,j) + k(N,j+1)),2)) - 0.5)*T(N+1 + j*(N+1)) + (2.0*(pow(k(N,j),2)/pow((k(N,j) + k(N,j+1)),2)))*T(N + j*(N+1)) + 0.5*T(N + 1+(j+1)*(N+1));
    dKdk_uo((j+2)*(N+1)-1) = (-2.0*(pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N,j+1)),2) + pow(k(N,j+2),2)/pow((k(N,j+2) + k(N,j+1)),2)))*T(N + (j+1)*(N+1)) + (2.0*(pow(k(N-1,j+1),2)/pow((k(N-1,j+1) + k(N,j+1)),2)))*T(N + j*(N+1)) + (2.0*pow(k(N,j+2),2)/pow((k(N,j+2) + k(N,j+1)),2))*T(N + 1+(j+1)*(N+1));
    dKdk_uo((j+2)*(N+1)) = (-0.5 - 2.0*(pow(k(N,j+2),2)/pow((k(N,j+2) + k(N,j+1)),2)))*T(N+1 + (j+1)*(N+1)) + 0.5*T(N+1 + j*(N+1)) + (2.0*(pow(k(N,j+2),2)/pow((k(N,j+2) + k(N,j+2)),2)))*T(N + (j+1)*(N+1));
    
    dcdk((j+1)*N) = dot(lambda, dKdk_ub);
    dcda_m((j+1)*N) = penal*(65.0-0.2)*pow(a((j+1)*N),penal-1)*dcdk((j+1)*N);
end

%Neumann links boven en onder dirichlett
for i = 2:1:(0.3*N-1)
    dKdk_unlb = sparse((N+1)*(N+1),1);
    dKdk_unlb(i) = (-0.5 - 2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2)))*T(i) + 0.5*T(i+1) + (2.0*pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2))*T(i + (N+1));
    dKdk_unlb(i+1) = (-0.5 - 2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2)))*T(i+1) + 0.5*T(i) + (2.0*pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2))*T(i+1 + (N+1));
	dKdk_unlb(i + (N+1)) = (-2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2)))*T(i) + (2.0*pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2))*T(i+1 + (N+1));
	dKdk_unlb(i+1 + (N+1)) = (-2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2)))*T(i+1) + (2.0*(pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1));

    dcdk(i) = dot(lambda, dKdk_unlb);
    dcda_m(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
end

for i = (0.7*N+2):1:(N-1)
    dKdk_unlo = sparse((N+1)*(N+1),1);
    dKdk_unlo(i) = (-0.5 - 2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2)))*T(i) + 0.5*T(i+1) + (2.0*pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2))*T(i + (N+1));
    dKdk_unlo(i+1) = (-0.5 - 2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2)))*T(i+1) + 0.5*T(i) + (2.0*pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2))*T(i+1 + (N+1));
	dKdk_unlo(i + (N+1)) = (-2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2)))*T(i) + (2.0*pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2))*T(i+1 + (N+1));
	dKdk_unlo(i+1 + (N+1)) = (-2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2)))*T(i+1) + (2.0*(pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1));

    dcdk(i) = dot(lambda, dKdk_unlo);
	dcda_m(i) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
    
end


%Neumann rechts
for i = 2:1:(0.3*N-1)
    dKdk_unlb = sparse((N+1)*(N+1),1);
    
    dKdk_unlb(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2))*T(i + (N)*(N+1));
    dKdk_unlb(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2))*T(i+1 + (N)*(N+1));
    dKdk_unlb(i + (N)*(N+1)) = (-2.0*(pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)) - 0.5)*T(i + (N)*(N+1)) + (2.0*(pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N)*(N+1));
    dKdk_unlb(i+1 + (N)*(N+1)) = (-2.0*(pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)) - 0.5)*T(i+1 + (N)*(N+1)) + (2.0*(pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));

    dcdk(N*(N-1)+i) = dot(lambda, dKdk_unlb);
    dcda_m(N*(N-1)+i) = penal*(65.0-0.2)*pow(a(N*(N-1)+i),penal-1)*dcdk(N*(N-1)+i);
    
end

for i = (0.7*N+2):1:(N-1)
    
    dKdk_unlo = sparse((N+1)*(N+1),1);
   
    dKdk_unlo(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2))*T(i + (N)*(N+1));	
    dKdk_unlo(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2))*T(i+1 + (N)*(N+1));	
    dKdk_unlo(i + (N)*(N+1)) = (-2.0*(pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)) - 0.5)*T(i + (N)*(N+1)) + (2.0*(pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + 0.5*T(i+1 + (N)*(N+1));
	dKdk_unlo(i+1 + (N)*(N+1)) = (-2.0*(pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)) - 0.5)*T(i+1 + (N)*(N+1)) + (2.0*(pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + 0.5*T(i + (N)*(N+1));

    dcdk(N*(N-1)+i) = dot(lambda, dKdk_unlo);
	dcda_m(N*(N-1)+i) = penal*(65.0-0.2)*pow(a(N*(N-1)+i),penal-1)*dcdk(N*(N-1)+i);
end

%Links dirichlet
for i = (0.3*N+1):1:0.7*N
    dKdk_udl = sparse((N+1)*(N+1),1);
    dKdk_udl(i+(N+1)) = (-2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1)) + (2.0*(pow(k(i-1,1),2)/pow((k(i-1,1) + k(i,1)),2)))*T(i) + (2.0*pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2))*T(i+1 + (N+1));
    dKdk_udl(i+(N+1)+1) = (-2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2) + pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i+1 + (N+1)) + (2.0*(pow(k(i+1,1),2)/pow((k(i+1,1) + k(i,1)),2)))*T(i+1) + (2.0*(pow(k(i,2),2)/pow((k(i,2) + k(i,1)),2)))*T(i + (N+1));
    %(i)
    %(i+1)
    dcdk(i) = dot(lambda, dKdk_udl);
    dcda_m(i) =  penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i);
    
end

%rechts dirichlet
for i = (0.3*N+1):1:0.7*N
    dKdk_udr = sparse((N+1)*(N+1),1);
    
    dKdk_udr(i + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*pow(k(i-1,N),2)/pow((k(i-1,N) + k(i,N)),2))*T(i + (N)*(N+1));
	dKdk_udr(i+1 + (N-1)*(N+1)) = (-2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2) + pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2)))*T(i+1 + (N-1)*(N+1)) + (2.0*(pow(k(i,N-1),2)/pow((k(i,N-1) + k(i,N)),2)))*T(i + (N-1)*(N+1)) + (2.0*pow(k(i+1,N),2)/pow((k(i+1,N) + k(i,N)),2))*T(i+1 + (N)*(N+1));

    dcdk(i + (N-1)*N) = dot(lambda, dKdk_udr);
    dcda_m(i + (N-1)*N) = penal*(65.0-0.2)*pow(a(i),penal-1)*dcdk(i + (N-1)*N);
    
end


%Vakje met 1 hoekpunt van 293K, linkerkant, boven dirichlett boundary condition
dKdk_uhlb = sparse((N+1)*(N+1),1);
dKdk_uhlb(0.3*N) = (-0.5 - 2.0*(pow(k(0.3*N - 1,1),2)/pow((k(0.3*N - 1,1) + k(0.3*N,1)),2)))*T(0.3*N) + 0.5*T(0.3*N+1) + (2.0*pow(k(0.3*N - 1,1),2)/pow((k(0.3*N - 1,1) + k(0.3*N,1)),2))*T(0.3*N  + (N+1));
%(0.3*N+1)
dKdk_uhlb(0.3*N +(N+1)) = (-2.0*(pow(k(0.3*N - 1,1),2)/pow((k(0.3*N - 1,1) + k(0.3*N ,1)),2) + pow(k(0.3*N ,2),2)/pow((k(0.3*N,2) + k(0.3*N,1)),2)))*T(0.3*N + (N+1)) + (2.0*(pow(k(0.3*N - 1,1),2)/pow((k(0.3*N - 1,1) + k(0.3*N,1)),2)))*T(0.3*N) + (2.0*pow(k(0.3*N ,2),2)/pow((k(0.3*N ,2) + k(0.3*N ,1)),2))*T(0.3*N +1+ (N+1));
dKdk_uhlb(0.3*N +(N+1)+ 1) = (-2.0*(pow(k(0.3*N+1,1),2)/pow((k(0.3*N+1,1) + k(0.3*N,1)),2) + pow(k(0.3*N,2),2)/pow((k(0.3*N,2) + k(0.3*N,1)),2)))*T(0.3*N + 1+(N+1)) + (2.0*(pow(k(0.3*N+1,1),2)/pow((k(0.3*N+1,1) + k(0.3*N,1)),2)))*T(0.3*N+1) + (2.0*(pow(k(0.3*N,2),2)/pow((k(0.3*N,2) + k(0.3*N ,1)),2)))*T(0.3*N + (N+1));
dcdk(0.3*N) = dot(lambda, dKdk_uhlb);
dcda_m(0.3*N) = penal*(65.0-0.2)*pow(a(0.3*N),penal-1)*dcdk(0.3*N);

%Vakje met 1 hoekpunt van 293K, linkerkant, onder dirichlett boundary condition
dKdk_uhlo = sparse((N+1)*(N+1),1);
dKdk_uhlo(0.7*N+2) = (-0.5 - 2.0*(pow(k(0.7*N+2,1),2)/pow((k(0.7*N+2,1) + k(0.7*N+1,1)),2)))*T(0.7*N+2) + 0.5*T(0.7*N+1) + (2.0*pow(k(0.7*N+2,1),2)/pow((k(0.7*N+2,1) + k(0.7*N+1,1)),2))*T(0.7*N+2 + (N+1));
%(0.7*N+1)
dKdk_uhlo(0.7*N+1+(N+1)) = (-2.0*(pow(k(0.7*N,1),2)/pow((k(0.7*N,1) + k(0.7*N+1,1)),2) + pow(k(0.7*N+1,2),2)/pow((k(0.7*N+1,2) + k(0.7*N+1,1)),2)))*T(0.7*N + 1+(N+1)) + (2.0*(pow(k(0.7*N,1),2)/pow((k(0.7*N,1) + k(0.7*N+1,1)),2)))*T(0.7*N+1) + (2.0*pow(k(0.7*N+1,2),2)/pow((k(0.7*N+1,2) + k(0.7*N+1,1)),2))*T(0.7*N+2 + (N+1));
dKdk_uhlo(0.7*N +(N+1)+2) = (-2.0*(pow(k(0.7*N+2,1),2)/pow((k(0.7*N+2,1) + k(0.7*N+1,1)),2) + pow(k(0.7*N+1,2),2)/pow((k(0.7*N+1,2) + k(0.7*N+1,1)),2)))*T(0.7*N+2 + (N+1)) + (2.0*(pow(k(0.7*N+2,1),2)/pow((k(0.7*N+2,1) + k(0.7*N+1,1)),2)))*T(0.7*N+2) + (2.0*(pow(k(0.7*N+1,2),2)/pow((k(0.7*N+1,2) + k(0.7*N+1,1)),2)))*T(0.7*N + 1+(N+1));

dcdk(0.7*N+1) = dot(lambda, dKdk_uhlo);
dcda_m(0.7*N+1) = penal*(65.0-0.2)*pow(a(0.7*N+1),penal-1)*dcdk(0.7*N+1);


%Vakje met 1 hoekpunt van 293K, rechterkant, boven dirichlett boundary condition
dKdk_uhrb = sparse((N+1)*(N+1),1);


dKdk_uhrb(0.3*N  + (N-1)*(N+1)) = (-2.0*(pow(k(0.3*N ,N-1),2)/pow((k(0.3*N ,N-1) + k(0.3*N,N)),2) + pow(k(0.3*N - 1,N),2)/pow((k(0.3*N - 1,N) + k(0.3*N ,N)),2)))*T(0.3*N + (N-1)*(N+1)) + (2.0*(pow(k(0.3*N ,N-1),2)/pow((k(0.3*N,N-1) + k(0.3*N ,N)),2)))*T(0.3*N +1+(N-1)*(N+1)) + (2.0*pow(k(0.3*N - 1,N),2)/pow((k(0.3*N - 1,N) + k(0.3*N ,N)),2))*T(0.3*N  + (N)*(N+1));
dKdk_uhrb(0.3*N+1 + (N-1)*(N+1)) = (-2.0*(pow(k(0.3*N ,N-1),2)/pow((k(0.3*N ,N-1) + k(0.3*N,N)),2) + pow(k(0.3*N+1,N),2)/pow((k(0.3*N+1,N) + k(0.3*N ,N)),2)))*T(0.3*N + 1+(N-1)*(N+1)) + (2.0*(pow(k(0.3*N ,N-1),2)/pow((k(0.3*N ,N-1) + k(0.3*N ,N)),2)))*T(0.3*N + (N-1)*(N+1)) + (2.0*pow(k(0.3*N+1,N),2)/pow((k(0.3*N+1,N) + k(0.3*N,N)),2))*T(0.3*N + 1+(N)*(N+1));
dKdk_uhrb(0.3*N + (N)*(N+1)) = (-2.0*(pow(k(0.3*N - 1,N),2)/pow((k(0.3*N - 1,N) + k(0.3*N,N)),2)) - 0.5)*T(0.3*N  + (N)*(N+1)) + (2.0*(pow(k(0.3*N - 1,N),2)/pow((k(0.3*N - 1,N) + k(0.3*N ,N)),2)))*T(0.3*N + (N-1)*(N+1)) + 0.5*T(0.3*N +1+ (N)*(N+1));
	

dcdk((N-1)*N + 0.3*N) = dot(lambda, dKdk_uhrb);
dcda_m((N-1)*N + 0.3*N) =  penal*(65.0-0.2)*pow(a((N-1)*N + 0.3*N),penal-1)*dcdk((N-1)*N + 0.3*N);

%Vakje met 1 hoekpunt van 293K, rechterkant, onder dirichlett boundary condition
dKdk_hro = sparse((N+1)*(N+1),1);
%(N*(N+1) + 0.7*N)

dKdk_hro(0.7*N + (N-1)*(N+1)+1) = (-2.0*(pow(k(0.7*N+1,N-1),2)/pow((k(0.7*N+1,N-1) + k(0.7*N+1,N)),2) + pow(k(0.7*N,N),2)/pow((k(0.7*N,N) + k(0.7*N+1,N)),2)))*T(0.7*N + (N-1)*(N+1)+1) + (2.0*(pow(k(0.7*N+1,N-1),2)/pow((k(0.7*N+1,N-1) + k(0.7*N+1,N)),2)))*T(0.7*N+1 + (N-1)*(N+1)+1) + (2.0*pow(k(0.7*N+1-1,N),2)/pow((k(0.7*N,N) + k(0.7*N+1,N)),2))*T(0.7*N + (N)*(N+1)+1);
dKdk_hro(0.7*N+1 + (N-1)*(N+1)+1) = (-2.0*(pow(k(0.7*N+1,N-1),2)/pow((k(0.7*N+1,N-1) + k(0.7*N+1,N)),2) + pow(k(0.7*N+2,N),2)/pow((k(0.7*N+2,N) + k(0.7*N+1,N)),2)))*T(0.7*N+1 +1+ (N-1)*(N+1)) + (2.0*(pow(k(0.7*N+1,N-1),2)/pow((k(0.7*N+1,N-1) + k(0.7*N+1,N)),2)))*T(0.7*N +1+ (N-1)*(N+1)) + (2.0*pow(k(0.7*N+2,N),2)/pow((k(0.7*N+2,N) + k(0.7*N+1,N)),2))*T(0.7*N+2 + (N)*(N+1));
dKdk_hro(0.7*N+1 + (N)*(N+1)+1) = (-2.0*(pow(k(0.7*N+2,N),2)/pow((k(0.7*N+2,N) + k(0.7*N+1,N)),2)) - 0.5)*T(0.7*N+1 +1+ (N)*(N+1)) + (2.0*(pow(k(0.7*N+2,N),2)/pow((k(0.7*N+2,N) + k(0.7*N+1,N)),2)))*T(0.7*N+1 +1+ (N-1)*(N+1)) + 0.5*T(0.7*N + 1+(N)*(N+1));

dcdk(0.7*N + (N-1)*N+1) = dot(lambda, dKdk_hro);
dcda_m(0.7*N + (N-1)*N+1) = penal*(65.0-0.2)*pow(a(0.7*N + N*(N-1)+1),penal-1)*dcdk(0.7*N + (N-1)*N+1);

end









