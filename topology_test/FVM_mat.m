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
    h*(i-1);
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



ll = zeros(dp^2,dp^2);
N = dp-1;

%dirichlet
for i = (0.3*N+1):1:(0.7*N+1)
    ll(i,i) = 1;
    ll(dp^2-i+1, dp^2-i+1) = 1;
end

%hoekpunten
%linksboven
ll(1,1) = -0.5*(k(1,1)+k(1,1));
ll(1,2) = 0.5*k(1,1);
ll(1,dp+1) = 0.5*k(1,1);
%linksonder
ll(dp,dp) = -0.5*(k(N,1)+k(N,1));
ll(dp,dp-1) = 0.5*k(N,1);
ll(dp,2*dp) = 0.5*k(N,1);
%rechtsboven
ll((dp-1)*dp+1,(dp-1)*dp+1) = -0.5*(k(1,N)+k(1,N));
ll((dp-1)*dp+1,(dp-1)*dp+2) = 0.5*k(1,N);
ll((dp-1)*dp+1,dp*(dp-2)+1) = 0.5*k(1,N);
%rechtsonder
ll(dp^2,dp^2) = -0.5*(k(N,N)+k(N,N));
ll(dp^2,(dp-1)*dp) = 0.5*k(N-1,N-1);
ll(dp^2,dp^2-1) = 0.5*k(N-1,N-1);

%Neumann links&rechts
i = 1;
while i<0.3*N
    ll(i+1,i+1) = -0.5*k(i,1) -0.5*k(i+1,1) - ((k(i,1) + k(i+1,1))/2.0);
    ll(i+1,i) = 0.5*k(i,1);
    ll(i+1,i+1+1) = 0.5*k(i+1,1);
    ll(i+1,i+N+1+1) = (k(i,1) + k(i+1,1))/2.0;

    ll(N-i+1,N-i+1) = -0.5*k(N-i+1,1) -0.5*k(N-i,1) - ((k(N-i+1,1) + k(N-i,1))/2.0);
    ll(N-i+1,N-i-1+1) = 0.5*k(N-i,1);
    ll(N-i+1,N-i+1+1) = 0.5*k(N-i+1,1);
    ll(N-i+1,N-i+N+1+1) = (k(N-i+1,1) + k(N-i,1))/2.0;
    
    i = i+1;
end


i = 1;
while i<0.3*N
    ll(N*(N+1)+i+1,N*(N+1)+i+1) = -0.5*k(i-1+1,N-1+1) -0.5*k(i+1,N-1+1) - ((k(i-1+1,N-1+1) + k(i+1,N-1+1))/2.0);
    ll(N*(N+1)+i+1,N*(N+1)+i-1+1) = 0.5*k(i-1+1,N-1+1);
    ll(N*(N+1)+i+1,N*(N+1)+i+1+1) = 0.5*k(i+1,N-1+1);
    ll(N*(N+1)+i+1,N*(N+1)+i-(N+1)+1) = (k(i-1+1,N-1+1) + k(i+1,N-1+1))/2.0;

    ll(N*(N+2)-i+1,N*(N+2)-i+1) = -0.5*k(N-i+1,N-1+1) -0.5*k(N-1-i+1,N-1+1) - ((k(N-i+1,N-1+1) + k(N-1-i+1,N-1+1))/2.0);
    ll(N*(N+2)-i+1,N*(N+2)-i-1+1) = 0.5*k(N-1-i+1,N-1+1);
    ll(N*(N+2)-i+1,N*(N+2)-i+1+1) = 0.5*k(N-i+1,N-1+1);
    ll(N*(N+2)-i+1,N*(N+2)-i-(N+1)+1) = (k(N-i+1,N-1+1) + k(N-1-i+1,N-1+1))/2.0;

    i = i+1;
end

%i = 2;
%while i<(N+1)
for i = 2:1:N
   
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+1) = -0.5*k(1,i-1) -0.5*k(1,i) - ((k(1,i-1) + k(1,i))/2.0);
    ll((i-1)*(N+1)+1,(i-1)*(N+1)+1+1) = (k(1,i-1) + k(1,i))/2.0;
    ll((i-1)*(N+1)+1,(i-2)*(N+1)+1) = 0.5*k(1,i-1);
    ll((i-1)*(N+1)+1,(i)*(N+1)+1) = 0.5*k(1,i);

    ll(i*(N+1),i*(N+1)) = -0.5*k(N,i-1) -0.5*k(N,i) - ((k(N,i-1) + k(N,i))/2.0);
    ll(i*(N+1),i*(N+1)-1) = (k(N,i-1) + k(N,i))/2.0;
    ll(i*(N+1),(i-1)*(N+1)) = 0.5*k(N,i-1);
    ll(i*(N+1),(i+1)*(N+1)) = 0.5*k(N,i);
    
end

for j = 2:1:N
    for i = 2:1:N
        
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i) = -((k(i-1,j-1)+k(i-1,j))/2.0 + (k(i,j-1)+k(i,j))/2.0 + (k(i-1,j-1)+k(i,j-1))/2.0 + (k(i-1,j)+k(i,j))/2.0);
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i-1) = (k(i-1,j-1)+k(i-1,j))/2.0;
        ll((j-1)*(N+1)+i,(j-1)*(N+1)+i+1) = (k(i,j-1)+k(i,j))/2.0;
        ll((j-1)*(N+1)+i,(j-2)*(N+1)+i) = (k(i-1,j-1)+k(i,j-1))/2.0;
        ll((j-1)*(N+1)+i,j*(N+1)+i) = (k(i-1,j)+k(i,j))/2.0;
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

