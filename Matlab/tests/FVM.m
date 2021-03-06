clear all
close all


%gegeven
q = 2/(0.01*0.01*0.001);

N = 20; %veelvoud van 10 nemen
d = 0.01;
h = 0.01/N; 

dp = N+1; %discretisatiepunten waarin de temperatuur gemeten zal worden

%% set up Kmatrix (dimensie is 1 kleiner dan matrix met temperaturen)
k = 0.4;

kmat = 0.4*ones(dp-1,dp-1); % in de optimalisatiestap worden de k waardes op elke node aangepast



% maak grote kmatrix (bevat de gemiddeldes tussen elke node)
[s1,s2] = size(kmat);


bigkmat = zeros(s1*2+1,s2*2+1);
% bigkmat opvullen met kmat
for i = 2:2:(s1*2+1)
    for j = 2:2:(s2*2+1)
        bigkmat(i,j) = kmat(i/2,j/2);  
    end
end
% gemiddeldes nemen in bigkmat
% speciale gevallen:
% 1. inwendige punten (punten die tussen de waardes van de originele kmat
% liggen)
% 2. randpunten (die geen hoekpunt zijn)
% 3. hoekpunten?


% 1.

for i = 2:2:(s1*2) %(tijdens elke iteratie het gemiddelde nemen tussen element [i,j] en [i+2,j] en tussen [i,j] en [i,j+2] (element rechts en onder)
    for j = 2:2:(s2*2)
        if (i==s1*2)&&(j~=s2*2)
            bigkmat(i,j+1) = (bigkmat(i,j)+bigkmat(i,j+2))/2;
        elseif(j==s2*2)&&(i~=s1*2)
            bigkmat(i+1,j) = (bigkmat(i,j)+bigkmat(i+2,j))/2;
        elseif(i<=s1*2-2)&&(j<=s2*2-2)
            bigkmat(i+1,j) = (bigkmat(i,j)+bigkmat(i+2,j))/2;
            bigkmat(i,j+1) = (bigkmat(i,j)+bigkmat(i,j+2))/2;
        end
        
    end
end

% 2.
% eerste en laatste kolom
for i = 2:2:2*s1
    bigkmat(i,1) = bigkmat(i,2);
    bigkmat(i,(2*s2+1)) = bigkmat(i,(2*s2));
end

% eerste en laatste rij
for j = 2:2:2*s2
    bigkmat(1,j) = bigkmat(2,j);
    bigkmat((2*s1+1),j) = bigkmat((2*s1),j);
end

% de matrix bevat nu alle k(i-1/2,j), k(i+1/2,j),k(i,j-1/2) en k(i,j+1/2)
% die nodig zijn in de berekening van de flux in elke node. Beschouw de
% matrix als grid, nu voor elke "0" de FVM vergelijkingen opstellen en in
% een andere matrix steken

% test #nullen == dp*dp

% counter = 0;
% for i = 1:(s1*2+1)
%     for j = 1:(s2*2+1)
%         if (bigkmat(i,j)==0)
%             counter = counter+1;
%         end
%     end
% end
% 
% test = counter - dp^2


%% fill in fvm matrix

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



LL = zeros(dp^2,dp^2);

%col1 hardcoden in LL

%dirichletRVW in LL
dirll = eye(counter2);
i = 0.003/h+1;
j = 0.007/h+1;
LL(i:j,i:j)=dirll;
LL((dp*dp-j+1):(dp*dp-i+1),(dp*dp-j+1):(dp*dp-i+1)) = dirll;

%hoekpunten
%linskboven (1,1), LL(1,:)
LL(1,1) = -0.5*(bigkmat(2,1)+bigkmat(1,2));
LL(1,2) = 0.5* bigkmat(2,1);%south
LL(1,1+dp) = 0.5*bigkmat(1,2);%east
%linksonder (dp,1), LL(dp,:)
LL(dp,dp) = -0.5*(bigkmat((2*dp-1)-1,1)+bigkmat(2*dp-1,2));
LL(dp,dp-1) = 0.5*bigkmat((2*dp-1)-1,1);%north
LL(dp,dp+dp) = 0.5*bigkmat(2*dp-1,2);%east
%rechtsboven (1,dp), 
LL(dp*dp-dp+1,dp*dp-dp+1) = -0.5*(bigkmat(1,dp*2-1-1)+bigkmat(2,dp*2-1));
LL(dp*(dp-1)+1,dp*(dp-1)+1+1)=0.5*bigkmat(2,dp*2-1);%south
LL(dp*(dp-1)+1,dp*dp-dp+1-dp)= 0.5*bigkmat(1,dp*2-1-1);%west
%rechtsonder (dp,dp), LL(dp*dp,:)
LL(dp*dp,dp*dp) = -0.5*(bigkmat(dp*2-1,dp*2-1-1)+bigkmat(dp*2-1-1,dp*2-1));
LL(dp*dp,dp*dp-1) = 0.5*bigkmat(dp*2-1-1,dp*2-1);%north
LL(dp*dp,dp*dp-dp)=0.5*bigkmat(dp*2-1,dp*2-1-1);%west

%Neumann links & rechts
for j = 2:(dp-1)
    if(h*(j-1)<0.003)||(h*(j-1)>0.007)
        %LL(1:dp,1:2*dp)
        LL(j,j) = -(bigkmat(j*2-1,2)+0.5*bigkmat(j*2-1-1,1)+0.5*bigkmat(j*2-1+1,1));%east,north,south
        LL(j,j+dp) = bigkmat(j*2-1,2); %east
        LL(j,j-1) = 0.5*bigkmat(j*2-1-1,1);%north
        LL(j,j+1) = 0.5*bigkmat(j*2-1+1,1);%south
        
        %LL(dp*(dp-1)+1:dp*dp,dp*(dp-2)+1:dp*dp)
        LL(dp*dp-dp+j,dp*dp-dp+j) = -(bigkmat(j*2-1,dp*2-1-1)+0.5*bigkmat(j*2-1-1,2*dp-1)+0.5*bigkmat(j*2-1+1,2*dp-1));%west,north,south
        LL(dp*dp-dp+j,dp*dp-dp+j-dp) = bigkmat(j*2-1,dp*2-1-1); %west
        LL(dp*dp-dp+j,dp*dp-dp+j-1) = 0.5*bigkmat(j*2-1-1,2*dp-1);%north
        LL(dp*dp-dp+j,dp*dp-dp+j+1) = 0.5*bigkmat(j*2-1+1,2*dp-1);%south
        
    end
end

%Neumann boven & onder

for j = 1:(dp-2)
    k = j+1;
    %boven
    LL(j*dp+1,j*dp+1) = -(bigkmat(2,2*k-1)+0.5*bigkmat(1,2*k-1-1)+0.5*bigkmat(1,2*k-1+1));%south,west,east
    LL(j*dp+1,j*dp+1+1) =bigkmat(2,2*k-1);%south
    LL(j*dp+1,j*dp+1-dp) = 0.5*bigkmat(1,2*k-1-1);%west
    LL(j*dp+1,j*dp+1+dp) = 0.5*bigkmat(1,2*k-1+1);%east
    
    %onder
    LL(j*dp+dp,j*dp+dp)= -(bigkmat(2*dp-1-1,2*k-1)+0.5*bigkmat(2*dp-1,2*k-1-1)+0.5*bigkmat(2*dp-1,2*k-1+1));%north,west,east
    LL(j*dp+dp,j*dp+dp-1) = bigkmat(2*dp-1-1,2*k-1);%north
    LL(j*dp+dp,j*dp+dp-dp) = 0.5*bigkmat(2*dp-1,2*k-1-1);%west
    LL(j*dp+dp,j*dp+dp+dp) = 0.5*bigkmat(2*dp-1,2*k-1+1);%east
end

%gewonekolommen

for i = 2:(dp-1)
    for j = 2:(dp-1)
        k = (i-1)*dp+j;
        i2 = 2*i-1;
        j2 = 2*j-1;
        LL(k,k) = -(bigkmat(i2+1,j2)+bigkmat(i2-1,j2)+bigkmat(i2,j2+1)+bigkmat(i2,j2-1));%east,west,south,north 
        LL(k,k+dp)=bigkmat(i2+1,j2);%east
        LL(k,k-dp)=bigkmat(i2-1,j2);%west
        LL(k,k+1)=bigkmat(i2,j2+1);%south
        LL(k,k-1)=bigkmat(i2,j2-1);%north
    end
end



sol = LL\RL;

solmat = zeros(dp,dp);

for i = 1:dp
    solmat(:,i)=sol(i*dp-dp+1:i*dp);
end

surf(solmat)










