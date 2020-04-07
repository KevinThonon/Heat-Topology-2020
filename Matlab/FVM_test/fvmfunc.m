function [solmat] = fvmfunc(N,i)
%gegeven
%q = 2/(0.01*0.01*0.001);
%2e orde testcase
%q = 10000/(0.01*0.01);
%veelvoud van 10 nemen
d = 0.01;
h = 0.01/N; 


dp = N+1; %discretisatiepunten waarin de temperatuur gemeten zal worden


X = linspace(0,0.01,dp);
Y = linspace(0,0.01,dp);
%% set up Kmatrix (dimensie is 1 kleiner dan matrix met temperaturen)
%k = 0.4*65+0.5*0.2;

kmat = ones(dp-1,dp-1); % in de optimalisatiestap worden de k waardes op elke node aangepast

% k1 = dp/2-0.5;
% k2 = dp/2+0.5;
% kmat((k1):k1,:) = 1;
% kmat(k2:(k2),:) = 1;

bigkmat = bigKmat(kmat);



%% fill in fvm matrix

% %LL*T = RL
% % T is de kolomvector met de te berekenen variabelen:
% % [HP1,N,N,D,D,D,N,N,HP2,N,X,X,X,X,X,X,X,N,N,X,X,X,X,X,X,X,N,...,HP3,N,N,D,D,D,N,N,HP4]


[LL,RL] = FVM_mat(i,h,dp,bigkmat);

[solmat,sol] = solver(dp,LL,RL);

surf(X,Y,solmat)
end

