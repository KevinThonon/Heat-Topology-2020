N=100;

m = strcat('metal_0');
metal = strcat(m,'.txt');
fileID = fopen(metal,'r');
formatSpec = '%f';
metal_0 = fscanf(fileID,formatSpec);

m = strcat('metal_1');
metal = strcat(m,'.txt');
fileID = fopen(metal,'r');
formatSpec = '%f';
metal_1 = fscanf(fileID,formatSpec);

metal_mat_w = zeros(N,N);
metal_mat = zeros(N,N);

dp = N;
for i = 1:dp
    metal_mat_w(:,i)=metal_0(i*dp-dp+1:i*dp);
    metal_mat(:,i)=metal_1(i*dp-dp+1:i*dp);
    
end



