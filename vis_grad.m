N = 50;
fileID = fopen('gradient_3.txt','r');
formatSpec = '%f';
gradient = fscanf(fileID,formatSpec);

fileID = fopen('metal_360.txt','r');
formatSpec = '%f';
metal = fscanf(fileID,formatSpec);

fileID = fopen('temperature_512.txt','r');
formatSpec = '%f';
temperature = fscanf(fileID,formatSpec);

gradient_mat = zeros(N,N);
    metal_mat = zeros(N,N);
    temperature_mat = zeros(N+1,N+1);

dp = N;
for i = 1:dp
    gradient_mat(:,i)=gradient(i*dp-dp+1:i*dp);
    metal_mat(:,i)=metal(i*dp-dp+1:i*dp);
    
end

 dp = N+1;
    for i = 1:dp
        temperature_mat(:,i)=temperature(i*dp-dp+1:i*dp);
    end
    
%surface(gradient_mat,'FaceColor','interp')
surface(metal_mat,'FaceColor','interp')
%surface(temperature_mat)