N = 100;
fileID = fopen('gradient_1.txt','r');
formatSpec = '%f';
gradient = fscanf(fileID,formatSpec);

fileID = fopen('metal_60.txt','r');
formatSpec = '%f';
metal = fscanf(fileID,formatSpec);

fileID = fopen('temperature_1.txt','r');
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
metal_mat = meshrefine(metal_mat,2);
surface(metal_mat,'FaceColor','interp')
%surface(temperature_mat)