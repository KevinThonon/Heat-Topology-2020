N = 70;

for iterations = 1:100

gradient = 'gradient.txt';
fileID = fopen(gradient,'r');
formatSpec = '%f';
gradient = fscanf(fileID,formatSpec);

metal = 'metal.txt';
fileID = fopen(metal,'r');
formatSpec = '%f';
metal = fscanf(fileID,formatSpec);

temperature = 'temperature.txt';
fileID = fopen(temperature,'r');
formatSpec = '%f';
temperature = fscanf(fileID,formatSpec);


gradient_mat = zeros(N,N);
metal_mat = zeros(N,N);

dp = N;
for i = 1:dp
    gradient_mat(:,i)=gradient(i*dp-dp+1:i*dp);
    metal_mat(:,i)=metal(i*dp-dp+1:i*dp);
    
end

temperature_mat = zeros(N+1,N+1);

dp = N+1;
for i = 1:dp
    temperature_mat(:,i)=temperature(i*dp-dp+1:i*dp);
end


surface(gradient_mat)

hold on
surface(temperature_mat)
hold on

surface(metal_mat)

pause(1)

end




