N = 100;

for iterations = 1:50
    
    iter = int2str(iterations);

    g = strcat('gradient_',iter);
    gradient = strcat(g,'.txt');
    fileID = fopen(gradient,'r');
    formatSpec = '%f';
    gradient = fscanf(fileID,formatSpec);
    
  
    m = strcat('metal_',iter);
    metal = strcat(m,'.txt');
    fileID = fopen(metal,'r');
    formatSpec = '%f';
    metal = fscanf(fileID,formatSpec);
    
    t = strcat('temperature_',iter);
    temperature = strcat(t,'.txt');
    fileID = fopen(temperature,'r');
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
    
   
    metal_mat = meshrefine(metal_mat,2);
    surface(metal_mat,'FaceColor','interp')
    iterations
    pause(0.1)
    
    
end




