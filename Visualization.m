N = 10;

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
temperature_mat = zeros(N,N);
metal_mat = zeros(N,N);

for i = 1:N
    for j = 1:N
         gradient_mat(i,j) = gradient((i-1)*N+j);
         temperature_mat(i,j) = temperature((i-1)*N+j);
         metal_mat(i,j) = metal(i + N*(j-1));
    end
end


%surf(gradient_mat,'FaceColor','interp')

%surf(temperature_mat,'FaceColor','interp')

surf(metal_mat,'FaceColor','interp')



