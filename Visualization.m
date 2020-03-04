fileID = fopen('results.txt','r');
formatSpec = '%f';
sizeA = [1 Inf];
sol = fscanf(fileID,formatSpec,sizeA);

solmat = zeros(dp,dp);

for i = 1:dp
    solmat(:,i)=sol(i*dp-dp+1:i*dp);
end

surf(solmat,'FaceColor','interp')