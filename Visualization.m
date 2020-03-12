fileID = fopen('results.txt','r');
formatSpec = '%f';
size_sol = [1 Inf];
sol = fscanf(fileID,formatSpec,size_sol);

solmat = zeros(121,121);

for i = 1:dp
    solmat(:,i)=sol(i*dp-dp+1:i*dp);
end

surf(solmat,'FaceColor','interp')