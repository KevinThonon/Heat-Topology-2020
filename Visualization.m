f = fopen(metal.txt);
data = textscan(f,'%s');
fclose(f);
variable = str2double(data{1}(2:2:end));