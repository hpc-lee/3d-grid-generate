clc;
clear all;
close all;

nz = 200;
num_of_step = nz-1;

%for i=1:num_of_step
%  step(i) = 1;
%end

% Gradient Grid
% incre_ratio = 1.1;
% max_ratio = 3;
% incre_layer = ceil(log(max_ratio)/log(incre_ratio));
% 
% for i=1:10
%  step(i) = 1;
% end
% 
% for i=11:10+incre_layer
%  step(i) = step(i-1)*incre_ratio;
% end
% 
% for i=11+incre_layer:num_of_step
%  step(i) = step(10+incre_layer);
% end

% % Gradient Grid
incre_ratio = 1.1;
max_ratio = 3;
incre_layer = ceil(log(max_ratio)/log(incre_ratio));

for i=num_of_step:-1:num_of_step-9 
    step(i) = 1;
end
for i=num_of_step-10:-1:num_of_step-9-incre_layer
    step(i) = step(i+1)*incre_ratio;
end

for i=num_of_step-10-incre_layer:-1:1
    step(i) = step(num_of_step-9-incre_layer);
end


sum_step = sum(step);

% normalization step 
step_nor=step/sum_step;

arc_len(1) = 0;
for i=2:nz
  arc_len(i) = arc_len(i-1)+step_nor(i-1);
end
if((arc_len(nz)-1)>1e-8)
  error("step set is error, pelase check and reset");
end
% creat step file
file_name = './arc_len_file1.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of points\n'); 
fprintf(fid,'%d\n',nz);
fprintf(fid,'# arc_len\n'); 
for i=1:nz
  fprintf(fid,'%.9e \n',arc_len(i));
end
