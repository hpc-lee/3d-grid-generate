clc;
clear all;
close all;

nz = 301;
num_of_step = nz-1;

flag_flip = 1;
% 
for i=1:num_of_step
 step(i) = 1;
end

% % Gradient Grid
% incre_layer = 12;
% max_ratio = 3;
% incre_ratio = exp((log(max_ratio)/incre_layer));
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

if(flag_flip)
  step=flip(step);
end

sum_step = sum(step);

% normalization step 
step_nor=step/sum_step;
sum_step_nor = sum(step_nor);
if((sum_step_nor-1)>1e-8)
  error("step set is error, pelase check and reset");
end
% creat step file
file_name = '../step_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step_nor(i));
end
