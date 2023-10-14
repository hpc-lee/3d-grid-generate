clc;
clear all;
close all;

num_of_step = 50;
dh = -10;
for i=1:num_of_step
  step(i) = dh;
end
% creat step file
file_name = '../step_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step(i));
end
