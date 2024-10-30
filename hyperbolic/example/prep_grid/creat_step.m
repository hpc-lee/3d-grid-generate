clc;
clear all;
close all;

nz=301;
num_of_step = nz-1;
dh=-10;
for i=1:num_of_step
 step(i) = dh;
end
% for i=1:10
%   step(i) = -3;
% end
% for i=11:22
%   step(i) = step(i-1)*1.096;
% end
% for i=23:num_of_step
%   step(i) = step(22);
% end
% creat step file
file_name = '../step_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step(i));
end
