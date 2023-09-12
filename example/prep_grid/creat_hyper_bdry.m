% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

nx1 = 250;
num_pml = 00;
nx = nx1 + 2*num_pml; 

dx = 10;
origin_x = 0;
origin_z = 0;
a_z=0.2;
H_z=0.2*nx*dx;

bz = zeros(nx,2);
for i=1:nx1
    bz(i+num_pml,1) = origin_x + (i-1)*dx;
    bz(i+num_pml,2) = origin_z;

end

if  flag_topo_z
    for i = 1:nx1
       bz(i+num_pml,2) = bz(i+num_pml,2) + H_z*exp(-((i-1)/(nx1-1)-0.5)^2/a_z^2);
    end
end

[bz] = extend_abs_layer(bz,dx,nx,num_pml);


if flag_printf
    figure(1)   
    plot(bz(:,1),bz(:,2));
    axis equal;
end

% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%g\n',nx);
fprintf(fid,'# bz coords\n'); 
for i=1:nx
  fprintf(fid,'%g %g\n',bz(i,1),bz(i,2));
end
