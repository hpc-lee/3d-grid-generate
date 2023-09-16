clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

num_pml = 0;
nx1 = 100;
nx = nx1 + 2*num_pml; 

ny1 = 110;
ny = ny1 + 2*num_pml; 

dx = 10;
dy = 10;
origin_x = 0;
origin_y = 0;
origin_z = 0;

bz = zeros(ny,nx,3);
% NOTE: 
% Conventional right-hand 
% coordinate system
% coord x incremental with i
% coord y incremental with j
for j=1:ny1
  for i=1:nx1
    bz(j+num_pml,i+num_pml,1) = origin_x + (i-1)*dx;
    bz(j+num_pml,i+num_pml,2) = origin_y + (j-1)*dy;
    bz(j+num_pml,i+num_pml,3) = origin_z;
  end
end

if  flag_topo_z
  point_x= origin_x + floor(nx1/2)*dx; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 40*dx;
  H = 30*dx;
  for j = 1:ny1
    for i = 1:nx1
      r1 = sqrt((bz(j+num_pml,i+num_pml,1)-point_x)^2 + (bz(j+num_pml,i+num_pml,2)-point_y)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bz(j+num_pml,i+num_pml,3) = bz(j+num_pml,i+num_pml,3) + topo;
    end
  end
end
[bz] = extend_abs_layer(bz,dx,dy,nx,ny,num_pml);
A = 0.00001;
[bz] = arc_strech(A,bz);

if flag_printf
    figure(1)   
    plot3(permute(bz(:,:,1),[2,1,3]),permute(bz(:,:,2),[2,1,3]),permute(bz(:,:,3),[2,1,3]));
    hold on;
    plot3(bz(:,:,1),bz(:,:,2),bz(:,:,3));
    axis equal;
end

% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# ny number\n'); 
fprintf(fid,'%d\n',ny);
fprintf(fid,'# bz coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz(j,i,1),bz(j,i,2),bz(j,i,3));
  end
end
