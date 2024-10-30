% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 0;

num_pml = 00;
nx1 = 2000;
ny1 = 2000;

nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
dx = 10;
dy = 10;

nz = 1000;
dz = 10;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(ny,nx,3);
bz2 = zeros(ny,nx,3);

for j=1:ny1
  for i=1:nx1
    bz2(j+num_pml,i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(j+num_pml,i+num_pml,2) = origin_y + (j-1)*dy;
    bz2(j+num_pml,i+num_pml,3) = origin_z;
  end
end

if flag_topo_z
  point_x= origin_x + floor(nx1/2)*dx; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 400;
  H = 600;
  for j = 1:ny1
    for i = 1:nx1
      r1 = (bz2(j+num_pml,i+num_pml,1)-point_x)^2 + (bz2(j+num_pml,i+num_pml,2)-point_y)^2;
      topo = H * exp(-r1/L^2);
      
      bz2(j+num_pml,i+num_pml,3) = bz2(j+num_pml,i+num_pml,3) - topo;
    end
  end
end

[bz2] = extend_abs_layer(bz2,dx,dy,nx,ny,num_pml);

for j=1:ny
  for i=1:nx
    bz1(j,i,1) = bz2(j,i,1);
    bz1(j,i,2) = bz2(j,i,2);
    bz1(j,i,3) = origin_z - (nz-1)*dz;
  end
end

% A = 0.000001;
% [bz2] = arc_stretch(A,bz2);

if flag_printf
%     figure(1)   
%     plot3(bz1(:,:,1),bz1(:,:,2),bz1(:,:,3));
%     hold on;
%     plot3(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3));
%     axis equal;

   
end

% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx ny number\n'); 
fprintf(fid,'%d %d\n',nx,ny);
fprintf(fid,'# bz1 coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz1(j,i,1),bz1(j,i,2),bz1(j,i,3));
  end
end
fprintf(fid,'# bz2 coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz2(j,i,1),bz2(j,i,2),bz2(j,i,3));
  end
end

% 
% if flag_printf
%    print(gcf,'model1.png','-r300','-dpng');
% end
