clc;
clear all;
close all;

flag_printf = 1;
flag_topo_y = 1;

nx = 100;
nz = 110;

dx = 10;
dz = 10;
origin_x = 0;
origin_y = 0;
origin_z = 0;

by = zeros(nz,nx,3);
% NOTE: 
% Conventional right-hand 
% coordinate system
% coord y incremental with j
% coord z incremental with k
for k=1:nz
  for i=1:nx
    by(k,i,1) = origin_x + (i-1)*dx;
    by(k,i,2) = origin_y;
    by(k,i,3) = origin_z + (k-1)*dz;
  end
end

if  flag_topo_y
  point_x = origin_x + floor(nx/2)*dx; 
  point_z = origin_z + floor(nz/2)*dz; 
  L = 20*dx;
  H = 10*dx;
  for k = 1:nz
    for i = 1:nx
      r1 = sqrt((by(k,i,1)-point_x)^2 + (by(k,i,3)-point_z)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      by(k,i,2) = by(k,i,2) + topo;
    end
  end
end
A = 0.00001;
[by] = arc_strech(A,by);

if flag_printf
    figure(1)   
    plot3(permute(by(:,:,1),[2,1,3]),permute(by(:,:,2),[2,1,3]),permute(by(:,:,3),[2,1,3]));
    hold on;
    plot3(by(:,:,1),by(:,:,2),by(:,:,3));
    axis equal;
end

% trans direction y to z
% need keep Right-Handed Coordinate System
ny = nz;
bz = zeros(ny,nx,3);
bz(:,:,1) = by(:,:,1);
bz(:,:,2) = by(:,:,3);
bz(:,:,3) = by(:,:,2);
if flag_printf
    figure(2)   
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
