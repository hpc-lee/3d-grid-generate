clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 1;

ny = 400;
nz = 300;

dy = 10;
dz = 10;
origin_x = 0;
origin_y = 0;
origin_z = 0;

bx = zeros(nz,ny,3);
% NOTE: 
% Conventional right-hand 
% coordinate system
% coord x incremental with i
% coord z incremental with k
for k=1:nz
  for j=1:ny
    bx(k,j,1) = origin_x;
    bx(k,j,2) = origin_y + (j-1)*dy;
    bx(k,j,3) = origin_z + (k-1)*dz;
  end
end

if  flag_topo_x
  point_y = origin_y + floor(ny/2)*dy; 
  point_z = origin_z + floor(nz/2)*dz; 
  L = 0.3*ny*dy;
  H = 0.2*ny*dy;
  for k = 1:nz
    for j = 1:ny
      r1 = sqrt((bx(k,j,2)-point_y)^2 + (bx(k,j,3)-point_z)^2);
      topo = 0;
      if(r1 < L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bx(k,j,1) = bx(k,j,1) + topo;
    end
  end
end
A = 0.00001;

if flag_printf
    figure(1)   
    plot3(permute(bx(:,:,1),[2,1,3]),permute(bx(:,:,2),[2,1,3]),permute(bx(:,:,3),[2,1,3]));
    hold on;
    plot3(bx(:,:,1),bx(:,:,2),bx(:,:,3));
    axis equal;
end

% trans direction x to z
% need keep Right-Handed Coordinate System
%nx = nz;
%bz = zeros(ny,nx,3);
%bx = permute(bx,[2,1,3]);
%bz(:,:,1) = bx(:,:,3);
%bz(:,:,2) = bx(:,:,2);
%bz(:,:,3) = bx(:,:,1);
%if flag_printf
%    figure(2)   
%    plot3(permute(bz(:,:,1),[2,1,3]),permute(bz(:,:,2),[2,1,3]),permute(bz(:,:,3),[2,1,3]));
%    hold on;
%    plot3(bz(:,:,1),bz(:,:,2),bz(:,:,3));
%    axis equal;
%end


% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# ny number\n'); 
fprintf(fid,'%d\n',ny);
fprintf(fid,'# nz number\n'); 
fprintf(fid,'%d\n',nz);
fprintf(fid,'# bx coords\n'); 
for k=1:nz
  for j=1:ny
    fprintf(fid,'%.9e %.9e %.9e\n',bx(k,j,1),bx(k,j,2),bx(k,j,3));
  end
end
