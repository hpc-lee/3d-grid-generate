% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 1;

num_pml = 00;
ny1 = 300;
nz1 = 200;

ny = ny1+2*num_pml;
nz = nz1+2*num_pml;
dy = 100;
dz = 100;

nx = 400;
dx = 100;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bx1 = zeros(nz,ny,3);
bx2 = zeros(nz,ny,3);

for k=1:nz1
  for j=1:ny1
    bx2(k+num_pml,j+num_pml,1) = origin_x;
    bx2(k+num_pml,j+num_pml,2) = origin_y + (j-1)*dy;
    bx2(k+num_pml,j+num_pml,3) = origin_z + (k-1)*dz;
  end
end

if flag_topo_x
  point_y= origin_y + floor(ny1/2)*dy; 
  point_z= origin_z + floor(nz1/2)*dz; 
  L = 0.2*ny*dy;
  H = 0.15*nx*dx;
  for k = 1:nz1
    for j = 1:ny1
      r1 = sqrt((bx2(k+num_pml,j+num_pml,2)-point_y)^2 + (bx2(k+num_pml,j+num_pml,3)-point_z)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bx2(k+num_pml,j+num_pml,1) = bx2(k+num_pml,j+num_pml,1) - topo;
    end
  end
end

[bx2] = extend_abs_layer(bx2,dy,dz,ny,nz,num_pml);

for k=1:nz
  for j=1:ny
    bx1(k,j,1) = origin_x - (nx-1)*dx;
    bx1(k,j,2) = bx2(k,j,2);
    bx1(k,j,3) = bx2(k,j,3);
  end
end

A = 0.000001;
[bx2] = arc_stretch(A,bx2);

if flag_printf
    figure(1)   
    plot3(bx1(:,:,1),bx1(:,:,2),bx1(:,:,3));
    hold on;
    plot3(bx2(:,:,1),bx2(:,:,2),bx2(:,:,3));
    axis equal;
 
    figure(2)
    surf(bx2(:,:,1),bx2(:,:,2),bx2(:,:,3),'edgecolor','none');
%     light;
    axis equal;
end

% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# ny nz number\n'); 
fprintf(fid,'%d %d\n',ny,nz);
fprintf(fid,'# bx1 coords\n'); 
for k=1:nz
  for j=1:ny
    fprintf(fid,'%.9e %.9e %.9e\n',bx1(k,j,1),bx1(k,j,2),bx1(k,j,3));
  end
end
fprintf(fid,'# bx2 coords\n'); 
for k=1:nz
  for j=1:ny
    fprintf(fid,'%.9e %.9e %.9e\n',bx2(k,j,1),bx2(k,j,2),bx2(k,j,3));
  end
end

% 
% if flag_printf
%    print(gcf,'model1.png','-r300','-dpng');
% end
