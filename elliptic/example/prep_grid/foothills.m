% creat 3D boundary data file
% nx ny nz x1(left) x2(right) y1(front) y2(back) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
topo = importdata("topo_new.mat");

num_pml = 20;
nx1=501;
ny1=501;
nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
nz=301;


dx = 10;
dy = 10;
dz = 10;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(nx,ny,3);
bz2 = zeros(nx,ny,3);
by1 = zeros(nx,nz,3);
by2 = zeros(nx,nz,3);
bx1 = zeros(ny,nz,3);
bx2 = zeros(ny,nz,3);

i1 = 1;
i2 = i1 + nx1;
j1 = 501;
j2 = j1 + ny1;
topo_local = topo(i1:i2,j1:j2);

dx=10;
dy=10;
dz=10;
org_x = (i1-1)*dx;
org_y = (j1-1)*dy;

for j=1:ny1
  for i=1:nx1
    bz2(i+num_pml,j+num_pml,1) = org_x + (i-1)*dx;
    bz2(i+num_pml,j+num_pml,2) = org_y + (j-1)*dy;
    bz2(i+num_pml,j+num_pml,3) = topo_local(i,j);
  end
end
[bz2] = extend_abs_layer(bz2,dx,dy,nx,ny,num_pml);

for j=1:ny
  for i=1:nx
    bz1(j,i,1) = bz2(j,i,1);
    bz1(j,i,2) = bz2(j,i,2);
    bz1(j,i,3) = -(nz-1)*dz;
  end
end


for j=1:ny
  dz1 = (bz2(1,j,3)-bz1(1,j,3))/(nz-1);
  dz2 = (bz2(nx,j,3)-bz1(nx,j,3))/(nz-1);
  for k=1:nz
    bx1(j,k,1) = bz1(1,j,1);
    bx1(j,k,2) = bz1(1,j,2);
    bx1(j,k,3) = bz1(1,j,3) + (k-1)*dz1;

    bx2(j,k,1) = bz1(nx,j,1);
    bx2(j,k,2) = bz1(nx,j,2);
    bx2(j,k,3) = bz1(nx,j,3) + (k-1)*dz2;
  end
end

for i=1:nx
  dz1 = (bz2(i,1,3)-bz1(i,1,3))/(nz-1);
  dz2 = (bz2(i,ny,3)-bz1(i,ny,3))/(nz-1);
  for k=1:nz
    by1(i,k,1) = bz1(i,1,1);
    by1(i,k,2) = bz1(i,1,2);
    by1(i,k,3) = bz1(i,1,3) + (k-1)*dz1;

    by2(i,k,1) = bz1(i,ny,1);
    by2(i,k,2) = bz1(i,ny,2);
    by2(i,k,3) = bz1(i,ny,3) + (k-1)*dz2;
  end
end

if flag_printf

    figure(1)   
    plot3(bx1(:,:,1),bx1(:,:,2),bx1(:,:,3));
    hold on;
    plot3(bx2(:,:,1),bx2(:,:,2),bx2(:,:,3));
    axis equal;
 
    figure(2)   
    plot3(by1(:,:,1),by1(:,:,2),by1(:,:,3));
    hold on;
    plot3(by2(:,:,1),by2(:,:,2),by2(:,:,3));
    axis equal;
 
    figure(3)   
    plot3(bz1(:,:,1),bz1(:,:,2),bz1(:,:,3));
    hold on;
    plot3(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3));
    axis equal;

    figure(4)
    surf(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3),'edgecolor','none');
%     light;
    axis equal;
end

% creat data file
export_bdry;
