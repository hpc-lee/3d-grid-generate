% creat 3D boundary data file
% nx ny nz x1(left) x2(right) y1(front) y2(back) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 0;
flag_topo_y = 0;
flag_topo_z = 1;

num_pml = 20;
nx1 = 501;
ny1 = 501;

nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
nz = 301;

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

for j=1:ny1
  for i=1:nx1
    bz1(i+num_pml,j+num_pml,1) = origin_x + (i-1)*dx;
    bz1(i+num_pml,j+num_pml,2) = origin_y + (j-1)*dy;
    bz1(i+num_pml,j+num_pml,3) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,j+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,j+num_pml,2) = origin_y + (j-1)*dy;
    bz2(i+num_pml,j+num_pml,3) = origin_z;
  end
end

if flag_topo_z
  point_x= origin_x + floor(nx1/2)*dx; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 0.3*nx*dx;
  H = 0.2*nx*dx;
  for j = 1:ny1
    for i = 1:nx1
      r1 = sqrt((bz2(i+num_pml,j+num_pml,1)-point_x)^2 + (bz2(i+num_pml,j+num_pml,2)-point_y)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bz2(i+num_pml,j+num_pml,3) = bz2(i+num_pml,j+num_pml,3) + topo;
    end
  end
end

[bz1] = extend_abs_layer(bz1,dx,dy,nx,ny,num_pml);
[bz2] = extend_abs_layer(bz2,dx,dy,nx,ny,num_pml);

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

A = 0.0001;
%[by1,by2,bz1,bz2] = arc_strech_xi(A,by1,by2,bz1,bz2);
%[bx1,bx2,bz1,bz2] = arc_strech_et(A,bx1,bx2,bz1,bz2);
%[bx1,bx2,by1,by2] = arc_strech_zt(A,bx1,bx2,by1,by2);

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
