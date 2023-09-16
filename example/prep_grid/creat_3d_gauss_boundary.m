% creat 3D boundary data file
% nx ny nz x1(left) x2(right) y1(front) y2(back) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 0;
flag_topo_y = 1;
flag_topo_z = 0;

num_pml = 00;
nx1 = 200;
ny1 = 210;

nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
nz = 101;

dx = 100;
dy = 100;
dz = 100;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(ny,nx,3);
bz2 = zeros(ny,nx,3);
by1 = zeros(nz,nx,3);
by2 = zeros(nz,nx,3);
bx1 = zeros(nz,ny,3);
bx2 = zeros(nz,ny,3);

for j=1:ny1
  for i=1:nx1
    bz1(j+num_pml,i+num_pml,1) = origin_x + (i-1)*dx;
    bz1(j+num_pml,i+num_pml,2) = origin_y + (j-1)*dy;
    bz1(j+num_pml,i+num_pml,3) = origin_z - (nz-1)*dz;

    bz2(j+num_pml,i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(j+num_pml,i+num_pml,2) = origin_y + (j-1)*dy;
    bz2(j+num_pml,i+num_pml,3) = origin_z;
  end
end

if flag_topo_z
  point_x= origin_x + floor(nx1/2)*dx; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 40*dx;
  H = 30*dz;
  for j = 1:ny1
    for i = 1:nx1
      r1 = sqrt((bz2(j+num_pml,i+num_pml,1)-point_x)^2 + (bz2(j+num_pml,i+num_pml,2)-point_y)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bz2(j+num_pml,i+num_pml,3) = bz2(j+num_pml,i+num_pml,3) + topo;
    end
  end
end

[bz1] = extend_abs_layer(bz1,dx,dy,nx,ny,num_pml);
[bz2] = extend_abs_layer(bz2,dx,dy,nx,ny,num_pml);

for j=1:ny
  dz1 = (bz2(j,1,3)-bz1(j,1,3))/(nz-1);
  dz2 = (bz2(j,nx,3)-bz1(j,nx,3))/(nz-1);
  for k=1:nz
    bx1(k,j,1) = bz1(j,1,1);
    bx1(k,j,2) = bz1(j,1,2);
    bx1(k,j,3) = bz1(j,1,3) + (k-1)*dz1;

    bx2(k,j,1) = bz1(j,nx,1);
    bx2(k,j,2) = bz1(j,nx,2);
    bx2(k,j,3) = bz1(j,nx,3) + (k-1)*dz2;
  end
end

for i=1:nx
  dz1 = (bz2(1,i,3)-bz1(1,i,3))/(nz-1);
  dz2 = (bz2(ny,i,3)-bz1(ny,i,3))/(nz-1);
  for k=1:nz
    by1(k,i,1) = bz1(1,i,1);
    by1(k,i,2) = bz1(1,i,2);
    by1(k,i,3) = bz1(1,i,3) + (k-1)*dz1;

    by2(k,i,1) = bz1(ny,i,1);
    by2(k,i,2) = bz1(ny,i,2);
    by2(k,i,3) = bz1(ny,i,3) + (k-1)*dz2;
  end
end

if flag_topo_x
  point_z= origin_z - floor(nz/2)*dz; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 40*dz;
  H = 20*dx;
  for k = 1:nz
    for j = 1:ny1
      r1 = sqrt((bx2(k,j+num_pml,2)-point_y)^2 + (bx2(k,j+num_pml,3)-point_z)^2);
      topo = 0;
      if(r1 < L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      bx2(k,j+num_pml,1) = bx2(k,j+num_pml,1) + topo;
    end
  end
end

if flag_topo_y
  point_z= origin_z - floor(nz/2)*dz; 
  point_x= origin_x + floor(nx1/2)*dx; 
  L = 40*dz;
  H = 20*dy;
  for k = 1:nz
    for i = 1:nx1
      r1 = sqrt((by2(k,i+num_pml,1)-point_x)^2 + (by2(k,i+num_pml,3)-point_z)^2);
      topo = 0;
      if(r1 < L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      by2(k,i+num_pml,2) = by2(k,i+num_pml,2) + topo;
    end
  end
end
A = 0.000001;
 [bx1] = arc_strech(A,bx1);
 [bx2] = arc_strech(A,bx2);
 [by1] = arc_strech(A,by1);
 [by2] = arc_strech(A,by2);
 [bz1] = arc_strech(A,bz1);
 [bz2] = arc_strech(A,bz2);

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
