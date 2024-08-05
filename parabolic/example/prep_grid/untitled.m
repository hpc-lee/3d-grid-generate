clc;
clear all;
close all;
topo = load('./topo_srtm.mat');
x = topo.xp';
y = topo.yp';
z = topo.zp';

ny1 = size(x,1);
nx1 = size(x,2);


flag_printf = 1;

num_pml = 0;

nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
nz = 600;
dz = 50;

dx=50;
dy=50;

bz1 = zeros(ny,nx,3);
bz2 = zeros(ny,nx,3);
by1 = zeros(nz,nx,3);
by2 = zeros(nz,nx,3);
bx1 = zeros(nz,ny,3);
bx2 = zeros(nz,ny,3);

bz1(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,1) = x;
bz1(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,2) = y;
bz1(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,3) = -30000;

bz2(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,1) = x;
bz2(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,2) = y;
bz2(1+num_pml:ny1+num_pml,1+num_pml:nx1+num_pml,3) = z;

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


A = 0.000001;
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
