% creat 3D boundary data file
% nx ny nz x1(left) x2(right) y1(front) y2(back) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;

nx = 200;
ny = 400;
nz = 200;

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
pi=3.1415926;
theta=45/180*pi;
for k=1:nz
  for j=1:ny
    bx1(k,j,1) = origin_x- (nz-k)*dz*tan(theta);
    bx1(k,j,2) = origin_y + (j-1)*dy;
    bx1(k,j,3) = origin_z - (nz-k)*dz;

    bx2(k,j,1) = origin_x + (nx-1)*dx- (nz-k)*dz*tan(theta);
    bx2(k,j,2) = origin_y + (j-1)*dy;
    bx2(k,j,3) = origin_z - (nz-k)*dz;

  end
end

if 0
for k = 1:nz
    for j = 1:ny
        r1 = sqrt((bx1(k,j,2)-10.5e3).^2 + (bx1(k,j,3)+7.5e3).^2);
        r2 = sqrt((bx1(k,j,2)-30.5e3).^2 + (bx1(k,j,3)+7.5e3).^2);
        fxy = 0;
        if(r1 <3e3)
            fxy = 600 * (1+cos(pi*r1/3e3));
        end
        
        if(r2 <3e3)
            fxy = 600 * (1+cos(pi*r2/3e3));
        end
        
        bx1(k,j,1)=bx1(k,j,1)+fxy;
    end
end
end

for j=1:ny
  dx1 = (bx2(1,j,1)-bx1(1,j,1))/(nx-1);
  dx2 = (bx2(nz,j,1)-bx1(nz,j,1))/(nx-1);
  for i=1:nx
    bz1(j,i,1) = bx1(1,j,1) + (i-1)*dx1;
    bz1(j,i,2) = bx1(1,j,2);
    bz1(j,i,3) = bx1(1,j,3);

    bz2(j,i,1) = bx1(nz,j,1) + (i-1)*dx2;
    bz2(j,i,2) = bx1(nz,j,2);
    bz2(j,i,3) = bx1(nz,j,3);
  end
end

for k=1:nz
  dx1 = (bx2(k,1,1) -bx1(k,1,1))/(nx-1);
  dx2 = (bx2(k,ny,1)-bx1(k,ny,1))/(nx-1);
  for i=1:nx
    by1(k,i,1) = bx1(k,1,1) + (i-1)*dx1;
    by1(k,i,2) = bx1(k,1,2);
    by1(k,i,3) = bx1(k,1,3);

    by2(k,i,1) = bx1(k,ny,1) + (i-1)*dx2;
    by2(k,i,2) = bx1(k,ny,2);
    by2(k,i,3) = bx1(k,ny,3);
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

    pcolor(bx1(:,:,2),bx1(:,:,3),bx1(:,:,1));
%     light;
    shading flat
    colorbar
%     axis image
    axis equal;
end

% creat data file
export_bdry;
