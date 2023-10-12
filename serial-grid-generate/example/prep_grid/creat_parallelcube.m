clc;
clear
close all;

flag_printf = 1;
angle = 60;
theta = angle/180*pi;
slope = tan(theta);

nx = 310;
ny = 310;
nz = 300;

dx = 40;
dy = 40;
dz = 40;

origin_x = 0;
origin_y = 0;
origin_z = 0;

for j=1:ny
    for i=1:nx
      bz2(j,i,1) = origin_x + (i-1)*dx;
      bz2(j,i,2) = origin_y + (j-1)*dy;
      bz2(j,i,3) = origin_z;

      bz1(j,i,1) = origin_x + (i-1)*dx - (nz-1)*dz/slope;
      bz1(j,i,2) = origin_y + (j-1)*dy;
      bz1(j,i,3) = origin_z - (nz-1)*dz;
    end
end

for j=1:ny
    for k=1:nz
        bx1(k,j,1) = bz1(j,1,1) + (k-1)*dz/slope;
        bx1(k,j,2) = origin_y + (j-1)*dy;
        bx1(k,j,3) = bz1(j,1,3) + (k-1)*dz;

        bx2(k,j,1) = bz1(j,nx,1) + (k-1)*dz/slope;
        bx2(k,j,2) = origin_y + (j-1)*dy;
        bx2(k,j,3) = bz1(j,nx,3) + (k-1)*dz;
    end
end

for i=1:nx
    for k=1:nz
        by1(k,i,1) = origin_x + (i-1)*dx + bz1(1,i,1) + (k-1)*dz/slope;
        by1(k,i,2) = origin_y;
        by1(k,i,3) = bz1(1,i,3) + (k-1)*dz;

        by2(k,i,1) = origin_x + (i-1)*dx + bz1(ny,i,1) + (k-1)*dz/slope;
        by2(k,i,2) = origin_y + (ny-1)*dy;
        by2(k,i,3) = bz1(ny,i,3) + (k-1)*dz;
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
end

% creat data file
export_bdry;
