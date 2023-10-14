clc;
clear
close all;

flag_printf = 1;
angle = 30;
theta = angle/180*pi;
slope = tan(theta);
nx = 1010;
nz = 1000;
dx = 40;
dz = 40;
origin_x = 0;
origin_z = 0;
for i=1:nx
    bz2(i,1) = origin_x + (i-1)*dx;
    bz2(i,2) = origin_z;

    bz1(i,1) = origin_x + (i-1)*dx + (nz-1)*dz*cos(theta);
    bz1(i,2) = origin_z - (nz-1)*dz*sin(theta);
end

for k=1:nz
    bx1(k,1) = bz1(1,1) - (k-1)*dz*cos(theta);
    bx1(k,2) = bz1(1,2) + (k-1)*dz*sin(theta);

    bx2(k,1) = bz1(nx,1) - (k-1)*dz*cos(theta);
    bx2(k,2) = bz1(nx,2) + (k-1)*dz*sin(theta);
end

if flag_printf
    figure(1)   
    plot(bx1(:,1),bx1(:,2));
    hold on;
    plot(bx2(:,1),bx2(:,2));
    plot(bz1(:,1),bz1(:,2));
    plot(bz2(:,1),bz2(:,2));
    axis equal;
end

% creat data file
export_bdry;
