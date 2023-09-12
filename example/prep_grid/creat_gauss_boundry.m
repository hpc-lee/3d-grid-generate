% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 1;
flag_topo_z = 1;

nx1 = 200;
nz = 200;
num_pml = 0;
nx = nx1 + 2*num_pml; 

dx = 10;
dz = 10;
origin_x = 0;
origin_z = 0;
a_x=0.20;
H_x=-0.2*nx*dx;
a_z=0.2;
H_z=0.3*nz*dz;

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
for i=1:nx1
    bz1(i+num_pml,1) = origin_x + (i-1)*dx;
    bz1(i+num_pml,2) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,2) = origin_z;
end

if  flag_topo_z
%     point = origin_x + floor(nx/2)*dx;
%     L = 50*dx;
%     H = 20*dz;
%     for i = 1:nx1
%         r1 = sqrt((bz2(i+num_pml,1)-point)^2);
%         topo = 0;
%         if(r1 < L)
%             topo = 0.5*H * (1+cos(pi*r1/L));
%         end
%         bz2(i+num_pml,2)=bz2(i+num_pml,2)+topo;
%     end
    for i = 1:nx1
       bz2(i+num_pml,2) = bz2(i+num_pml,2) + H_z*exp(-((i-1)/(nx1-1)-0.5)^2/a_z^2);
%        bz2(i+num_pml,2) = bz2(i+num_pml,2) + (i-1)*dx;
%        bz1(i+num_pml,2) = bz1(i+num_pml,2) + H*exp(-((i-1)/(nx1-1)-0.5)^2/a^2);
    end
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);

dz1 = (bz2(1,2)-bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1);
    bx1(k,2) = bz1(1,2) + (k-1)*dz1;

    bx2(k,1) = bz1(nx,1);
    bx2(k,2) = bz1(nx,2) + (k-1)*dz2;
end

if  flag_topo_x
  for k=1:nz
      topo_x(k) = H_x*exp(-((k-1)/(nz-1)-0.5)^2/a_x^2);
  end
  topo_x = topo_x - topo_x(1);
  for k=1:nz
      bx2(k,1) = bx2(k,1) - topo_x(k);
  end
  for k=1:nz
%       bx1(k,1) = bx1(k,1) + topo_x(k);
  end
end

A=0.0001;
[bx1]=arc_strech(A,bx1);
[bx2]=arc_strech(A,bx2);

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
