close all; 
clear all;
clc; 

addmypath;

flag_printf=1;
num_pml = 20;
%
[Vp,SegyTraceH,SegyH]=ReadSegy('../../../velocity.segy');
[row, col] = size(Vp);

dz=10;
dx=15;
% pcolor(Vp(1:5:row,1:5:col));

% obtain topography
for i = 1 : col
    for j = 1:row
        a = Vp(j,i);
        if (abs(a-4000.0)>0.1)
            break;
        end
    end
    z(i)=j*dz;
    x(i)=i*dx;
end

a=800:1400; % important topo zone
% a=900:1100;
% a=300:900;
x = x(a);
z = z(a);
z = smooth(z,10);

% plot(x(a),z(a),'r');
% hold on;
% plot(x(a),z1(a),'b');
% axis equal;

nx1 = length(a);
nx = nx1 + 2*num_pml;

nz = nx1; % optional
% nz = 100;

origin_x = min(x);
origin_z = min(z);

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
for i=1:nx1
    bz1(i+num_pml,1) = x(i);
    bz1(i+num_pml,2) = origin_z-(nz-1)*dz;

    bz2(i+num_pml,1) = x(i);
    bz2(i+num_pml,2) = z(i);
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);
dz1 = (bz2(1,2)-bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1);
    bx1(k,2) = bz1(1,2) + (k-1) * dz1;

    bx2(k,1) = bz1(nx,1);
    bx2(k,2) = bz1(nx,2) + (k-1) * dz2;
end

if flag_printf
    figure(1)   
    plot(bx1(:,1),bx1(:,2));
    hold on;
    plot(bx2(:,1),bx2(:,2));
    plot(bz1(:,1),bz1(:,2));
    plot(bz2(:,1),bz2(:,2));
    axis equal;
    set(gcf,'color','w');
end

% creat data file
export_bdry;
