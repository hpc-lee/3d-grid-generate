% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_sample = 0;

topo = importdata("foothills_3d_topo.mat");

num_pml = 0;
nx1=701;
ny1=601;
nx = nx1+2*num_pml;
ny = ny1+2*num_pml;

i1 = 1;
i2 = i1 + nx1-1;
j1 = 1;
j2 = j1 + ny1-1;
topo_local = topo(i1:2:2*i2,j1:2:2*j2);

dx=20;
dy=20;
org_x = (i1-1)*dx;
org_y = (j1-1)*dy;

bz1 = zeros(nx,ny,3);
bz2 = zeros(nx,ny,3);

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
    bz1(i,j,1) = bz2(i,j,1);
    bz1(i,j,2) = bz2(i,j,2);
    bz1(i,j,3) = -3000;
  end
end

figure(2)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz2(11,:,1)/1e3,bz2(11,:,2)/1e3,bz2(11,:,3)/1e3,'r',LineWidth=3);

plot3(6, 2, 1+0.1, 'rv',MarkerSize=20,LineWidth=3);
plot3(6, 4, 1+0.1, 'rv',MarkerSize=20,LineWidth=3);
plot3(6, 6, 1+0.1, 'rv',MarkerSize=20,LineWidth=3);
plot3(6, 8, 1+0.1, 'rv',MarkerSize=20,LineWidth=3);
plot3(6, 10, 1+0.1, 'rv',MarkerSize=20,LineWidth=3);

text(6, 2+0.15, 1+0.5,'R1',Color='r',FontSize=20,FontWeight='bold');
text(6, 4+0.15, 1+0.5,'R2',Color='r',FontSize=20,FontWeight='bold');
text(6, 6+0.15, 1+0.5,'R3',Color='r',FontSize=20,FontWeight='bold');
text(6, 8+0.15, 1+0.5,'R4',Color='r',FontSize=20,FontWeight='bold');
text(6, 10+0.15, 1+0.5,'R5',Color='r',FontSize=20,FontWeight='bold');

plot3(10, 7, 0.5,'rp',MarkerSize=20,LineWidth=3);
axis equal;
set(gcf,'Position',[0,0,3000,1000]);
cid = colorbar;
cid.Label.String='(km)';
% caxis([0.19,1.6]);

set(gcf,'color','white');
xlabel('X axis (km)',rotation=20,position=[5.0 -2.0]);
ylabel('Y axis (km)',rotation=-25,position=[-2 4]);
zlabel('Z axis (km)',rotation=0,position=[3.5 17.0]);
view(-40,20);
grid off;
set(gca, FontWeight='bold', FontSize=20);

% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx ny number\n'); 
fprintf(fid,'%d %d\n',nx,ny);
fprintf(fid,'# bz1 coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz1(i,j,1),bz1(i,j,2),bz1(i,j,3));
  end
end
fprintf(fid,'# bz2 coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz2(i,j,1),bz2(i,j,2),bz2(i,j,3));
  end
end

if flag_printf
   print(gcf,'model5.png','-r400','-dpng');
end

