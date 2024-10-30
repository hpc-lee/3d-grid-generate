% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

num_pml = 00;
nx1 = 401;
ny1 = 401;

nx = nx1+2*num_pml;
ny = ny1+2*num_pml;
dx = 10;
dy = 10;

nz = 201;
dz = 10;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(nx,ny,3);
bz2 = zeros(nx,ny,3);

for j=1:ny1
  for i=1:nx1
    bz2(i+num_pml,j+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,j+num_pml,2) = origin_y + (j-1)*dy;
    bz2(i+num_pml,j+num_pml,3) = origin_z;
  end
end

if flag_topo_z
  point_x= origin_x + floor(nx1/2)*dx; 
  point_y= origin_y + floor(ny1/2)*dy; 
  L = 400;
  H = 600;
  for j = 1:ny1
    for i = 1:nx1
      r1 = (bz2(i+num_pml,j+num_pml,1)-point_x)^2 + (bz2(i+num_pml,j+num_pml,2)-point_y)^2;
      topo = H * exp(-r1/L^2);
      
      bz2(i+num_pml,j+num_pml,3) = bz2(i+num_pml,j+num_pml,3) - topo;
    end
  end
end

% [bz2] = extend_abs_layer(bz2,dx,dy,nx,ny,num_pml);

for j=1:ny
  for i=1:nx
    bz1(i,j,1) = bz2(i,j,1);
    bz1(i,j,2) = bz2(i,j,2);
    bz1(i,j,3) = origin_z - (nz-1)*dz;
  end
end

% A = 0.000001;
% [bz2] = arc_stretch_x(A,bz2);

% figure(1)   
% %plot3(bz1(:,:,1),bz1(:,:,2),bz1(:,:,3));
% surf(bz1(:,:,1),bz1(:,:,2),bz1(:,:,3));
% hold on;
% %plot3(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3));
% surf(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3));
% axis equal;
% shading interp;

figure(2)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
hold on;

plot3(bz2(150,200,1)/1e3, bz2(150,200,2)/1e3, bz2(150,200,3)/1e3+0.1, 'rv',MarkerSize=15);
plot3(bz2(200,200,1)/1e3, bz2(200,200,2)/1e3, bz2(200,200,3)/1e3+0.1, 'rv',MarkerSize=15);
plot3(bz2(250,200,1)/1e3, bz2(250,200,2)/1e3, bz2(250,200,3)/1e3+0.1, 'rv',MarkerSize=15);
plot3(bz2(300,200,1)/1e3, bz2(300,200,2)/1e3, bz2(300,200,3)/1e3+0.1, 'rv',MarkerSize=15);
plot3(bz2(350,200,1)/1e3, bz2(350,200,2)/1e3, bz2(350,200,3)/1e3+0.1, 'rv',MarkerSize=15);
text(bz2(150,200,1)/1e3-0.1, bz2(150,200,2)/1e3, bz2(150,200,3)/1e3+0.2,'R1',FontSize=15,FontWeight='bold');
text(bz2(200,200,1)/1e3-0.1, bz2(200,200,2)/1e3, bz2(200,200,3)/1e3+0.2,'R2',FontSize=15,FontWeight='bold');
text(bz2(250,200,1)/1e3-0.1, bz2(250,200,2)/1e3, bz2(250,200,3)/1e3+0.2,'R3',FontSize=15,FontWeight='bold');
text(bz2(300,200,1)/1e3-0.1, bz2(300,200,2)/1e3, bz2(300,200,3)/1e3+0.2,'R4',FontSize=15,FontWeight='bold');
text(bz2(350,200,1)/1e3-0.1, bz2(350,200,2)/1e3, bz2(350,200,3)/1e3+0.2,'R5',FontSize=15,FontWeight='bold');
plot3(1,1,0.1,'rp',MarkerSize=15);

axis equal;
cid = colorbar;
cid.Label.String='(km)';
caxis([-0.6,0]);
set(gca,'layer','top');
set(gcf,'color','white');
set(gcf,'Position',[0,0,1920,1200]);
set(gca, FontWeight='bold', FontSize=15);
xlabel('x-axis (km)',FontSize=15,FontWeight='bold',rotation=-40,position=[2.5 -0.8]);
ylabel('y-axis (km)',FontSize=15,FontWeight='bold',rotation=50, position=[4.7 1.5]);
zlabel('z-axis (km)',FontSize=15,FontWeight='bold',rotation=0,position=[-0.4 0.2]);
view(40,60);

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
   print(gcf,'model1.png','-r300','-dpng');
end
