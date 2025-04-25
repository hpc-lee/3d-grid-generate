% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

nx = 401;
ny = 401;

dx = 10;
dy = 10;

nz = 201;
dz = 10;

origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(nx,ny,3);
bz2 = zeros(nx,ny,3);

for j=1:ny
  for i=1:nx
    bz2(i,j,1) = origin_x + (i-1)*dx;
    bz2(i,j,2) = origin_y + (j-1)*dy;
    bz2(i,j,3) = origin_z;
  end
end

if flag_topo_z
  x1=1000;
  y1=2000;
  x2=3000;
  y2=2000;
  L = 200;
  H = 500;
  for j = 1:ny
    for i = 1:nx
      r1 = (bz2(i,j,1)-x1)^2 + (bz2(i,j,2)-y1)^2;
      r2 = (bz2(i,j,1)-x2)^2 + (bz2(i,j,2)-y2)^2;

      topo = H * exp(-r1/L^2) - H * exp(-r2/L^2);
        
      bz2(i,j,3) = bz2(i,j,3) + topo;
    end
  end
end

% smooth alone two direction
for i=1:nx
 bz2(i,:,1) = smooth(bz2(i,:,1),20);
 bz2(i,:,2) = smooth(bz2(i,:,2),20);
 bz2(i,:,3) = smooth(bz2(i,:,3),20);
end

for j=1:ny
bz2(:,j,1) = smooth(bz2(:,j,1),20);
bz2(:,j,2) = smooth(bz2(:,j,2),20);
bz2(:,j,3) = smooth(bz2(:,j,3),20);
end

bz1(:,:,1) = bz2(:,:,1);
bz1(:,:,2) = bz2(:,:,2);
bz1(:,:,3) = -2000;

% A = 0.000001;
% [bz2] = arc_stretch_x(A,bz2);


figure(1)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz2(:,201,1)/1e3,bz2(:,201,2)/1e3,bz2(:,201,3)/1e3,'r',LineWidth=2);

plot3(1.0, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.5, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(2.0, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(2.5, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(3.0, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
text(1.0, 1.0-0.1, +0.1+0.2,'R1',Color='r',FontSize=20,FontWeight='bold');
text(1.5, 1.0-0.1, +0.1+0.2,'R2',Color='r',FontSize=20,FontWeight='bold');
text(2.0, 1.0-0.1, +0.1+0.2,'R3',Color='r',FontSize=20,FontWeight='bold');
text(2.5, 1.0-0.1, +0.1+0.2,'R4',Color='r',FontSize=20,FontWeight='bold');
text(3.0, 1.0-0.1, +0.1+0.2,'R5',Color='r',FontSize=20,FontWeight='bold');
plot3(3.0,2.5,0.1,'rp',MarkerSize=20,LineWidth=2);

grid off;
axis equal;
%caxis([-0.5,0]);
cid = colorbar;
cid.Label.String='(km)';
set(gca,'layer','top');
set(gcf,'color','white');
set(gcf,'Position',[0,0,2200,1200]);
set(gca, FontWeight='bold', FontSize=15);
xlabel('X axis (km)',FontWeight='bold',rotation=-25,position=[2.7 -1.2]);
ylabel('Y axis (km)',FontWeight='bold',rotation=35, position=[5.1 1.0]);
zlabel('Z axis (km)',FontWeight='bold',rotation=0,position=[-0.6 0.7]);
view(40,30);

grid off;
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
