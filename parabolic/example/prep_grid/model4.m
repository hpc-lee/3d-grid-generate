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
  point_x= origin_x + floor(nx/2)*dx; 
  point_y= origin_y + floor(ny/2)*dy; 
  L = 400;
  H = 2000;
  for j = 1:ny
    for i = 1:nx
      r1 = (bz2(i,j,1)-point_x)^2 + (bz2(i,j,2)-point_y)^2;

        topo = H * exp(-r1/L^2);
        if(topo>500)
            topo =500;
        end
        
      bz2(i,j,3) = bz2(i,j,3) - topo;
    end
  end
end

% % smooth alone two direction
for i=1:nx
  bz2(i,:,1) = smooth(bz2(i,:,1),10);
  bz2(i,:,2) = smooth(bz2(i,:,2),10);
  bz2(i,:,3) = smooth(bz2(i,:,3),10);
end

for j=1:ny
 bz2(:,j,1) = smooth(bz2(:,j,1),10);
 bz2(:,j,2) = smooth(bz2(:,j,2),10);
 bz2(:,j,3) = smooth(bz2(:,j,3),10);
end

x=linspace(-500,4500,401);
y=linspace(-500,4500,401);
for i=1:nx
  for j=1:ny
    bz1(i,j,1) = x(i);
    bz1(i,j,2) = y(j);
    bz1(i,j,3) = -2000;
  end
end

% A = 0.000001;
% [bz2] = arc_stretch_x(A,bz2);

% for j=2:ny
%     k(j) = (bz2(201,j,3)-bz2(201,j-1,3))/(bz2(201,j,2)-bz2(201,j-1,2));
% end
% figure(3)
% plot(k);
% 
% figure(2)
% plot(bz2(201,:,2)/1e3,bz2(201,:,3)/1e3,'r',LineWidth=2);
% axis equal;

figure(1)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz2(201,:,1)/1e3,bz2(201,:,2)/1e3,bz2(201,:,3)/1e3,'r',LineWidth=2);

plot3(1.8, 1.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.8, 1.5, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.8, 2.0, -0.3, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.8, 2.5, +0.1, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.8, 3.0, +0.1, 'rv',MarkerSize=20,LineWidth=2);
text(1.8-0.1, 1.0-0.1, +0.1+0.4,'R1',Color='r',FontSize=20,FontWeight='bold');
text(1.8-0.1, 1.5-0.1, +0.1+0.4,'R2',Color='r',FontSize=20,FontWeight='bold');
text(1.8-0.1, 2.0-0.1, -0.3+0.4,'R3',Color='r',FontSize=20,FontWeight='bold');
text(1.8-0.1, 2.5-0.1, +0.1+0.4,'R4',Color='r',FontSize=20,FontWeight='bold');
text(1.8-0.1, 3.0-0.1, +0.1+0.4,'R5',Color='r',FontSize=20,FontWeight='bold');
plot3(1,1.3,0.1,'rp',MarkerSize=20,LineWidth=2);

grid off;
axis equal;
caxis([-0.5,0]);
cid = colorbar;
cid.Label.String='(km)';
set(gca,'layer','top');
set(gcf,'color','white');
set(gcf,'Position',[0,0,1920,1200]);
set(gca, FontWeight='bold', FontSize=15);
xlabel('X axis (km)',FontWeight='bold',rotation=-40,position=[2.5 -0.8]);
ylabel('Y axis (km)',FontWeight='bold',rotation=50, position=[4.7 1.5]);
zlabel('Z axis (km)',FontWeight='bold',rotation=0,position=[-0.4 0.2]);
view(40,60);

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
   print(gcf,'model4.png','-r300','-dpng');
end
