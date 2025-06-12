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
bz = zeros(nx,ny,3);

for j=1:ny
  for i=1:nx
    bz(i,j,1) = origin_x + (i-1)*dx;
    bz(i,j,2) = origin_y + (j-1)*dy;
    bz(i,j,3) = origin_z;
  end
end

if flag_topo_z
  point_x= origin_x + floor(nx/2)*dx; 
  point_y= origin_y + floor(ny/2)*dy; 
  L = 200;
  H = 2000;
  for j = 1:ny
    for i = 1:nx
      r1 = (bz(i,j,1)-point_x)^2 + (bz(i,j,2)-point_y)^2;

        topo = H * exp(-r1/L^2);
        if(topo>500)
            topo =500;
        end
        
      bz(i,j,3) = bz(i,j,3) - topo;
    end
  end
end

% smooth alone two direction
for i=1:nx
  bz(i,:,1) = smooth(bz(i,:,1),20);
  bz(i,:,2) = smooth(bz(i,:,2),20);
  bz(i,:,3) = smooth(bz(i,:,3),20);
end

for j=1:ny
  bz(:,j,1) = smooth(bz(:,j,1),20);
  bz(:,j,2) = smooth(bz(:,j,2),20);
  bz(:,j,3) = smooth(bz(:,j,3),20);
end

% A = 0.000001;
% [bz] = arc_stretch(A,bz);


figure(1)
surf(bz(:,:,1)/1e3,bz(:,:,2)/1e3,bz(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz(:,201,1)/1e3,bz(:,201,2)/1e3,bz(:,201,3)/1e3,'r',LineWidth=3);

grid off;
axis equal;
%caxis([-0.5,0]);
cid = colorbar;
cid.Label.String='(km)';
set(gca,'layer','top');
set(gcf,'color','white');
set(gcf,'Position',[0,0,1900,1000]);
set(gca, FontWeight='bold', FontSize=15);
xlabel('X axis (km)',FontWeight='bold',rotation=-35,position=[2.7 -0.7]);
ylabel('Y axis (km)',FontWeight='bold',rotation=50, position=[4.8 1.5]);
zlabel('Z axis (km)',FontWeight='bold',rotation=0,position=[-0.4 0.2]);
view(40,70);

grid off;
% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# ny number\n'); 
fprintf(fid,'%d\n',ny);
fprintf(fid,'# bz coords\n'); 
for j=1:ny
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',bz(i,j,1),bz(i,j,2),bz(i,j,3));
  end
end

if flag_printf
   print(gcf,'model3.png','-r400','-dpng');
end

