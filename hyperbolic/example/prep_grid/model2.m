clc;
clear all;
close all;

flag_printf = 1;

nx = 401;
ny = 301;

p1 = 1;
p2 = 81;
p3 = 131;
p4 = 271;
p5 = 321;
p6 = 401;

dh = 10;
origin_x = 0;
origin_y = 0;
origin_z = 0;

bz = zeros(nx,ny,3);

for j=1:ny
  for i=p1:p2
    bz(i,j,1) = origin_x + (i-1)*dh;
    bz(i,j,2) = origin_y + (j-1)*dh;
    bz(i,j,3) = origin_z;
  end
  for i=p2+1:p3
    bz(i,j,1) = origin_x + bz(p2,j,1) + cos(0.5*pi*90/90)*(i-p2)*dh;
    bz(i,j,2) = origin_y + (j-1)*dh;
    bz(i,j,3) = origin_z + bz(p2,j,3) + sin(0.5*pi*90/90)*(i-p2)*dh;
  end
  for i=p3+1:p4
    bz(i,j,1) = origin_x + bz(p3,j,1) + (i-p3)*dh;
    bz(i,j,2) = origin_y + (j-1)*dh;
    bz(i,j,3) = origin_z + bz(p3,j,3);
  end
  for i=p4+1:p5
    bz(i,j,1) = origin_x + bz(p4,j,1) + cos(0.5*pi*90/90)*(i-p4)*dh;
    bz(i,j,2) = origin_y + (j-1)*dh;
    bz(i,j,3) = origin_z + bz(p4,j,3) - sin(0.5*pi*90/90)*(i-p4)*dh;
  end
  for i=p5+1:p6
    bz(i,j,1) = origin_x + bz(p5,j,1) + (i-p5)*dh;
    bz(i,j,2) = origin_y + (j-1)*dh;
    bz(i,j,3) = origin_z + bz(p5,j,3);
  end
end

% for i=1:nx
%   bz(i,:,1) = smooth(bz(i,:,1),40);
%   bz(i,:,2) = smooth(bz(i,:,2),40);
%   bz(i,:,3) = smooth(bz(i,:,3),40);
% end
% 
% for j=1:ny
%  bz(:,j,1) = smooth(bz(:,j,1),40);
%  bz(:,j,2) = smooth(bz(:,j,2),40);
%  bz(:,j,3) = smooth(bz(:,j,3),40);
% end

figure(1)
surf(bz(:,:,1)/1e3,bz(:,:,2)/1e3,bz(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz(:,151,1)/1e3,bz(:,151,2)/1e3,bz(:,151,3)/1e3,'r',LineWidth=3);

%plot3(1.0, 0.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
%plot3(1.0, 1.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
%plot3(1.0, 2.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
%plot3(1.5, 1.5, -0.4, 'rv',MarkerSize=20,LineWidth=2);
%plot3(1.5, 0.5, -0.4, 'rv',MarkerSize=20,LineWidth=2);
%text(1.0, 0.5-0.1, +0.0+0.2,'R1',Color='r',FontSize=20,FontWeight='bold');
%text(1.0, 1.5-0.1, +0.0+0.2,'R2',Color='r',FontSize=20,FontWeight='bold');
%text(1.0, 2.5-0.1, +0.0+0.2,'R3',Color='r',FontSize=20,FontWeight='bold');
%text(1.5, 1.5-0.1, -0.4+0.2,'R4',Color='r',FontSize=20,FontWeight='bold');
%text(1.5, 0.5-0.1, -0.4+0.2,'R5',Color='r',FontSize=20,FontWeight='bold');
%plot3(1,1,0.0,'rp',MarkerSize=20,LineWidth=2);

grid off;
axis equal;
% caxis([0,0.6]);
cid = colorbar;
cid.Label.String='(km)';
%
set(gcf,'color','white');
set(gcf,'Position',[0,0,1900,1000]);
xlabel('X axis (km)',FontWeight='bold',rotation=-10,position=[2.0 -0.5]);
ylabel('Y axis (km)',FontWeight='bold',rotation=60, position=[3.4 1.0]);
zlabel('Z axis (km)',FontWeight='bold',rotation=0,position=[-0.3 0.9]);
set(gca,'layer','top');
set(gca, FontWeight='bold', FontSize=15);
view(20,30);

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
   print(gcf,'model2.png','-r400','-dpng');
 end

