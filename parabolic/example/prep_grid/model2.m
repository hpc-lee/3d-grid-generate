% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;

nx = 401;
ny = 301;
dh = 10;
origin_x = 0;
origin_y = 0;
origin_z = 0;
bz1 = zeros(nx,ny,3);
bz2 = zeros(nx,ny,3);

p1 = 1;
p2 = 101;
p3 = 151;
p4 = 251;
p5 = 301;
p6 = 401;

for j=1:ny
  for i=p1:p2
    bz2(i,j,1) = origin_x + (i-1)*dh;
    bz2(i,j,2) = origin_y + (j-1)*dh;
    bz2(i,j,3) = origin_z;
  end
  for i=p2+1:p3
    bz2(i,j,1) = origin_x + bz2(p2,j,1) + cos(0.5*pi*70/90)*(i-p2)*dh;
    bz2(i,j,2) = origin_y + (j-1)*dh;
    bz2(i,j,3) = origin_z + bz2(p2,j,3) - sin(0.5*pi*70/90)*(i-p2)*dh;
  end
  for i=p3+1:p4
    bz2(i,j,1) = origin_x + bz2(p3,j,1) + (i-p3)*dh;
    bz2(i,j,2) = origin_y + (j-1)*dh;
    bz2(i,j,3) = origin_z + bz2(p3,j,3);
  end
  for i=p4+1:p5
    bz2(i,j,1) = origin_x + bz2(p4,j,1) + cos(0.5*pi*70/90)*(i-p4)*dh;
    bz2(i,j,2) = origin_y + (j-1)*dh;
    bz2(i,j,3) = origin_z + bz2(p4,j,3) + sin(0.5*pi*70/90)*(i-p4)*dh;
  end
  for i=p5+1:p6
    bz2(i,j,1) = origin_x + bz2(p5,j,1) + (i-p5)*dh;
    bz2(i,j,2) = origin_y + (j-1)*dh;
    bz2(i,j,3) = origin_z + bz2(p5,j,3);
  end
end

% smooth alone two direction
for i=1:nx
  bz2(i,:,1) = smooth(bz2(i,:,1),40);
  bz2(i,:,2) = smooth(bz2(i,:,2),40);
  bz2(i,:,3) = smooth(bz2(i,:,3),40);
end

for j=1:ny
 bz2(:,j,1) = smooth(bz2(:,j,1),40);
 bz2(:,j,2) = smooth(bz2(:,j,2),40);
 bz2(:,j,3) = smooth(bz2(:,j,3),40);
end

min(bz2(:,1,1))
max(bz2(:,1,1))
% min(bz2(:,1,2))
% max(bz2(:,1,2))
x=linspace(-500,3800,401);
y=linspace(0,3000,301);
for i=1:nx
  for j=1:ny
    bz1(i,j,1) = x(i);
    bz1(i,j,2) = y(j);
    bz1(i,j,3) = -2000;
  end
end


% A = 0.000001;
% [bz2] = arc_stretch_x(A,bz2);


figure(1)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
hold on;
plot3(bz2(:,151,1)/1e3,bz2(:,151,2)/1e3,bz2(:,151,3)/1e3,'r',LineWidth=2);

plot3(1.0, 0.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.0, 1.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.0, 2.5, +0.0, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.5, 1.5, -0.4, 'rv',MarkerSize=20,LineWidth=2);
plot3(1.5, 0.5, -0.4, 'rv',MarkerSize=20,LineWidth=2);
text(1.0, 0.5-0.1, +0.0+0.2,'R1',Color='r',FontSize=20,FontWeight='bold');
text(1.0, 1.5-0.1, +0.0+0.2,'R2',Color='r',FontSize=20,FontWeight='bold');
text(1.0, 2.5-0.1, +0.0+0.2,'R3',Color='r',FontSize=20,FontWeight='bold');
text(1.5, 1.5-0.1, -0.4+0.2,'R4',Color='r',FontSize=20,FontWeight='bold');
text(1.5, 0.5-0.1, -0.4+0.2,'R5',Color='r',FontSize=20,FontWeight='bold');
plot3(1,1,0.0,'rp',MarkerSize=20,LineWidth=2);

grid off;
axis equal;
% caxis([0,0.6]);
cid = colorbar;
cid.Label.String='(km)';
%
set(gcf,'color','white');
set(gcf,'Position',[0,0,2200,1000]);
xlabel('X axis (km)',FontWeight='bold',rotation=-25,position=[2.5 -1.3]);
ylabel('Y axis (km)',FontWeight='bold',rotation=40, position=[4.5 0.5]);
zlabel('Z axis (km)',FontWeight='bold',rotation=0,position=[-0.3 0.3]);
set(gca,'layer','top');
set(gca, FontWeight='bold', FontSize=15);
view(40,30);

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
   print(gcf,'model2.png','-r400','-dpng');
end
