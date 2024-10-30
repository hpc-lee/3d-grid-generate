% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_sample = 0;

topo = importdata("topo_new.mat");

num_pml = 20;
nx1=501;
ny1=501;
nx = nx1+2*num_pml;
ny = ny1+2*num_pml;

i1 = 1;
i2 = i1 + nx1;
j1 = 501;
j2 = j1 + ny1;
topo_local = topo(i1:i2,j1:j2);

nz = 301;
dx=10;
dy=10;
dz=10;
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
    bz1(i,j,3) = -(nz-1)*dz;
  end
end

if flag_sample
    nx_new = 641;
    ny_new = 641;
    bz2_new = zeros(nx_new,ny_new,3);
    x_line = linspace(bz2(1,1,1),bz2(nx,1,1),nx_new);
    y_line = linspace(bz2(1,1,2),bz2(1,ny,2),ny_new);  
    for j=1:ny_new
      for i=1:nx_new
        bz2_new(i,j,1) = x_line(i);
        bz2_new(i,j,2) = y_line(j);
      end
    end
    bz2_new(:,:,3) = griddata(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3),bz2_new(:,:,1),bz2_new(:,:,2),'cubic');

    for j=1:ny_new
      for i=1:nx_new
        bz1_new(i,j,1) = bz2_new(i,j,1);
        bz1_new(i,j,2) = bz2_new(i,j,2);
        bz1_new(i,j,3) = -(nz-1)*dz;
      end
    end
end

% figure(1)
% surf(bz2(:,:,1),bz2(:,:,2),bz2(:,:,3));
% hold on;
% surf(bz1(:,:,1),bz1(:,:,2),bz1(:,:,3));
% shading interp;

figure(2)
surf(bz2(:,:,1)/1e3,bz2(:,:,2)/1e3,bz2(:,:,3)/1e3);
shading interp;
% hold on;
% plot3(bz2(121,401,1)/1e3, bz2(121,401,2)/1e3, bz2(121,401,3)/1e3+0.1, 'rv',MarkerSize=15);
% plot3(bz2(201,401,1)/1e3, bz2(201,401,2)/1e3, bz2(201,401,3)/1e3+0.1, 'rv',MarkerSize=15);
% plot3(bz2(281,401,1)/1e3, bz2(281,401,2)/1e3, bz2(281,401,3)/1e3+0.1, 'rv',MarkerSize=15);
% plot3(bz2(361,401,1)/1e3, bz2(361,401,2)/1e3, bz2(361,401,3)/1e3+0.1, 'rv',MarkerSize=15);
% plot3(bz2(441,401,1)/1e3, bz2(441,401,2)/1e3, bz2(441,401,3)/1e3+0.1, 'rv',MarkerSize=15);
% 
% text(bz2(121,401,1)/1e3-0.1, bz2(121,401,2)/1e3, bz2(121,401,3)/1e3+0.15,'R1',FontSize=15,FontWeight='bold');
% text(bz2(201,401,1)/1e3-0.1, bz2(201,401,2)/1e3, bz2(201,401,3)/1e3+0.15,'R2',FontSize=15,FontWeight='bold');
% text(bz2(281,401,1)/1e3-0.1, bz2(281,401,2)/1e3, bz2(281,401,3)/1e3+0.15,'R3',FontSize=15,FontWeight='bold');
% text(bz2(361,401,1)/1e3-0.1, bz2(361,401,2)/1e3, bz2(361,401,3)/1e3+0.15,'R4',FontSize=15,FontWeight='bold');
% text(bz2(441,401,1)/1e3-0.1, bz2(441,401,2)/1e3, bz2(441,401,3)/1e3+0.15,'R5',FontSize=15,FontWeight='bold');
% 
% plot3(bz2(221,221,1)/1e3, bz2(221,221,2)/1e3, bz2(221,221,3)/1e3+0.25,'rp',MarkerSize=15);
axis equal;
set(gcf,'Position',[0,0,1920,1200]);
cid = colorbar;
cid.Label.String='(km)';
caxis([0.2,1.5]);
set(gca, FontWeight='bold', FontSize=15);
set(gcf,'color','white');
xlabel('x-axis (km)',FontSize=15,FontWeight='bold',rotation=-50,position=[3 4.0]);
ylabel('y-axis (km)',FontSize=15,FontWeight='bold',rotation=20,position=[5.9 7]);
zlabel('z-axis (km)',FontSize=15,FontWeight='bold',rotation=0,position=[-2.8 6.3]);
view(60,30)

if flag_sample == 0
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
     print(gcf,'foothills.png','-r300','-dpng');
  end
end

if flag_sample == 1
  % creat data file
  file_name = '../data_file_3d.txt';
  fid=fopen(file_name,'w'); % Output file name 
  fprintf(fid,'# nx ny number\n'); 
  fprintf(fid,'%d %d\n',nx_new,ny_new);
  fprintf(fid,'# bz1 coords\n'); 
  for j=1:ny_new
    for i=1:nx_new
      fprintf(fid,'%.9e %.9e %.9e\n',bz1_new(i,j,1),bz1_new(i,j,2),bz1_new(i,j,3));
    end
  end
  fprintf(fid,'# bz2 coords\n');
  for j=1:ny_new
    for i=1:nx_new
      fprintf(fid,'%.9e %.9e %.9e\n',bz2_new(i,j,1),bz2_new(i,j,2),bz2_new(i,j,3));
    end
  end
  
  if flag_printf
     print(gcf,'foothills.png','-r300','-dpng');
  end
end
  
if flag_printf
   print(gcf,'foothills.png','-r300','-dpng');
end
