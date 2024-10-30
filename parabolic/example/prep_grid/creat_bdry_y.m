% creat boundary data file bz1 bz2

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_y = 1;

num_pml = 00;
nx1 = 400;
nz1 = 200;

nx = nx1+2*num_pml;
nz = nz1+2*num_pml;
dx = 100;
dz = 100;

ny = 300;
dy = 100;

origin_x = 0;
origin_y = 0;
origin_z = 0;
by1 = zeros(nz,nx,3);
by2 = zeros(nz,nx,3);

for k=1:nz1
  for i=1:nx1
    by2(k+num_pml,i+num_pml,1) = origin_x + (i-1)*dx;
    by2(k+num_pml,i+num_pml,2) = origin_y;
    by2(k+num_pml,i+num_pml,3) = origin_z + (k-1)*dz;
  end
end

if flag_topo_y
  point_x= origin_x + floor(nx1/2)*dx; 
  point_z= origin_z + floor(nz1/2)*dz; 
  L = 0.2*nx*dx;
  H = 0.2*ny*dy;
  for k = 1:nz1
    for i = 1:nx1
      r1 = sqrt((by2(k+num_pml,i+num_pml,1)-point_x)^2 + (by2(k+num_pml,i+num_pml,3)-point_z)^2);
      topo = 0;
      if(r1 <L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      
      by2(k+num_pml,i+num_pml,2) = by2(k+num_pml,i+num_pml,2) - topo;
    end
  end
end

[by2] = extend_abs_layer(by2,dx,dz,nx,nz,num_pml);

for k=1:nz
  for i=1:nx
    by1(k,i,1) = by2(k,i,1);
    by1(k,i,2) = origin_y - (ny-1)*dy;
    by1(k,i,3) = by2(k,i,3);
  end
end

A = 0.000001;
[by2] = arc_stretch(A,by2);

if flag_printf
    figure(1)   
    plot3(by1(:,:,1),by1(:,:,2),by1(:,:,3));
    hold on;
    plot3(by2(:,:,1),by2(:,:,2),by2(:,:,3));
    axis equal;
 
    figure(2)
    surf(by2(:,:,1),by2(:,:,2),by2(:,:,3),'edgecolor','none');
%     light;
    axis equal;
end

% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx nz number\n'); 
fprintf(fid,'%d %d\n',nx,nz);
fprintf(fid,'# by1 coords\n'); 
for k=1:nz
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',by1(k,i,1),by1(k,i,2),by1(k,i,3));
  end
end
fprintf(fid,'# by2 coords\n'); 
for k=1:nz
  for i=1:nx
    fprintf(fid,'%.9e %.9e %.9e\n',by2(k,i,1),by2(k,i,2),by2(k,i,3));
  end
end

% 
% if flag_printf
%    print(gcf,'model1.png','-r300','-dpng');
% end
