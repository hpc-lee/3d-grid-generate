clc;
clear all;
close all;

% p1        p2        p5            p6
% ...........         ...............
%           .         .
%           .         .
%           ...........    
%           p3        p4

flag_printf = 1;

nx = 401;
ny = 301;

p1 = 1;
p2 = 101;
p3 = 151;
p4 = 251;
p5 = 301;
p6 = 401;

dh = 1;
origin_x = 0;
origin_y = 0;
origin_z = 0;

bz = zeros(ny,nx,3);
% NOTE: x direction tight
% so coord x incremental with x
for j=1:ny
  for i=p1:p2
    bz(j,i,1) = origin_x + (i-1)*dh;
    bz(j,i,2) = origin_y + (j-1)*dh;
    bz(j,i,3) = origin_z;
  end
  for i=p2+1:p3
    bz(j,i,1) = origin_x + bz(j,p2,1) + cos(0.5*pi*90/90)*(i-p2)*dh;
    bz(j,i,2) = origin_y + (j-1)*dh;
    bz(j,i,3) = origin_z + bz(j,p2,3) - sin(0.5*pi*90/90)*(i-p2)*dh;
  end
  for i=p3+1:p4
    bz(j,i,1) = origin_x + bz(j,p3,1) + (i-p3)*dh;
    bz(j,i,2) = origin_y + (j-1)*dh;
    bz(j,i,3) = origin_z + bz(j,p3,3);
  end
  for i=p4+1:p5
    bz(j,i,1) = origin_x + bz(j,p4,1) + cos(0.5*pi*90/90)*(i-p4)*dh;
    bz(j,i,2) = origin_y + (j-1)*dh;
    bz(j,i,3) = origin_z + bz(j,p4,3) + sin(0.5*pi*90/90)*(i-p4)*dh;
  end
  for i=p5+1:p6
    bz(j,i,1) = origin_x + bz(j,p5,1) + (i-p5)*dh;
    bz(j,i,2) = origin_y + (j-1)*dh;
    bz(j,i,3) = origin_z + bz(j,p5,3);
  end
end

% for j = 1:301
%     bz(j,:,3) = smooth(bz(j,:,3),20);
%     bz(j,:,1) = smooth(bz(j,:,1),20);
% end

if flag_printf
    figure(1)   
    surf(bz(:,:,1),bz(:,:,2),bz(:,:,3));
    axis equal;
    shading interp;
    set(gcf,'Position',[200,200,650,400]);
    xlabel('X axis (m)');
%     ylabel('Y axis (m)');
    zlabel('Z axis (m)');
    colorbar;
    view(10,30);
    %xlim([-100,700]);
    %ylim([-400,100]);
%     title('concave step model',FontSize=15);
    set(gca,'layer','top');
    set(gca,'FontSize',10,FontWeight='bold');
    set(gcf,'color','white','renderer','painters');
end

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
    fprintf(fid,'%.9e %.9e %.9e\n',bz(j,i,1),bz(j,i,2),bz(j,i,3));
  end
end

if flag_printf
   print(gcf,'model3.png','-r300','-dpng');
 end

