function [bz] = extend_abs_layer(bz,dx,dy,nx,ny,num_pml)

coef = 0.7;
A = -2;
% exponential stretching function
% xi interval is [0,1]
% (exp(A*xi)-1)/(exp(A)-1)


% x1 left

% first along pml normal direction(x),
% extend pml layer. alone normal,
% Height(Z) is basically equal
for i=num_pml:-1:1
  for j=1:ny
    dz = bz(j,i+2,3)-bz(j,i+1,3);
    slope = coef*dz/dx;
    bz(j,i,1) = bz(j,i+1,1) - dx;
    bz(j,i,2) = bz(j,i+1,2);
    bz(j,i,3) = bz(j,i+1,3) - slope*dx;
  end
end


% x2 right
for i=nx-num_pml+1:nx
  for j=1:ny
    dz = bz(j,i-1,3)-bz(j,i-2,3);
    slope = coef*dz/dx;
    bz(j,i,1) = bz(j,i-1,1) + dx;
    bz(j,i,2) = bz(j,i-1,2);
    bz(j,i,3) = bz(j,i-1,3) + slope*dx;
  end
end

% y1 front
for i=1:nx
  for j=num_pml:-1:1
    dz = bz(j+2,i,3)-bz(j+1,i,3);
    slope = coef*dz/dy;
    bz(j,i,1) = bz(j+1,i,1);
    bz(j,i,2) = bz(j+1,i,2) - dy;
    bz(j,i,3) = bz(j+1,i,3) - slope*dy;
  end
end

% y2 back
for i=1:nx
  for j=ny-num_pml+1:ny
    dz = bz(j-1,i,3)-bz(j-2,i,3);
    slope = coef*dz/dy;
    bz(j,i,1) = bz(j-1,i,1);
    bz(j,i,2) = bz(j-1,i,2) + dy;
    bz(j,i,3) = bz(j-1,i,3) + slope*dy;
  end
end

% along normal direction, reduce height.
% inner pml layer reduce more.
% this function is reduce pml 
% horizontal direction(y) angle.


%z2_avg = sum(bz2(:,1,3))/ny;
%z1_avg = sum(bz1(:,1,3))/ny;
%for j=1:ny
%  H2 = bz2(j,1,3) - z2_avg; 
%  H1 = bz1(j,1,3) - z1_avg; 
%  for i=num_pml:-1:1
%    xi = (i-1)/(num_pml-1); 
%    c(i) = (exp(A*xi)-1)/(exp(A)-1);
%
%    bz2(j,i,3) = bz2(j,1,3) - H2 + c(i)*H2;
%    bz1(j,i,3) = bz1(j,1,3) - H1 + c(i)*H1;
%  end
%end
%
%z2_avg = sum(bz2(:,nx,3))/ny;
%z1_avg = sum(bz1(:,nx,3))/ny;
%for j=1:ny
%  H2 = bz2(j,nx,3) - z2_avg; 
%  H1 = bz1(j,nx,3) - z1_avg; 
%  for i=nx-num_pml+1:nx
%    m = i-nx+num_pml;
%    xi = (m-1)/(num_pml-1); 
%    c(m) = (exp(-A*xi)-1)/(exp(-A)-1);
%
%    bz2(j,i,3) = bz2(j,nx-num_pml,3) - c(m)*H2;
%    bz1(j,i,3) = bz1(j,nx-num_pml,3) - c(m)*H1;
%  end
%end
%
%z2_avg = sum(bz2(1,:,3))/nx;
%z1_avg = sum(bz1(1,:,3))/nx;
%for i=1:nx
%  H2 = bz2(1,i,3) - z2_avg; 
%  H1 = bz1(1,i,3) - z1_avg; 
%  for j=num_pml:-1:1
%    et = (j-1)/(num_pml-1); 
%    c(j) = (exp(A*et)-1)/(exp(A)-1);
%
%    bz2(j,i,3) = bz2(1,i,3) - H2 + c(j)*H2;
%    bz1(j,i,3) = bz1(1,i,3) - H1 + c(j)*H1;
%  end
%end
%
%z2_avg = sum(bz2(ny,:,3))/nx;
%z1_avg = sum(bz1(ny,:,3))/nx;
%for i=1:nx
%  H2 = bz2(ny,i,3) - z2_avg; 
%  H1 = bz1(ny,i,3) - z1_avg; 
%  for j=ny-num_pml+1:ny
%    m = j-ny+num_pml;
%    et = (m-1)/(num_pml-1); 
%    c(m) = (exp(-A*et)-1)/(exp(-A)-1);
%
%    bz2(j,i,3) = bz2(ny-num_pml,i,3) - c(m)*H2;
%    bz1(j,i,3) = bz1(ny-num_pml,i,3) - c(m)*H1;
%  end
%end
