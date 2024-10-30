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
    dz = bz(i+2,j,3)-bz(i+1,j,3);
    slope = coef*dz/dx;
    bz(i,j,1) = bz(i+1,j,1) - dx;
    bz(i,j,2) = bz(i+1,j,2);
    bz(i,j,3) = bz(i+1,j,3) - slope*dx;
  end
end


% x2 right
for i=nx-num_pml+1:nx
  for j=1:ny
    dz = bz(i-1,j,3)-bz(i-2,j,3);
    slope = coef*dz/dx;
    bz(i,j,1) = bz(i-1,j,1) + dx;
    bz(i,j,2) = bz(i-1,j,2);
    bz(i,j,3) = bz(i-1,j,3) + slope*dx;
  end
end

% y1 front
for i=1:nx
  for j=num_pml:-1:1
    dz = bz(i,j+2,3)-bz(i,j+1,3);
    slope = coef*dz/dy;
    bz(i,j,1) = bz(i,j+1,1);
    bz(i,j,2) = bz(i,j+1,2) - dy;
    bz(i,j,3) = bz(i,j+1,3) - slope*dy;
  end
end

% y2 back
for i=1:nx
  for j=ny-num_pml+1:ny
    dz = bz(i,j-1,3)-bz(i,j-2,3);
    slope = coef*dz/dy;
    bz(i,j,1) = bz(i,j-1,1);
    bz(i,j,2) = bz(i,j-1,2) + dy;
    bz(i,j,3) = bz(i,j-1,3) + slope*dy;
  end
end
