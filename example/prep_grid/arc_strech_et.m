function  [bx1,bx2,bz1,bz2] = arc_strech_et(A,bx1,bx2,bz1,bz2)

% bx1
nz = size(bx1,1);
ny = size(bx1,2);
x = squeeze(bx1(:,:,1));
y = squeeze(bx1(:,:,2));
z = squeeze(bx1(:,:,3));

% cal second dim arc length
s = zeros(nz,ny);
for k=1:nz
  for j=2:ny
    s(k,j) = s(k,j-1) + sqrt((x(k,j)-x(k,j-1))^2 ...
            + (y(k,j)-y(k,j-1))^2 + (z(k,j)-z(k,j-1))^2);
  end
end
% normalized
u = zeros(nz,ny);
for k=1:nz
  for j=2:ny
    u(k,j) = s(k,j)/s(k,ny);
  end
end

for k=1:nz
  for j=2:ny-1
    et = (j-1)/(ny-1);
    r = (exp(A*et)-1)/(exp(A)-1);
    for m=1:ny-1
      if(r>=u(k,m) && r<u(k,m+1))
        n=m;
        break;
      end
    end
    bx1(k,j,1) = x(k,n) + (x(k,n+1)-x(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    bx1(k,j,2) = y(k,n) + (y(k,n+1)-y(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    bx1(k,j,3) = z(k,n) + (z(k,n+1)-z(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
  end
end

% bx2
nz = size(bx2,1);
ny = size(bx2,2);
x = squeeze(bx2(:,:,1));
y = squeeze(bx2(:,:,2));
z = squeeze(bx2(:,:,3));

% cal second dim arc length
s = zeros(nz,ny);
for k=1:nz
  for j=2:ny
    s(k,j) = s(k,j-1) + sqrt((x(k,j)-x(k,j-1))^2 ...
            + (y(k,j)-y(k,j-1))^2 + (z(k,j)-z(k,j-1))^2);
  end
end
% normalized
u = zeros(nz,ny);
for k=1:nz
  for j=2:ny
    u(k,j) = s(k,j)/s(k,ny);
  end
end

for k=1:nz
  for j=2:ny-1
    et = (j-1)/(ny-1);
    r = (exp(A*et)-1)/(exp(A)-1);
    for m=1:ny-1
      if(r>=u(k,m) && r<u(k,m+1))
        n=m;
        break;
      end
    end
    bx2(k,j,1) = x(k,n) + (x(k,n+1)-x(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    bx2(k,j,2) = y(k,n) + (y(k,n+1)-y(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    bx2(k,j,3) = z(k,n) + (z(k,n+1)-z(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
  end
end

% bz1
ny = size(bz1,1);
nx = size(bz1,2);
x = squeeze(bz1(:,:,1));
y = squeeze(bz1(:,:,2));
z = squeeze(bz1(:,:,3));

% cal second dim arc length
s = zeros(ny,nx);
for j=2:ny
  for i=1:nx
    s(j,i) = s(j-1,i) + sqrt((x(j,i)-x(j-1,i))^2 ...
            + (y(j,i)-y(j-1,i))^2 + (z(j,i)-z(j-1,i))^2);
  end
end
% normalized
u = zeros(ny,nx);
for j=2:ny
  for i=1:nx
    u(j,i) = s(j,i)/s(ny,i);
  end
end

for j=2:ny-1
  for i=1:nx
    et = (j-1)/(ny-1);
    r = (exp(A*et)-1)/(exp(A)-1);
    for m=1:ny-1
      if(r>=u(m,i) && r<u(m+1,i))
        n=m;
        break;
      end
    end
    bz1(j,i,1) = x(n,i) + (x(n+1,i)-x(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bz1(j,i,2) = y(n,i) + (y(n+1,i)-y(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bz1(j,i,3) = z(n,i) + (z(n+1,i)-z(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
  end
end

% bz2
ny = size(bz2,1);
nx = size(bz2,2);
x = squeeze(bz2(:,:,1));
y = squeeze(bz2(:,:,2));
z = squeeze(bz2(:,:,3));

% cal second dim arc length
s = zeros(ny,nx);
for j=2:ny
  for i=1:nx
    s(j,i) = s(j-1,i) + sqrt((x(j,i)-x(j-1,i))^2 ...
            + (y(j,i)-y(j-1,i))^2 + (z(j,i)-z(j-1,i))^2);
  end
end
% normalized
u = zeros(ny,nx);
for j=2:ny
  for i=1:nx
    u(j,i) = s(j,i)/s(ny,i);
  end
end

for j=2:ny-1
  for i=1:nx
    et = (j-1)/(ny-1);
    r = (exp(A*et)-1)/(exp(A)-1);
    for m=1:ny-1
      if(r>=u(m,i) && r<u(m+1,i))
        n=m;
        break;
      end
    end
    bz2(j,i,1) = x(n,i) + (x(n+1,i)-x(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bz2(j,i,2) = y(n,i) + (y(n+1,i)-y(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bz2(j,i,3) = z(n,i) + (z(n+1,i)-z(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
  end
end

