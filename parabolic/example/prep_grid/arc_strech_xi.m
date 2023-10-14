function  [by1,by2,bz1,bz2] = arc_strech_xi(A,by1,by2,bz1,bz2)

% by1
nz = size(by1,1);
nx = size(by1,2);
x = squeeze(by1(:,:,1));
y = squeeze(by1(:,:,2));
z = squeeze(by1(:,:,3));

% cal second dim arc length
s = zeros(nz,nx);
for k=1:nz
  for i=2:nx
    s(k,i) = s(k,i-1) + sqrt((x(k,i)-x(k,i-1))^2 ...
            + (y(k,i)-y(k,i-1))^2 + (z(k,i)-z(k,i-1))^2);
  end
end
% normalized
u = zeros(nz,nx);
for k=1:nz
  for i=2:nx
    u(k,i) = s(k,i)/s(k,nx);
  end
end

for k=1:nz
  for i=2:nx-1
    xi = (i-1)/(nx-1);
    r = (exp(A*xi)-1)/(exp(A)-1);
    for m=1:nx-1
      if(r>=u(k,m) && r<u(k,m+1))
        n=m;
        break;
      end
    end
    by1(k,i,1)  = x(k,n) + (x(k,n+1)-x(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    by1(k,i,2)  = y(k,n) + (y(k,n+1)-y(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    by1(k,i,3)  = z(k,n) + (z(k,n+1)-z(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
  end
end

% by2
nz = size(by2,1);
nx = size(by2,2);
x = squeeze(by2(:,:,1));
y = squeeze(by2(:,:,2));
z = squeeze(by2(:,:,3));

% cal second dim arc length
s = zeros(nz,nx);
for k=1:nz
  for i=2:nx
    s(k,i) = s(k,i-1) + sqrt((x(k,i)-x(k,i-1))^2 ...
            + (y(k,i)-y(k,i-1))^2 + (z(k,i)-z(k,i-1))^2);
  end
end
% normalized
u = zeros(nz,nx);
for k=1:nz
  for i=2:nx
    u(k,i) = s(k,i)/s(k,nx);
  end
end

for k=1:nz
  for i=2:nx-1
    xi = (i-1)/(nx-1);
    r = (exp(A*xi)-1)/(exp(A)-1);
    for m=1:nx-1
      if(r>=u(k,m) && r<u(k,m+1))
        n=m;
        break;
      end
    end
    by2(k,i,1) = x(k,n) + (x(k,n+1)-x(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    by2(k,i,2) = y(k,n) + (y(k,n+1)-y(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
    by2(k,i,3) = z(k,n) + (z(k,n+1)-z(k,n))*(r-u(k,n))/(u(k,n+1)-u(k,n));
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
for j=1:ny
  for i=2:nx
    s(j,i) = s(j,i-1) + sqrt((x(j,i)-x(j,i-1))^2 ...
            + (y(j,i)-y(j,i-1))^2 + (z(j,i)-z(j,i-1))^2);
  end
end
% normalized
u = zeros(ny,nx);
for j=1:ny
  for i=2:nx
    u(j,i) = s(j,i)/s(j,nx);
  end
end

for j=1:ny
  for i=2:nx-1
    xi = (i-1)/(nx-1);
    r = (exp(A*xi)-1)/(exp(A)-1);
    for m=1:nx-1
      if(r>=u(j,m) && r<u(j,m+1))
        n=m;
        break;
      end
    end
    bz1(j,i,1) = x(j,n) + (x(j,n+1)-x(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bz1(j,i,2) = y(j,n) + (y(j,n+1)-y(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bz1(j,i,3) = z(j,n) + (z(j,n+1)-z(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
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
for j=1:ny
  for i=2:nx
    s(j,i) = s(j,i-1) + sqrt((x(j,i)-x(j,i-1))^2 ...
            + (y(j,i)-y(j,i-1))^2 + (z(j,i)-z(j,i-1))^2);
  end
end
% normalized
u = zeros(ny,nx);
for j=1:ny
  for i=2:nx
    u(j,i) = s(j,i)/s(j,nx);
  end
end

for j=1:ny
  for i=2:nx-1
    xi = (i-1)/(nx-1);
    r = (exp(A*xi)-1)/(exp(A)-1);
    for m=1:nx-1
      if(r>=u(j,m) && r<u(j,m+1))
        n=m;
        break;
      end
    end
    bz2(j,i,1) = x(j,n) + (x(j,n+1)-x(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bz2(j,i,2) = y(j,n) + (y(j,n+1)-y(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bz2(j,i,3) = z(j,n) + (z(j,n+1)-z(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
  end
end
