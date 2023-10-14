function  [bx1,bx2,by1,by2] = arc_strech_zt(A,bx1,bx2,by1,by2)

% bx1
nz = size(bx1,1);
ny = size(bx1,2);
x = squeeze(bx1(:,:,1));
y = squeeze(bx1(:,:,2));
z = squeeze(bx1(:,:,3));

s = zeros(nz,ny);
for k=2:nz
  for j=1:ny
    s(k,j) = s(k-1,j) + sqrt((x(k,j)-x(k-1,j))^2 ...
            + (y(k,j)-y(k-1,j))^2 + (z(k,j)-z(k-1,j))^2);
  end
end
% normalized
u = zeros(nz,ny);
for k=2:nz
  for j=1:ny
    u(k,j) = s(k,j)/s(nz,j);
  end
end

for k=2:nz-1
  for j=1:ny
    zt = (k-1)/(nz-1);
    r = (exp(A*zt)-1)/(exp(A)-1);
    for m=1:nz-1
      if(r>=u(m,j) && r<u(m+1,j))
        n=m;
        break;
      end
    end
    bx1(k,j,1) = x(n,j) + (x(n+1,j)-x(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
    bx1(k,j,2) = y(n,j) + (y(n+1,j)-y(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
    bx1(k,j,3) = z(n,j) + (z(n+1,j)-z(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
  end
end

% bx2
nz = size(bx2,1);
ny = size(bx2,2);
x = squeeze(bx2(:,:,1));
y = squeeze(bx2(:,:,2));
z = squeeze(bx2(:,:,3));

s = zeros(nz,ny);
for k=2:nz
  for j=1:ny
    s(k,j) = s(k-1,j) + sqrt((x(k,j)-x(k-1,j))^2 ...
            + (y(k,j)-y(k-1,j))^2 + (z(k,j)-z(k-1,j))^2);
  end
end
% normalized
u = zeros(nz,ny);
for k=2:nz
  for j=1:ny
    u(k,j) = s(k,j)/s(nz,j);
  end
end

for k=2:nz-1
  for j=1:ny
    zt = (k-1)/(nz-1);
    r = (exp(A*zt)-1)/(exp(A)-1);
    for m=1:nz-1
      if(r>=u(m,j) && r<u(m+1,j))
        n=m;
        break;
      end
    end
    bx2(k,j,1) = x(n,j) + (x(n+1,j)-x(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
    bx2(k,j,2) = y(n,j) + (y(n+1,j)-y(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
    bx2(k,j,3) = z(n,j) + (z(n+1,j)-z(n,j))*(r-u(n,j))/(u(n+1,j)-u(n,j));
  end
end

% by1
nz = size(by1,1);
nx = size(by1,2);
x = squeeze(by1(:,:,1));
y = squeeze(by1(:,:,2));
z = squeeze(by1(:,:,3));

s = zeros(nz,nx);
for k=2:nz
  for i=1:nx
    s(k,i) = s(k-1,i) + sqrt((x(k,i)-x(k-1,i))^2 ...
            + (y(k,i)-y(k-1,i))^2 + (z(k,i)-z(k-1,i))^2);
  end
end
% normalized
u = zeros(nz,nx);
for k=2:nz
  for i=1:nx
    u(k,i) = s(k,i)/s(nz,i);
  end
end

for k=2:nz-1
  for i=1:nx
    zt = (k-1)/(nz-1);
    r = (exp(A*zt)-1)/(exp(A)-1);
    for m=1:nz-1
      if(r>=u(m,i) && r<u(m+1,i))
        n=m;
        break;
      end
    end
    by1(k,i,1)  = x(n,i) + (x(n+1,i)-x(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    by1(k,i,2)  = y(n,i) + (y(n+1,i)-y(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    by1(k,i,3)  = z(n,i) + (z(n+1,i)-z(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
  end
end

% by2
nz = size(by2,1);
nx = size(by2,2);
x = squeeze(by2(:,:,1));
y = squeeze(by2(:,:,2));
z = squeeze(by2(:,:,3));

s = zeros(nz,nx);
for k=2:nz
  for i=1:nx
    s(k,i) = s(k-1,i) + sqrt((x(k,i)-x(k-1,i))^2 ...
            + (y(k,i)-y(k-1,i))^2 + (z(k,i)-z(k-1,i))^2);
  end
end
% normalized
u = zeros(nz,nx);
for k=2:nz
  for i=1:nx
    u(k,i) = s(k,i)/s(nz,i);
  end
end

for k=2:nz-1
  for i=1:nx
    zt = (k-1)/(nz-1);
    r = (exp(A*zt)-1)/(exp(A)-1);
    for m=1:nz-1
      if(r>=u(m,i) && r<u(m+1,i))
        n=m;
        break;
      end
    end
    by2(k,i,1)  = x(n,i) + (x(n+1,i)-x(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    by2(k,i,2)  = y(n,i) + (y(n+1,i)-y(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    by2(k,i,3)  = z(n,i) + (z(n+1,i)-z(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
  end
end
