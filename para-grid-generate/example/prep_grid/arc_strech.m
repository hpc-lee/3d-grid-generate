function  [bdry] = arc_strech(A,bdry)

n2 = size(bdry,1);
n1 = size(bdry,2);
x_0 = squeeze(bdry(:,:,1));
y_0 = squeeze(bdry(:,:,2));
z_0 = squeeze(bdry(:,:,3));
x = squeeze(bdry(:,:,1));
y = squeeze(bdry(:,:,2));
z = squeeze(bdry(:,:,3));

% cal second dim arc length
s = zeros(n2,n1);
for j=1:n2
  for i=2:n1
    s(j,i) = s(j,i-1) + sqrt((x_0(j,i)-x_0(j,i-1))^2 + (y_0(j,i)-y_0(j,i-1))^2 + (z_0(j,i)-z_0(j,i-1))^2);
  end
end
% normalized
u = zeros(n2,n1);
for j=1:n2
  for i=2:n1
    u(j,i) = s(j,i)/s(j,n1);
  end
end

for j=1:n2
  for i=2:n1-1
    xi = (i-1)/(n1-1);
    r = (exp(A*xi)-1)/(exp(A)-1);
    for m=1:n1-1
      if(r>=u(j,m) && r<u(j,m+1))
        n=m;
        break;
      end
    end
    bdry(j,i,1)  = x_0(j,n) + (x_0(j,n+1)-x_0(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bdry(j,i,2)  = y_0(j,n) + (y_0(j,n+1)-y_0(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
    bdry(j,i,3)  = z_0(j,n) + (z_0(j,n+1)-z_0(j,n))*(r-u(j,n))/(u(j,n+1)-u(j,n));
  end
end

% then cal first dim arc length
s = zeros(n2,n1);
for j=2:n2
  for i=1:n1
    s(j,i) = s(j-1,i) + sqrt((x(j,i)-x(j-1,i))^2 + (y(j,i)-y(j-1,i))^2 + (z(j,i)-z(j-1,i))^2);
  end
end
% normalized
u = zeros(n2,n1);
for j=2:n2
  for i=1:n1
    u(j,i) = s(j,i)/s(n2,i);
  end
end

for j=2:n2-1
  for i=1:n1
    et = (j-1)/(n2-1);
    r = (exp(A*et)-1)/(exp(A)-1);
    for m=1:n2-1
      if(r>=u(m,i) && r<u(m+1,i))
        n=m;
        break;
      end
    end
    bdry(j,i,1) = x(n,i) + (x(n+1,i)-x(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bdry(j,i,2) = y(n,i) + (y(n+1,i)-y(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
    bdry(j,i,3) = z(n,i) + (z(n+1,i)-z(n,i))*(r-u(n,i))/(u(n+1,i)-u(n,i));
  end
end
