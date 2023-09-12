function [bdry] = arc_strech(A,bdry)

npoints = size(bdry,1);
x_0 = bdry(:,1);
z_0 = bdry(:,2);
% cal arc length
s=zeros(npoints,1);
for i=2:npoints
  s(i) = s(i-1) + sqrt((x_0(i)-x_0(i-1))^2 + (z_0(i)-z_0(i-1))^2);
end

% normalized 
u=zeros(npoints,1);
for i=2:npoints
  u(i) = s(i) / s(npoints); 
end

for i=2:npoints-1
  xi = (i-1)/(npoints-1);
  r = (exp(A*xi)-1)/(exp(A)-1);
  for m =1:npoints-1
    if(r>=u(m) && r<u(m+1))
      n = m;  
      break;
    end
  end
  bdry(i,1) = x_0(n) + (x_0(n+1)-x_0(n))*(r-u(n))/(u(n+1)-u(n));
  bdry(i,2) = z_0(n) + (z_0(n+1)-z_0(n))*(r-u(n))/(u(n+1)-u(n));
end
