function [bz] = extend_abs_layer(bz,dx,nx,num_pml)

coef = 0.7;
for i=num_pml:-1:1
    dz = bz(i+2,2)-bz(i+1,2);
    slope = coef*dz/dx;
    bz(i,1) = bz(i+1,1) - dx;
    bz(i,2) = bz(i+1,2) - slope*dx;

end

for i=nx-num_pml+1:nx
    dz = bz(i-1,2)-bz(i-2,2);
    slope2 = coef*dz/dx;
    bz(i,1) = bz(i-1,1) + dx;
    bz(i,2) = bz(i-1,2) + slope*dx;

end
%dz1 = bz1(i+2,2)-bz1(i+1,2);
%slope1 = coef*dz1/dx;
%bz1(i,1) = bz1(i+1,1) - dx;
%bz1(i,2) = bz1(i+1,2) - slope1*dx;
%dz1 = bz1(i-1,2)-bz1(i-2,2);
%slope1 = coef*dz1/dx;
%bz1(i,1) = bz1(i-1,1) + dx;
%bz1(i,2) = bz1(i-1,2) + slope1*dx;
