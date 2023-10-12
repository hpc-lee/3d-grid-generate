function [x,y,z] = gather_coord(output_dir,subs,subc,subt)

% load
fnm_coord=[output_dir,'/','coord.nc'];
    
if ~ exist(fnm_coord,'file')
   error([mfilename ': file ' fnm_coord 'does not exist']);
end

xyzc = nc_attget(fnm_coord,nc_global,'number_of_points');
xyzc = double(xyzc);

if(subt(1)<0)
  xt = abs(subt(1));
  flip_subs = xyzc(1) - subs(1) + 1;
  if(subc(1) == -1)
    xc = ceil((xyzc(1)-flip_subs+1)/xt);
  else
    xc = subc(1);
  end
  xs = subs(1)-(xc-1)*xt-1;
end
if(subt(1)>0)
  xt = subt(1);
  if(subc(1) == -1)
    xc = ceil((xyzc(1)-subs(1)+1)/xt);
  else
    xc = subc(1);
  end
  xs = subs(1)-1;
end

if(subt(2)<0)
  yt = abs(subt(2));
  flip_subs = xyzc(2) - subs(2) + 1;
  if(subc(2) == -1)
    yc = ceil((xyzc(2)-flip_subs+1)/yt);
  else
    yc = subc(2);
  end
  ys = subs(2)-(yc-1)*yt-1;
end
if(subt(2)>0)
  yt = subt(2);
  if(subc(2) == -1)
    yc = ceil((xyzc(2)-subs(2)+1)/yt);
  else
    yc = subc(2);
  end
  ys = subs(2)-1;
end

if(subt(3)<0)
  zt = abs(subt(3));
  flip_subs = xyzc(3) - subs(3) + 1;
  if(subc(3) == -1)
    zc = ceil((xyzc(3)-flip_subs+1)/zt);
  else
    zc = subc(3);
  end
  zs = subs(3)-(zc-1)*zt-1;
end
if(subt(3)>0)
  zt = subt(3);
  if(subc(3) == -1)
    zc = ceil((xyzc(3)-subs(3)+1)/zt);
  else
    zc = subc(3);
  end
  zs = subs(3)-1;
end

i1 = 1;
i2 = i1 + xc - 1;
j1 = 1;
j2 = j1 + yc - 1;
k1 = 1;
k2 = k1 + zc - 1;

x(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'x',[zs,ys,xs],[zc,yc,xc],[zt,yt,xt]);
y(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'y',[zs,ys,xs],[zc,yc,xc],[zt,yt,xt]);
z(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'z',[zs,ys,xs],[zc,yc,xc],[zt,yt,xt]);
