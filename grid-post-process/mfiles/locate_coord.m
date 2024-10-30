function coordinfo = locate_coord(parfnm,output_dir,subs,subc,subt)

% check parameter file exist
if ~ exist(parfnm,'file')
    error([mfilename ': file ' parfnm ' does not exist']);
end

% read parameters file
jsontext=fileread(parfnm);
par=jsondecode(jsontext);
num_of_grid = length(par.input_grids_info);

total_nx = 0;
total_ny = 0;
total_nz = 0;
if(strcmp(par.merge_direction,"x") == 1)
  for i=1:num_of_grid
    total_nx = total_nx + par.input_grids_info(i).number_of_grid_points(1);
    total_ny = par.input_grids_info(1).number_of_grid_points(2);
    total_nz = par.input_grids_info(1).number_of_grid_points(3);
  end
  total_nx = total_nx - num_of_grid + 1;
  ngijk=[total_nx,total_ny,total_nz];
end

if(strcmp(par.merge_direction,"y") == 1)
  for i=1:num_of_grid
    total_ny = total_ny + par.input_grids_info(i).number_of_grid_points(2);
    total_nx = par.input_grids_info(1).number_of_grid_points(1);
    total_nz = par.input_grids_info(1).number_of_grid_points(3);
  end
  total_ny = total_ny - num_of_grid + 1;
  ngijk=[total_nx,total_ny,total_nz];
end

if(strcmp(par.merge_direction,"z") == 1)
  for i=1:num_of_grid
    total_nz = total_nz + par.input_grids_info(i).number_of_grid_points(3);
    total_nx = par.input_grids_info(1).number_of_grid_points(1);
    total_ny = par.input_grids_info(1).number_of_grid_points(2);
  end
  total_nz = total_nz - num_of_grid + 1;
  ngijk=[total_nx,total_ny,total_nz];
end

if(par.flag_sample == 1)
  ngijk(1) = ngijk(1)*par.sample_factor_xi;
  ngijk(2) = ngijk(2)*par.sample_factor_et;
  ngijk(3) = ngijk(3)*par.sample_factor_zt;
end

gsubs = subs;
gsubt = subt;
gsubc = subc;

% reset count=-1 to total number
indx=find(subc==-1);
gsubc(indx)=ceil((ngijk(indx)-gsubs(indx)+1)./gsubt(indx));

gsube=gsubs+(gsubc-1).*gsubt;

% search the nc file headers to locate the threads/processors
coordprefix='coord';
coordlist=dir([output_dir,'/',coordprefix,'*.nc']);
n=0;
for i=1:length(coordlist)
    
    coordnm=[output_dir,'/',coordlist(i).name];
    xyzs=double(ncreadatt(coordnm, '/', 'global_index_of_first_physical_points'));
    xs=xyzs(1);
    ys=xyzs(2);
    zs=xyzs(3);
    xyzc=double(ncreadatt(coordnm, '/', 'count_of_physical_points'));
    xc=xyzc(1);
    yc=xyzc(2);
    zc=xyzc(3);
    xarray=[xs:xs+xc-1];
    yarray=[ys:ys+yc-1];
    zarray=[zs:zs+zc-1];
    if (length(find(xarray>=gsubs(1)-1 & xarray<=gsube(1)-1)) ~= 0 && ...
        length(find(yarray>=gsubs(2)-1 & yarray<=gsube(2)-1)) ~= 0 && ...
        length(find(zarray>=gsubs(3)-1 & zarray<=gsube(3)-1)) ~= 0)
        n=n+1;
        px(n)=str2num(coordlist(i).name( strfind(coordlist(i).name,'px' )+2 : ...
                                         strfind(coordlist(i).name,'_py')-1));
        py(n)=str2num(coordlist(i).name( strfind(coordlist(i).name,'py' )+2 : ...
                                         strfind(coordlist(i).name,'_pz')-1));
        pz(n)=str2num(coordlist(i).name( strfind(coordlist(i).name,'pz' )+2 : ...
                                         strfind(coordlist(i).name,'.nc')-1));
    end
    
end

nthd=0;
for ip=1:length(px)
    
    nthd=nthd+1;
    
    coordnm=[output_dir,'/',coordprefix,'_px',num2str(px(ip)),...
            '_py',num2str(py(ip)),'_pz',num2str(pz(ip)),'.nc'];
    xyzs=double(ncreadatt(coordnm, '/', 'global_index_of_first_physical_points'));
    xs=xyzs(1);
    ys=xyzs(2);
    zs=xyzs(3);
    xyzc=double(ncreadatt(coordnm, '/', 'count_of_physical_points'));
    xc=xyzc(1);
    yc=xyzc(2);
    zc=xyzc(3);
    xe=xs+xc-1;
    ye=ys+yc-1;
    ze=zs+zc-1;
    
    gxarray=gsubs(1):gsubt(1):gsube(1);
    gxarray=gxarray-1;
    gyarray=gsubs(2):gsubt(2):gsube(2);
    gyarray=gyarray-1;
    gzarray=gsubs(3):gsubt(3):gsube(3);
    gzarray=gzarray-1;
    
    coordinfo(nthd).thisid=[px(ip),py(ip),pz(ip)];
    i=find(gxarray>=xs & gxarray<=xe);
    j=find(gyarray>=ys & gyarray<=ye);
    k=find(gzarray>=zs & gzarray<=ze);
    coordinfo(nthd).indxs=[i(1),j(1),k(1)];
    coordinfo(nthd).indxe=[i(end),j(end),k(end)];
    coordinfo(nthd).indxc=coordinfo(nthd).indxe-coordinfo(nthd).indxs+1;
    
    coordinfo(nthd).subs=[ gxarray(i(1))-xs+1, ...
                           gyarray(j(1))-ys+1, ...
                           gzarray(k(1))-zs+1 ];
    coordinfo(nthd).subc=coordinfo(nthd).indxc;
    coordinfo(nthd).subt=gsubt;
    
    coordinfo(nthd).fnmprefix=coordprefix;
    
end

end

