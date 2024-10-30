function [x,y,z] = gather_coord(coordinfo,output_dir)

% check
if ~ exist(output_dir,'dir')
    error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
coordprefix='coord';
nthd=length(coordinfo);
for n=1:nthd
    
    n_i=coordinfo(n).thisid(1); n_j=coordinfo(n).thisid(2); n_k=coordinfo(n).thisid(3);
    i1=coordinfo(n).indxs(1); j1=coordinfo(n).indxs(2); k1=coordinfo(n).indxs(3);
    i2=coordinfo(n).indxe(1); j2=coordinfo(n).indxe(2); k2=coordinfo(n).indxe(3);
    subs=coordinfo(n).subs;
    subc=coordinfo(n).subc;
    subt=coordinfo(n).subt;
    fnm_coord=[output_dir,'/',coordprefix,'_px',num2str(n_i),'_py',num2str(n_j),'_pz',num2str(n_k),'.nc'];
    
    if ~ exist(fnm_coord,'file')
       error([mfilename ': file ' fnm_coord 'does not exist']);
    end

    %- x coord
    x(i1:i2,j1:j2,k1:k2)=ncread(fnm_coord,'x',subs,subc,subt);

    %- x coord
    y(i1:i2,j1:j2,k1:k2)=ncread(fnm_coord,'y',subs,subc,subt);

    %- z coord
    z(i1:i2,j1:j2,k1:k2)=ncread(fnm_coord,'z',subs,subc,subt);
end

x = squeeze(x);
y = squeeze(y);
z = squeeze(z);

end
