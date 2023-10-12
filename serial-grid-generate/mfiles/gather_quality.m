function [v] = gather_quality(qualityinfo,output_dir,varnm)

% check
if ~ exist(output_dir,'dir')
    error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
varprefix=varnm;
nthd=length(qualityinfo);
for n=1:nthd
    
    n_i=qualityinfo(n).thisid(1); n_j=qualityinfo(n).thisid(2); n_k=qualityinfo(n).thisid(3);
    i1=qualityinfo(n).indxs(1); j1=qualityinfo(n).indxs(2); k1=qualityinfo(n).indxs(3);
    i2=qualityinfo(n).indxe(1); j2=qualityinfo(n).indxe(2); k2=qualityinfo(n).indxe(3);
    subs=qualityinfo(n).subs;
    subc=qualityinfo(n).subc;
    subt=qualityinfo(n).subt;
    fnm_var=[output_dir,'/',varprefix,'_px',num2str(n_i),'_py',num2str(n_j),'_pz',num2str(n_k),'.nc'];
    
    if ~ exist(fnm_var,'file')
       error([mfilename ': file ' fnm_var 'does not exist']);
    end

    subs=fliplr(subs);subc=fliplr(subc);subt=fliplr(subt);

    %- x coord
    v(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_var,varnm,subs-1,subc,subt);

end

v = squeeze(v);

end
