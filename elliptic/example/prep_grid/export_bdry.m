% creat data file
file_name = '../data_file_3d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# ny number\n'); 
fprintf(fid,'%d\n',ny);
fprintf(fid,'# nz number\n'); 
fprintf(fid,'%d\n',nz);

fprintf(fid,'# bx1 coords\n'); 
for k=1:nz
    for j=1:ny
        fprintf(fid,'%.9e %.9e %.9e\n',bx1(j,k,1),bx1(j,k,2),bx1(j,k,3));
    end
end

fprintf(fid,'# bx2 coords\n'); 
for k=1:nz
    for j=1:ny
        fprintf(fid,'%.9e %.9e %.9e\n',bx2(j,k,1),bx2(j,k,2),bx2(j,k,3));
    end
end

fprintf(fid,'# by1 coords\n'); 
for k=1:nz
    for i=1:nx
        fprintf(fid,'%.9e %.9e %.9e\n',by1(i,k,1),by1(i,k,2),by1(i,k,3));
    end
end

fprintf(fid,'# by2 coords\n'); 
for k=1:nz
    for i=1:nx
        fprintf(fid,'%.9e %.9e %.9e\n',by2(i,k,1),by2(i,k,2),by2(i,k,3));
    end
end

fprintf(fid,'# bz1 coords\n'); 
for j=1:ny
    for i=1:nx
        fprintf(fid,'%.9e %.9e %.9e\n',bz1(i,j,1),bz1(i,j,2),bz1(i,j,3));
    end
end

fprintf(fid,'# bz2 coords\n'); 
for j=1:ny
    for i=1:nx
        fprintf(fid,'%.9e %.9e %.9e\n',bz2(i,j,1),bz2(i,j,2),bz2(i,j,3));
    end
end
