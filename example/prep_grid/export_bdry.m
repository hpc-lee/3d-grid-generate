% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# nz number\n'); 
fprintf(fid,'%d\n',nz);

fprintf(fid,'# bx1 coords\n'); 
for k=1:nz
  fprintf(fid,'%.9e %.9e\n',bx1(k,1),bx1(k,2));
end

fprintf(fid,'# bx2 coords\n'); 
for k=1:nz
  fprintf(fid,'%.9e %.9e\n',bx2(k,1),bx2(k,2));
end

fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz1(i,1),bz1(i,2));
end

fprintf(fid,'# bz2 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz2(i,1),bz2(i,2));
end
