function makeCamConfigFile(camParaCalib,fname)
% Expects a 1xN struct coming out of calibTsai.
% Writes a text file that can be read by the 3D tracking code.

fid=fopen(fname,'w');

fprintf(fid,'# Camera configuration file\n');
fprintf(fid,['# Generated ' date '\n\n']);

s=size(camParaCalib);

fprintf(fid,'%d\t\t# number of cameras\n\n',s(2));

for ii=1:s(2)
  fprintf(fid,'######## camera %d ########\n',ii-1);
  fprintf(fid,'%d\t\t\t# Npix_x\n',camParaCalib(ii).Npixw);
  fprintf(fid,'%d\t\t\t# Npix_y\n',camParaCalib(ii).Npixh);
  fprintf(fid,'%g\t\t# pixsize_x\n',camParaCalib(ii).wpix);
  fprintf(fid,'%g\t\t# pixsize_y\n',camParaCalib(ii).hpix);
  fprintf(fid,'%g\t\t# f_eff\n',camParaCalib(ii).f_eff);
  fprintf(fid,'%g\t\t# kr\n',0);
  fprintf(fid,'%g\t\t# kx\n',1);
  fprintf(fid,'%g\t\t# R\n',camParaCalib(ii).R');
  fprintf(fid,'%g\t\t# T\n',camParaCalib(ii).T);
  fprintf(fid,'%g\t\t# Rinv\n',camParaCalib(ii).Rinv');
  fprintf(fid,'%g\t\t# Tinv\n',camParaCalib(ii).Tinv);
  fprintf(fid,'\n');
end

fprintf(fid,'##### 3D match parameters #####\n\n');
fprintf(fid,'%g\t\t# mindist_pix\n', 0.1);
fprintf(fid,'%g\t\t# mindist_3D\n', 1);

disp('Don''t forget to specify the 3D matching parameters!');

fclose(fid);
end