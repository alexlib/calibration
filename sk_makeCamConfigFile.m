function sk_makeCamConfigFile(camParaCalib,fname)
% Expects a 1xN struct coming out of calibTsai.
% Writes a text file that can be read by the 3D tracking code.
% Now write the configuration file
s=size(camParaCalib);
ncams=s(2);
fid = fopen(fname, 'w');
fprintf(fid, '# PTV experiment configuration file\n');
fprintf(fid,['# Generated ' date '\n\n']);
fprintf(fid, '\n %d\t# ncams\n\n', ncams);
for icam = 1:ncams
    fprintf(fid, '######## cam %d ########\n', icam-1);
    fprintf(fid, '%d\t\t\t# Npix_x\n', camParaCalib(icam).Npixw);
    fprintf(fid, '%d\t\t\t# Npix_y\n', camParaCalib(icam).Npixh);
    fprintf(fid, '%11.8f\t\t# pixsize_x (mm)\n', camParaCalib(icam).wpix);
    fprintf(fid, '%11.8f\t\t# pixsize_y (mm)\n', camParaCalib(icam).hpix);
    fprintf(fid, '%12.8f\t\t# effective focal length (mm)\n', camParaCalib(icam).f_eff);
    % Note the sign change for k1, because its meaning is different in calib_Tsai and the stereomatching code
    fprintf(fid, '0\t\t# radial distortion kr (1/pixel)\n');
    fprintf(fid, '1\t\t# radial distortion kx (1/pixel)\n');
    % possible that this is backwards (rows v cols)...
%     fprintf(fid, '%12.8f\t# R\n', camParaCalib(icam).R); v1
%     fprintf(fid, '%12.8f\t# R\n', camParaCalib(icam).R'); v2
    fprintf(fid, '%12.8f\t# R\n', camParaCalib(icam).R');
    fprintf(fid, '%15.8f\t# T\n', camParaCalib(icam).T);
%     fprintf(fid, '%12.8f\t# Rinv\n', camParaCalib(icam).Rinv); v1
%     fprintf(fid, '%12.8f\t# Rinv\n', camParaCalib(icam).Rinv'); v2
    fprintf(fid, '%12.8f\t# Rinv\n', camParaCalib(icam).Rinv');
    fprintf(fid, '%15.8f\t# Tinv\n', camParaCalib(icam).Tinv);
    fprintf(fid, '\n');
end
fprintf(fid, '#### parameter for 3D matching ####');
fprintf(fid, '\n\n');
fprintf(fid, '0.1\t# mindist_pix (pixel)');
fprintf(fid, '\n');
fprintf(fid, '1\t# mindist_3D (mm)');
fclose(fid);
disp('Don''t forget to specify the 3D matching parameters!');

end