% Now rewrite the configuration file
% restoredefaultpath;clear;close all

addpath('N:\skramel\water-tunnel\7-1-2016\calibration\')
load camParaCalib.mat

ncams = 4;
Npix_x = 1024;
Npix_y = 1280;

fid = fopen('camconfig.txt', 'w');
fprintf(fid, '# PTV experiment configuration file\n');
fprintf(fid, '\n %d\t# ncams\n\n', ncams);
for icam = 1:ncams
	fprintf(fid, '######## cam %d ########\n', icam-1);
	fprintf(fid, '%d\t\t\t# Npix_x\n', Npix_x);
	fprintf(fid, '%d\t\t\t# Npix_y\n', Npix_y);
	fprintf(fid, '%11.8f\t\t# pixsize_x (mm)\n', camParaCalib(icam).wpix);
	fprintf(fid, '%11.8f\t\t# pixsize_y (mm)\n', camParaCalib(icam).hpix);
	fprintf(fid, '%12.8f\t\t# effective focal length (mm)\n', camParaCalib(icam).f_eff);
	% Note the sign change for k1, because its meaning is different in calib_Tsai and the stereomatching code
	fprintf(fid, '%15.8e\t\t# radial distortion kr (1/pixel)\n', -(camParaCalib(icam).k1));
    fprintf(fid, '0\t\t# radial distortion kx (1/pixel)\n');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(1,1))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(1,2))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(1,3))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(2,1))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(2,2))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(2,3))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(3,1))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(3,2))');
    fprintf(fid, '%12.8f\t# R\n', (camParaCalib(icam).R(3,3))');
    fprintf(fid, '%15.8f\t# T\n', camParaCalib(icam).T(1,1));
    fprintf(fid, '%15.8f\t# T\n', camParaCalib(icam).T(2,1));
    fprintf(fid, '%15.8f\t# T\n', camParaCalib(icam).T(3,1));
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(1,1))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(1,2))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(1,3))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(2,1))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(2,2))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(2,3))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(3,1))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(3,2))');
    fprintf(fid, '%12.8f\t# Rinv\n', (camParaCalib(icam).Rinv(3,3))');
	fprintf(fid, '%15.8f\t# Tinv\n', camParaCalib(icam).Tinv(1,1));
    fprintf(fid, '%15.8f\t# Tinv\n', camParaCalib(icam).Tinv(2,1));
    fprintf(fid, '%15.8f\t# Tinv\n', camParaCalib(icam).Tinv(3,1));
    fprintf(fid, '\n');
	
end

fprintf(fid, '#### parameter for 3D matching ####');
fprintf(fid, '\n\n');
fprintf(fid, '0.01\t# mindist_pix (pixel)');
fprintf(fid, '\n');
fprintf(fid, '1\t# mindist_3D (mm)');

fclose(fid);