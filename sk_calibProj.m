function X_proj = sk_calibProj(camParaCalib, X3D)
% Use the calibrated camera parameters to predict the particle position
% projected onto the image plane.
%
% inputs:
%   camParaCalib    --  calibrated camera parameters
%   X3D         --  test particle coordinates in 3D world system
%
% output:
%   X_proj      --  projected particle position on image plane (in pixels)
%
Xc = X3D * (camParaCalib.R)';
Xc(:,1) = Xc(:,1) + camParaCalib.T(1);
Xc(:,2) = Xc(:,2) + camParaCalib.T(2);
Xc(:,3) = Xc(:,3) + camParaCalib.T(3);
dummy = camParaCalib.f_eff./Xc(:,3);
Xu = Xc(:,1).*dummy;  % undistorted image coordinates
Yu = Xc(:,2).*dummy;
ru2 = Xu.*Xu + Yu.*Yu;
dummy = 1+camParaCalib.k1*ru2;
Xd = Xu.*dummy;
Yd = Yu.*dummy;
% iterate once
dummy = 1 + camParaCalib.k1*(Xd.*Xd + Yd.*Yd);
Xd = Xu.*dummy;
Yd = Yu.*dummy;
Np = size(X3D,1);
X_proj = zeros(Np,2);
X_proj(:,1) = Xd/camParaCalib.wpix + camParaCalib.Noffw + camParaCalib.Npixw/2;
X_proj(:,2) = camParaCalib.Npixh/2 - camParaCalib.Noffh - Yd/camParaCalib.hpix;

