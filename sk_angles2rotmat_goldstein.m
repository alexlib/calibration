function ori = sk_angles2rotmat_goldstein(eul)
% This function uses the Goldstein convention to turn an array containing
% the euler angles into a rotation/orientation matrix and returns it.

phi = eul(1);
th = eul(2);
psi = eul(3);

% Alternatively, could implement the D, C, B matrix convention. Might be
% easier to follow. 
o1 = [cos(psi)*cos(phi) - cos(th)*sin(psi)*sin(phi); ...
    -sin(psi)*cos(phi) - cos(th)*cos(psi)*sin(phi); ...
    sin(th)*sin(phi)];

o2 = [cos(psi)*sin(phi) + cos(th)*cos(phi)*sin(psi); ...
    -sin(psi)*sin(phi) + cos(th)*cos(phi)*cos(psi); ...
    -sin(th)*cos(phi)];

o3 = [sin(psi)*sin(th); ...
    cos(psi)*sin(th); ...
    cos(th)];
    
% Now we assemble the entire rotation matrix
ori = [o1,o2,o3];


% if round(det(ori)) == -1
%     % check for and fix left-handed coords.
%     ori(:,1)=-ori(:,1);
% end

end