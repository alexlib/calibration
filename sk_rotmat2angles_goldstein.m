function [eul1,eul2] = sk_rotmat2angles_goldstein(o)
% Convert a rotation matrix to Euler angles using Golstein convention (p153).
% Account for the Gimbal lock case and return two both possible euler angles that correspons to the rotation matrix. 

% amount residual at euler angle pole
rp=abs(o(3,3)-1);
rm=abs(o(3,3)+1);

% Boolean saying if the euler angle is too close to pole
boolp = rp<1e-11;
boolm = rm<1e-11;

if ( boolp || boolm )
    ps1 = 0; % Gimbal lock, psi is undefined
    if boolp
        th1 = rp;
        ph1 = atan2(o(1,2),o(1,1));
    else
        th1 = pi-rm; % same as -pi in our case
        ph1 = atan2(o(1,2),o(1,1));
    end
    
    % There is only one set of angles in this case
    ph2=ph1;
    th2=th1;
    ps2=ps1;
    
else
    % There are two possible euler representations except in the case
    % handled above.
    
    % first set of angles
    th1 = acos(o(3,3));
    ph1 = atan2(o(3,1)/sin(th1),-o(3,2)/sin(th1));
    ps1 = atan2(o(1,3)/sin(th1),o(2,3)/sin(th1));
    
    % second set of angles
    th2 = -th1;
    ph2 = atan2(o(3,1)/sin(th2),-o(3,2)/sin(th2));
    ps2 = atan2(o(1,3)/sin(th2),o(2,3)/sin(th2));
    
end

    % Return the angles
    eul1 = [ph1,th1,ps1];
    eul2 = [ph2,th2,ps2];

end