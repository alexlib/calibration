% restoredefaultpath;clear;close all

NAS = 'U';
addpath([NAS ':\skramel\codes\2-20-2017\calibration\']); % codes

% promptstr = strcat(promptstr, 'How many cameras in the PTV system?\n');
ncams = 4
gridsize = 2
% promptstr = 'Please indicate along which axis the mask is moving (x, y, or z)\n';
movingaxis = 'x'

if strcmp(lower(movingaxis), 'x')
    maskaxis=['y', 'z'];
end
if strcmp(lower(movingaxis), 'y')
    maskaxis=['x', 'z'];
end
if strcmp(lower(movingaxis), 'z')
    maskaxis=['x', 'y'];
end

%img_inversion = input('If the patterns are darker than background, you need to invert the image. \n Do you want to invert? (y/n)  ', 's');
img_inversion = 'y'
%ctr_finding = input('Which method of particle-center-finding do you want to use, Gaussianfit or COM? (g/c)  ', 's');
ctr_finding = 'gaussianfit';

nfig = figure;

icam = input('Please type camera-id of the camera you want to calibrate (1=camera 1, 2=camera 2,.., 4).\n');
promptstr = sprintf('\\nNow starting calibartion for camera %d of %d\\n', icam, ncams);
%promptstr = strcat(promptstr, 'How many calibration images for this camera?\n');
%nimgs = input(promptstr);
nimgs = 1
promptstr = sprintf('Please type the file name of the %s-axis positions of the mask.', movingaxis);
promptstr = strcat(promptstr, '\nNote that the order of the positions must match the order of the images to be processed.\n');
%posfname = input(promptstr, 's');
posfname='cal1_Apr10.txt';
axispos = myload(posfname, '#');
pointsfname = strcat('calibpointspos.cam',num2str(icam-1),'.dat');
if exist(pointsfname)
    promptstr = strcat('Default file for saving particle center position "', pointsfname, '" already exists. \n');
    promptstr = strcat(promptstr, 'Please type another file name. Type "o" for overwrite.\n');
    %response = input(promptstr, 's');
    response = 'o'
    if (length(response) > 1) | (lower(response(1)) ~= 'o')
        pointsfname = response;
    end
end

fidpos = fopen(pointsfname, 'w');
fprintf(fidpos, '# pimg_x (pix)\t pimg_y (pix) \t x (mm) \t y (mm) \t z (mm) \t Np \t sum(x) \t sum(y) \t sum(I) \n');
pimg = [];
pos3D = [];

str = sprintf('Please type the filename of calibration image (tiff image for Voth lab) %d\\n', nimgs);
fname = input(str, 's');
while ~exist(fname)
    str = strcat('File "', fname, '" does NOT exist!!! \n', str);
    fname = input(str, 's');
end
if lower(img_inversion(1)) == 'y'
    %Iimg = 255 - mean(double(mCINREAD(fname, 1, 1, 1)), 3);
    Iimg = 255 - imread(fname);
else
    %Iimg = mean(double(mCINREAD(fname, 1, 1, 1)), 3);
    Iimg = imread(fname);
end
[Npix_y Npix_x] = size(Iimg);
camParaknown.Npixh = Npix_y;
camParaknown.Npixw = Npix_x;
camParaknown.hpix = .014;    % pixel size (mm)   [basler A504k = .012, mikrotron mc1631 = .014]
camParaknown.wpix = .014;    % pixel size (mm)   [basler A504k = .012, mikrotron mc1631 = .014]
figure(nfig);
redo = 'y';
while strcmp(lower(redo), 'y')
    xc=[];
    yc=[];
    Ap=[];
    Ith = zeros(Npix_y, Npix_x);
%     h1 = subplot('position', [0.05 0.3 0.4 0.6]);
%     h2 = subplot('position', [0.05 0.08 0.4 0.15]);
%     h3 = subplot('position', [0.55 0.3 0.4 0.6]);
%     subplot(h1);
    figure(1)
    imagesc(uint8(Iimg));
    hold on;
    %	title(prnstr(fname));
    disp('Now please indicate the area of interest');
    nsubregion = input('How many sub-regions you want to divide the image to?  ');
    for isub = 1:nsubregion
        str = sprintf('Please choose region %d of %d using the mouse', isub, nsubregion);
        disp(str);
%         subplot(h1);
        figure(1)
        but = 0;
        while(but ~= 1)
            [xmin ymin but] = ginput(1);
        end
        xmin = max(floor(xmin), 1);
        ymin = max(floor(ymin), 1);
        but = 0;
        while(but ~= 1)
            [xmax ymax but] = ginput(1);
        end
        xmax = min(ceil(xmax), Npix_x);
        ymax = min(ceil(ymax), Npix_y);
        area = [xmin xmax ymin ymax];
        hold on
        plotrect(area, 'b--');
        hold off
        Nhist = hist(reshape(double(Iimg(ymin:ymax, xmin:xmax)), (xmax-xmin+1)*(ymax-ymin+1), 1), [0:255]);
%         subplot(h2);
        figure(2)
        semilogy([0:255], Nhist, 'b-');
        axis([0 255 1 10000]);
        th = input('Please choose threshold   ');
        Amin = input('Please choose the minimum particle image size   ');
        [xcsub ycsub Apsub Ithsub] = par_ctr(Iimg, th, Amin, ctr_finding, 'noshow', area);
        xc=[xc; xcsub];
        yc=[yc; ycsub];
        Ap = [Ap; Apsub];
        Ith = Ith+Ithsub; % Note Ith is a binary image
%         subplot(h3);
        figure(3)
        imagesc(Ith);
        hold on;
        plot(xc, yc, 'r+');
        hold off
    end
    redo = input('Do you want to re-process the image? (y/n)   ', 's');
end
Np = length(xc);

% manually remove some points if it's more convenient
rmman = input('Do you want to manually remove some possibly wrong particles? (y/n)  ', 's');
if (lower(rmman(1)) == 'y')
    ind = ones(Np,1);
    nrm = 0;
    disp('Please click the particle centers that you want to remove.');
    disp('Right click the mouse when you are done.');
    but = 1;
    hold on
    while but ~= 3
        [xrm yrm but] = ginput(1);
        if but == 1
            dist = (xc-xrm).^2+(yc-yrm).^2;
            [mindist irm] = min(dist);
            plot(xc(irm), yc(irm), 'ro');
            ind(irm) = 0;
            nrm = nrm+1;
        end
    end
    hold off
    xc = xc(logical(ind));
    yc = yc(logical(ind));
    Ap = Ap(logical(ind),:);
    Np = length(xc);
    disp(strcat(num2str(nrm), ' points have been removed. Please check the image again'));
    imagesc(Ith);
    hold on;
    plot(xc, yc, 'r+');
    hold off
end

% find three base points that defines a triad
disp('Now please indicate the three base point on the mask by click mouse on the thresholded image.');
disp('Starting with the first');
but = 0;
while but ~= 1
    [x0 y0 but] = ginput(1);
end
dist = (xc-x0).^2+(yc-y0).^2;
[mindist i0] = min(dist);
% subplot(h3);
figure(3)
hold on
plot(xc(i0), yc(i0), 'bo');
promptstr = sprintf('What are the indices of this point in %s and %s-dir?\\n', maskaxis(1), maskaxis(2));
temp = input(promptstr);
i0xind = temp(1);
i0yind = temp(2);

disp('Please indicate the second base point');
but = 0;
while but ~= 1
    [x1 y1 but] = ginput(1);
end
dist = (xc-x1).^2+(yc-y1).^2;
[mindist i1] = min(dist);
plot(xc(i1), yc(i1), 'bo');
promptstr = sprintf('What are the indices of this point in (ex. [0,0] %s and %s-dir?\\n', maskaxis(1), maskaxis(2));
temp = input(promptstr);
i1xind = temp(1);
i1yind = temp(2);

disp('Please indicate the third base point');
but = 0;
while but ~= 1
    [x2 y2 but] = ginput(1);
end
dist = (xc-x2).^2+(yc-y2).^2;
[mindist i2] = min(dist);
plot(xc(i2), yc(i2), 'bo');
% hold off
promptstr = sprintf('What are the indices of this point in %s and %s-dir?\\n', maskaxis(1), maskaxis(2));
temp = input(promptstr);
i2xind = temp(1);
i2yind = temp(2);

% Need to know a point close to the camera in order to
% compensate for the perspective projection distortion when find particle coordinates
disp('Please indicate a point that is nearly the closest to the camera');
but = 0;
while but ~= 1
    [xcam0 ycam0 but] = ginput(1);
end


% Now determine the point coordinates
% first, form two base vectors on the mask
e1 = [i1xind-i0xind, i1yind-i0yind];
e2 = [i2xind-i0xind, i2yind-i0yind];
% The prjection of these two vectors on image plane
e1p = [xc(i1)-xc(i0), yc(i1)-yc(i0)];
e2p = [xc(i2)-xc(i0), yc(i2)-yc(i0)];
e1pnorm = sum(e1p.^2);
e2pnorm = sum(e2p.^2);
e1pe2p = sum(e1p.*e2p);
d = (e1pnorm*e2pnorm - e1pe2p*e1pe2p);		% The denominator
% calaulte the coords of all points using the two base vectors
pind = zeros(Np, 2);
lenref = (xcam0-xc(i0))^2+(ycam0-yc(i0))^2;
for i=1:Np
    c = [xc(i)-xc(i0), yc(i)-yc(i0)];
    % The coefficient is used to compensate the effect of perspective projection
    coef = 1 + 0.0*(((xc(i)-xcam0)^2+(yc(i)-ycam0)^2) - sum(c.^2))/lenref;
    A = coef*(sum(c.*e1p)*e2pnorm - sum(c.*e2p)*e1pe2p)/d;
    B = coef*(sum(c.*e2p)*e1pnorm - sum(c.*e1p)*e1pe2p)/d;
    pind(i,:)=[A B];
end
% Now calculate the two components of dots' 3D coordinates on the mask plane
pmask = zeros(Np,2);
for i=1:Np
    pmask(i,:) = floor(e1*pind(i,1) +0.5) + floor(e2*pind(i,2)+0.5) + [i0xind i0yind];
end

% check to see if there is any inconsistency
ncoll = 0;
icoll = [];
for i=2:Np
    for j = 1:i-1
        if (pmask(i,1) == pmask(j,1)) & (pmask(i,2) == pmask(j,2))
            ncoll = ncoll+1;
            icoll = [icoll; [i j]];
        end
    end
end
if ncoll == 0
    disp('No confilcts found, but please still check particle coordinates.');
else
    str = sprintf('Some particle coordinates are probably wrong. %d conflicts found: \n', ncoll);
    for  ic = 1:ncoll
        str = strcat(str, sprintf('(%d, %d)\n', pmask(icoll(ic),:)));
    end
    disp(str);
end

change = 'y';
while lower(change) == 'y'
%     subplot(h3)
    figure(3)
    imagesc(Ith);
    hold on;
    plot(xc, yc, 'r+');
    plot(xc(i0), yc(i0), 'bo');
    plot(xc(i1), yc(i1), 'bo');
    plot(xc(i2), yc(i2), 'bo');
    for i=1:Np
        str = strcat('(',num2str(pmask(i,1)),',',num2str(pmask(i,2)),')');
        text(xc(i), yc(i), str, 'Color', 'g');
    end
    hold off
    change = input('Do you want to manually change coordinates? (y/n)  ', 's');
    if lower(change(1)) == 'y'
        disp('Please click the particle centers that you want to change manually.');
        disp('Right click the mouse when you are done.');
        but = 1;
        hold on
        while but ~= 3
            [xm ym but] = ginput(1);
            if but == 1
                dist = (xc-xm).^2+(yc-ym).^2;
                [mindist im] = min(dist);
                plot(xc(im), yc(im), 'ro');
                promptstr = sprintf('What are the NEW indices of this point in %s and %s-dir?\\n', maskaxis(1), maskaxis(2));
                temp = input(promptstr);
                pmask(im, :) = temp;
            end
        end
        hold off
        disp('The coords have been changed. Please check the image again');
    end
end

pmask = pmask*gridsize;
if strcmp(movingaxis, 'x')
    p3D = [axispos(nimgs)*ones(Np,1) pmask];
end
if strcmp(movingaxis, 'y')
    p3D = [pmask(:,1) axispos(nimgs)*ones(Np,1) pmask(:,2)];
end
if strcmp(movingaxis, 'z')
    p3D = [pmask axispos(nimgs)*ones(Np,1)];
end
% write the data to file
fprintf(fidpos, '# mask position %d, %d points\n', nimgs, Np);
fprintf(fidpos, '%12.7f\t%12.7f\t%12.7f\t%12.7f\t%12.7f\t%5d\t%12.5f\t%12.5f\t%12.5f\n', [xc yc p3D Ap]');
pimg = [pimg; [xc yc]];
pos3D = [pos3D; p3D];

%calibpointspos.cam.dat is not formatted properly, unknown effects-DB8,3]

fclose(fidpos);
%save PTVsetupWorkspace
% Now calibration
load camParaCalib-2-20-2017.mat;
camParaCalib(icam) = db_calib_Tsai_Co_planar(pimg, pos3D, camParaknown);

save camParaCalib-2-20-2017-fixed.mat camParaCalib