clear;close all;restoredefaultpath;

NAS3_user = 'U';
NAS3_data = 'V';

addpath(genpath([NAS3_user ':\skramel\codes\'])); % codes
addpath(genpath([NAS3_data ':\skramel\2-20-2017\'])); % calibration

% define paths
datapath = [NAS3_data ':\skramel\2-20-2017\'];
calipath = [NAS3_user ':\skramel\codes\2-20-2017\calibration\'];
filepath = [NAS3_data ':\skramel\2-20-2017\results\'];
savepath = [NAS3_user ':\skramel\codes\2-20-2017\calibration\'];

% define files (no file extensions)
filename = 'dynamic_calib_6';
caliname = 'dynamic_camParaCalib-2-20-2017';

camParaCalib = load([calipath caliname]);
camParaCalib = camParaCalib.camParaCalib;

tracked=read2_gdf([filepath filename '\' filename '_tracers_NEW_tracked_4cams_3.gdf']);
tracked(tracked(:,19)==1,:)=[];
% tracked(tracked(:,8)==0,:)=[];
% tracked(tracked(:,11)==0,:)=[];
% tracked(tracked(:,14)==0,:)=[];
% tracked(tracked(:,17)==0,:)=[];

clear pos3d2d
pos3d2d(:,1:3)=tracked(:,2:4);
pos3d2d(:,4)=sqrt(tracked(:,18));
pos3d2d(:,5)=99;
pos3d2d(:,6)=tracked(:,5).*450;
pos3d2d(:,7:14)=tracked(:,[6:7 9:10 12:13 15:16]);

for col=7:14
    pos3d2d(isnan(pos3d2d(:,col)),:)=[];
    pos3d2d(isinf(pos3d2d(:,col)),:)=[];
    pos3d2d(pos3d2d(:,col)<=0,:)=[];
    pos3d2d(pos3d2d(:,col)>=1280,:)=[];
end

ncams=size(camParaCalib,2);

%for our data with 200mm lenses and 1cm field of view, the radial
%distortion does not seem to affect the calibration quality at all, so it
%is set to zero to simplify things.  This code should work fine
%with these four lines commented out.  In the current form, it would leave
%the distortion as a fixed parameter, but it could be optimized as well.
for icam=1:ncams
    camParaCalib(icam).k1=0;
    camParaCalib(icam).k1star=1;
end
%Now we pick which 3d points to use.  It's important to pick 3d points that
%are found using the cameras we are calibrating.

%Pick the points to use for the fit.  Because the data is tracks, the
%positions are highly correlated and it is best to choose positions spaced
%throughout the file.
npointsfit=5000; %4e3 minimum number of points required for a fit
pc=[1 2 3 4];
pp=(1:size(pos3d2d,1))';     %pertinant particles
ncams=length(pc);
[n2d,ncols]=size(pp);
if n2d < npointsfit
    warning('not enough points to calibrate')
end
step = fix(n2d/npointsfit);
chosen=pp(1:step:n2d);
pos3d2d_chosen=pos3d2d(chosen(1:npointsfit),:);  %Not very elegant--it just needs to select npointsfit elements from throughout dist3D

%Plot the histogram of mismatches as determined by the matching code that
%created the pos3d2d file.
figure(1)
hist(pos3d2d(pp,4),80);
title('Matching error in pos3d2d data');
xlabel('error (mm)');

cam2d=zeros(npointsfit,2,ncams);

%Choose another set of points to check the calibration against at the end.
chosen_check=pp(10:step:n2d);
%chosen_check=(round(setdiff(1:n2d,chosen)))';
pos3d2d_chosen_check=pos3d2d(chosen_check(1:npointsfit),:);
cam2d_check=zeros(npointsfit,2,ncams);

%repack the 2d camera coordinate data into a cam2d array.
for icam = 1:ncams
    cam2d(:,:,icam)=pos3d2d_chosen(:,(pc(icam)*2)+5:(pc(icam)*2)+6);
    cam2d_check(:,:,icam)=pos3d2d_chosen_check(:,(pc(icam)*2)+5:(pc(icam)*2)+6);
end

%show the particles from camera 1
figure(2)
plot(cam2d(:,1,1),cam2d(:,2,1), '.g');
title('Points used as seen by the first camera')
xlabel('x position (pixels)')
ylabel('y position (pixels)')

for icam = 1:ncams
    cam2d(:,1,icam)=(cam2d(:,1,icam)-camParaCalib(pc(icam)).Npixw/2-camParaCalib(pc(icam)).Noffw)*camParaCalib(pc(icam)).wpix;
    cam2d(:,2,icam)=(-cam2d(:,2,icam)+camParaCalib(pc(icam)).Npixh/2-camParaCalib(pc(icam)).Noffh)*camParaCalib(pc(icam)).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
    cam2d_check(:,1,icam)= (cam2d_check(:,1,icam)-camParaCalib(pc(icam)).Npixw/2-camParaCalib(pc(icam)).Noffw)*camParaCalib(pc(icam)).wpix;
    cam2d_check(:,2,icam)=(-cam2d_check(:,2,icam)+camParaCalib(pc(icam)).Npixh/2-camParaCalib(pc(icam)).Noffh)*camParaCalib(pc(icam)).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
end

%Now stuff the calibration data into the arrays needed by fminsearch.
%calinitial and calarray contain the calib parameters for all 3 cameras.
%fminsearch wants an array with only the optimization parameters, so we pass
%those separately, and then stuff the optimizing parameters back into the full
%calarray inside gv_dynamic_fitfunc.  If you change which parameters are
%optimized you need to change gv_dynamic_fitfunc so that it stuffs the
%correct array elements.

calinitial=zeros(8,ncams);
for icam = 1:ncams
    calinitial(1:3,icam)=sk_rotmat2angles_goldstein(camParaCalib(pc(icam)).R);
    calinitial(4:6,icam)=camParaCalib(pc(icam)).T;
    calinitial(7,icam)=camParaCalib(pc(icam)).f_eff;
    calinitial(8,icam)=camParaCalib(pc(icam)).k1;
end
calarray=calinitial;

%plot initial distribution of mismatches--not necessary, but
%convenient to do here.  Should be similar to figure 1
[dist3, dist1,points3D]=sk_calc_ray_mismatch(calarray,cam2d,ncams);
figure(3);
hist(dist3,50);
title('initial mismatch distribution using initial calibration');
xlabel('mismatch (mm)')
allcams_initial_mismatches=mean(dist3);
display(strcat('before anything, the mean mismatch of the previous found 3d points is:',num2str(mean(pos3d2d_chosen(:,4)))))

display(strcat('before any calibration adjusting the mean mismatch is:',num2str(allcams_initial_mismatches)))

display('the first 10 chosen 3d particles are:')
pos3d2d_chosen(1:10,1:4)
display('after converting Rotation matrix to angle and back again the 3d points should be the same')
points3D(1:10,:)

%The parameters not in calopt are fixed.  Currently we only fit R and T,
%the rotation matrix and the position of the cameras (the rotation matrix
%is represented by 3 angles).
calopt=zeros(6,ncams-1);
for icam=2:ncams            %The first camera is fixed, the remaining two move.
    calopt(1:3,icam-1)=sk_rotmat2angles_goldstein(camParaCalib(pc(icam)).R);
    calopt(4:6,icam-1)=camParaCalib(pc(icam)).T;
    %calopt(7, 1)=camParaCalib(icam).f_eff;
    %calopt(8,1)=camParaCalib(icam).k1;
end

%gv_dynamic_fitfunc(calopt,calfixed,cam2d)
fmin_options.Display='final'; %'final' to display only the final mean distance, 'iter' to show the steps along the way.
fmin_options.MaxFunEvals=2000 %350;
calout = fminsearch(@(x) gv_dynamic_fitfunc(x,calarray,cam2d,ncams), calopt, fmin_options);

calarray(1:6,2:ncams)=calout(1:6,1:ncams-1);
[dist3, dist1]=sk_calc_ray_mismatch(calarray,cam2d,ncams);

%display results of initial optimization
figure(4);
hist(dist3,50);
title('mismatch distribution after first optimization');
xlabel('mismatch (mm)')
allcams=mean(dist3)
cammismatch1=zeros(ncams);
for n1=1:ncams
cammismatch1(n1)=mean(dist1(n1,:));
end


%choose good matches
mingoodmatch=mean(dist3)*0.75;%.2;%.04  %for now let's take the best 75% of matches
good_matches=(dist3 < mingoodmatch);
cam2d_good=cam2d(good_matches,:,:);


calopt=calarray(1:6,2:ncams);
%By using calopt=calarray(1:8,2:3) here (and using gv_tmp_dynamic_fitfunc) this
%can optimize all 8 parameters, but it doesn't give any better fit, so I
%am using the simpler 12 parameter fit as above but without the bad matches.

calout = fminsearch(@(x) sk_dynamic_fitfunc(x,calarray,cam2d_good,ncams), calopt, fmin_options);
calarray(1:6,2:ncams)=calout(1:6,1:ncams-1);

%Check the quality of the final calibration
[dist3, dist1]=sk_calc_ray_mismatch(calarray,cam2d_good,ncams);

figure(5);
hist(dist3,50);
title('mismatch distribution after optimization on good matches');
xlabel('mismatch (mm)')
allcams=mean(dist3)
cammismatch2=zeros(ncams);
for n1=1:ncams
cammismatch2(n1)=mean(dist1(n1,:));
end

[dist3, dist1]=sk_calc_ray_mismatch(calarray,cam2d_check,ncams);

figure(6);
hist(dist3,50);
title('final mismatch distribution of data not used for the optimization');
xlabel('mismatch (mm)')
allcams_check_with_mismatches=mean(dist3)


% Refill the camParaCalib structure and write it to a .cfg file.
%  (There turns out to be a very small change in the first camera even though it was not
%    dynamically calibrated since the rotation matrix is projected onto a matrix with
%     determinant exactly -1, so it changes just a bit.  The new R needs to be kept since the
%     dynamic calibration was done for this rotation matrix)
for icam = 1:ncams
    camParaCalib(pc(icam)).R=sk_angles2rotmat_goldstein(calarray(1:3,icam));
    camParaCalib(pc(icam)).Rinv=inv(camParaCalib(icam).R);
    camParaCalib(pc(icam)).T=calarray(4:6,icam);
    camParaCalib(pc(icam)).Tinv=camParaCalib(icam).Rinv * (-1* camParaCalib(icam).T);
    camParaCalib(pc(icam)).f_eff = calarray(7,icam);
    camParaCalib(pc(icam)).k1 = calarray(8,icam);
    camParaCalib(pc(icam)).err_x=cammismatch2(icam);
    camParaCalib(pc(icam)).err_y=cammismatch2(icam);
end

beep
beep 
beep
%fname_cfg=strcat(dist3D_path,'PTVSetup_optimized_',dist3D_stem,'.cfg');
%gv_write_calib_cfg(camParaCalib, ncams, fname_cfg);


save([savepath 'dynamic2_camParaCalib-2-20-2017'],'camParaCalib');

