%evaluate experimental molecule Data in biplane imaging mode --> 3D localisations
%data from two image channels ("top" and "bottom") are used 
%one is for instance defocused, the other in focus
%note that in this case the defocused channel must be defined as "BOTTOM"
%image!

clc
clear all; 
close all

global TS pxSize dist_max data_BG PSF_top_norm PSF_bottom_norm interp_incr x_data z_vec ux
%% user parameters 

pxSize=117;
dist_max=400;

%camera parameters
gain=100; 
amp=9.9; %electrons per count 
QE=0.95; %quantum efficiency
baseline=100; %offset level in ADUs; should be 100 for Andor ixon cams

r_sphere=2.0e-3/2; %radius of sphere
ratio=1/1.4; %fixed energy-ratio between both images ratio=energy_topimage/energy_bottomimage

PSFpath='C:\Users\q004aj\Desktop\PSFs\';

%PSF for top image stack: 
name_PSF_top='PSF5D_13x13_0-2-250nm_RI=1.45_dz=0_aberrfree_os3.mat';
%PSF for bottom image stack: 
name_PSF_bottom='PSF5D_13x13_0-2-250nm_RI=1.45_dz=-500_aberrfree_os3.mat';


%% loading PSF models

% ----PSF for top imgae--------

load([PSFpath name_PSF_top]); PSF_top=PSF5D;

if ndims(PSF_top)==3 %if "standard" 3D PSF is loaded (interpolation in MLE necessary -> slower)
    [nx0,ny0,nz0]=size(PSF_top); %size of model    clear PSF_tot;
    %PSF: normalization of energy in each z-slice
    energies_top=(trapz(trapz(PSF_top,1),2)); 
    %PSF_norm_SAF=PSF_SAF./repmat(energies_SAF,[size(PSF_SAF,1),size(PSF_SAF,2),1]);
    interp_incr=1; 
elseif ndims(PSF_top)==5 
    [nx0,ny0,nz0,nxi,nyi]=size(PSF_top); %size of model
    clear PSF5D;
    %PSF: normalization of energy in each z-slice
    energies_top=(trapz(trapz(PSF_top(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)),1),2)); 
    %PSF_norm_SAF=PSF_SAF./repmat(energies,[size(PSF_SAF,1),size(PSF_SAF,2),1,nxi,nyi]);
end
disp('PSF for top image loaded');

if exist('sigma')
    z_ini_info1=sigma/os; %if the loaded PSF-model contains a vector called "sigma"; sigma represents a gaussian-width-versus-z_ini-curve which provides a good initial z-estimate
else  %otherwise, take the center-value of the entire z-range as initial estimate for z
    z_ini_info1=round(size(PSF,3)/2); %constant initial z-estimate (frame-no of PSF-stack)
end

% -----loading PSF model for bottom image -----

load([PSFpath name_PSF_bottom]); PSF_bottom=PSF5D; 

if ndims(PSF_bottom)==3 %if "standard" 3D PSF is loaded (interpolation in MLE necessary -> slower)
    [nx0,ny0,nz0]=size(PSF_bottom); %size of model    clear PSF_tot;
    %PSF: normalization of energy in each z-slice
    energies_bottom=(trapz(trapz(PSF_bottom,1),2)); 
    %PSF_norm_UAF=PSF_UAF./repmat(energies,[size(PSF_UAF,1),size(PSF_UAF,2),1]);
    interp_incr=1; 
elseif ndims(PSF_bottom)==5 
    [nx0,ny0,nz0,nxi,nyi]=size(PSF_bottom); %size of model
    clear PSF5D;
    %PSF: normalization of energy in each z-slice
    energies_bottom=(trapz(trapz(PSF_bottom(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)),1),2)); 
    %PSF_norm_UAF=PSF_UAF./repmat(energies_UAF,[size(PSF_UAF,1),size(PSF_UAF,2),1,nxi,nyi]);
end
disp('PSF for bottom image loaded');

if exist('sigma')
    z_ini_info2=sigma/os; %if the loaded PSF-model contains a vector called "sigma"; sigma represents a gaussian-width-versus-z_ini-curve which provides a good initial z-estimate
else  %otherwise, take the center-value of the entire z-range as initial estimate for z
    z_ini_info2=round(size(PSF,3)/2); %constant initial z-estimate (frame-no of PSF-stack)
end

ux=ux*os; %correcting for oversampling factor

%----------------- combine models---------------

fw=0; %frame-width (width of border around the PSF-images, needed to allow for x-y-localization (shifts))
       
%create coordinate system 
x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
[X,Y]=ndgrid(x,y);
clear x_data;
x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx, ny]=size(X); %size of molecule image

uz=z_vec(2)-z_vec(1); %z-increment provided with PSF-models
PSF_tot=[squeeze(PSF_top(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2))) squeeze(PSF_bottom(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)))];


%% plot CRB for combined models user defined signal/bg-level

% signal=10e3; 
% background=500; 
% [CRBx,CRBy,CRBz]=fun_CRB(PSF_tot./repmat(trapz(trapz(PSF_tot,1),2),[size(PSF_tot,1) size(PSF_tot,2) 1]),ux,uz,signal,background,gain);
% 
% figure(1)
% plot(z_vec(1:end-1)*1e9,sqrt(CRBz));
% title(['mean \sigma_z: ' num2str(mean(sqrt(CRBz)),3) 'nm for sig/bg=' num2str(signal) '/' num2str(background)]);
% xlabel('z-pos/nm');
% ylabel('CRB /nm');
% grid on;

%% reading in experimental data (movie of blinking molecules)
% top and bottom stacks must be loaded separately; you can create them in ImageJ by
% selecting corresponding ROIs in the whole recorded stack (containing top &
% bottom parts) and pressing:  ctrl+shift+D

clear data_top data_bottom
[Name,Path,~]=uigetfile('*.tif','choose TOP and BOTTOM tif stacks','.\','MultiSelect','on');
    image_no=length(imfinfo([Path Name{1}]));
    for m=1:image_no
      disp(image_no-m);
      data_bottom(:,:,m)=flipud(double(imread([Path Name{1}],m)))-baseline; %data is instantly flipped to compensate for flip on camera
      data_top(:,:,m)=(double(imread([Path Name{2}],m)))-baseline; %weighting image energy with "ratio"
    end
disp('done');

% % EVER - background reduction
% % find minima in data stack
% V_min_top=min(data_top,[],3);
% V_min_bottom=min(data_bottom,[],3);
% V_BG1_top=(V_min_top+27.8)/0.9513 ; 
% V_BG1_bottom=(V_min_bottom+27.8)/0.9513 ; 
% 
% imagesc([V_BG1_top; V_BG1_bottom]); title('background');
% 
% data_corr_top=data_top-repmat(V_BG1_top,1,1,size(data_top,3));
% data_BG_top=V_BG1_top;
% 
% data_corr_bottom=data_bottom-repmat(V_BG1_bottom,1,1,size(data_bottom,3));
% data_BG_bottom=V_BG1_bottom;

%% creating coordinate transform between TOP and BOTTOM stacks 
%manual selection of point-pairs 

if isfile('transformpoints.mat') %if a coord. transform file is already contained in the current folder it is loaded
    load transformpoints.mat
    disp('transform matrix found and loaded.');
else  %otherwise, create a new transform matrix by clicking on matching points in the top and bottom images
    
    prompt='no transform matrix found. Create one by clicking on matched point pairs. Do you have specific calib. images? (e.g. containing beads) (y/n)';
    answer = inputdlg(prompt);
    
    if strcmp(answer{1},'n')
        img_frame=1; %frame to compare points
        calib_image1=data_top(:,:,img_frame);
        calib_image2=data_bottom(:,:,img_frame);
    else
       [Name,Path,~]=uigetfile('*.tif','choose TOP and BOTTOM calibration images','.\','MultiSelect','on');
       calib_image1=imread([Path Name{1}]);
       calib_image2=imread([Path Name{2}]);
    end
    
    cpselect(calib_image1/max(calib_image1(:)),calib_image2/max(calib_image2(:))); %save points as "movingPoints", "fixedPoints";
    movingPoints= cpcorr(movingPoints, fixedPoints, calib_image1, calib_image2); %fine-tuning using cross-correlation
    %------saving points to harddrive-execute the following line--------: 
    %save('transformpoints.mat','movingPoints','fixedPoints'); disp('transform matrix saved. ');
    %--------------------------------------------------------------------
end
   

%% load Thunderstorm locations found in the top image and calculating corresponding locations in the bottom image
%-> check if molecule images are accurately cropped from the image stacks

tmp= csvread('TS_top.csv',1,0);
TS=tmp(1:5000,:);%select parts of the entire data
%TS=tmp; 

%----data pre-filtering options
thresh_dist=2000; %minimum distance of a molecule to the next in nm (others are deleted in the function "RmNN")
thresh_sigma=1*max(TS(:,5)); %sigma of Gaussian fit in nm
thresh_uncert=1*max(TS(:,10)); %
TS_filt = RmNN(TS(end,2),TS,thresh_dist,thresh_sigma,thresh_uncert); %picking only those molecule images with no neighbours in the vicinity


%---------------calculating coords of locations in the bottom image---------

img_no=1; %choose test image from stack

tform = fitgeotrans(movingPoints,fixedPoints,'affine'); %finding transformation from selected points
%tform = fitgeotrans(fixedPoints,moving_pts_adj,'affine'); %finding transformation from selected points
%[tmp,RB]=imwarp(I1,tform);

TS_t=TS_filt; %initializing Thunderstorm coordinates for transformed (second) image

[xT,yT] = transformPointsForward(tform,TS_filt(:,3)/(ux*1e9),TS_filt(:,4)/(ux*1e9));
%[xT,yT] = transformPointsInverse(tform,TS_filt(:,3)/(ux*1e9),TS_filt(:,4)/(ux*1e9));
%[TS_t(:,3),TS_t(:,4)]=worldToIntrinsic(RB,xT,yT);

% ------test if correctly transformed-------------
% [xT,yT] = transformPointsForward(tform,moving_pts_adj(:,1),moving_pts_adj(:,2));
% plot(xT,yT,'b.'); hold on; 
% plot(fixedPoints(:,1),fixedPoints(:,2),'rx');
% hold off; 
% -------------------------------------------------

%<<<<<<<fine-adjust centering of bottom image by varying offset-values
offset=[0 0]*ux*1e9; %visually determined offset for crop-coordinates
TS_t(:,3:4)=[xT,yT]*(ux*1e9)-repmat(offset,size(TS_filt,1),1);

% %check transform by plotting overlaid points
% figure(1);
% plot(TS_filt(:,3)/(ux*1e9),TS_filt(:,4)/(ux*1e9),'b.')
% hold on;
% plot(TS_t(:,3),TS_t(:,4),'r.')
% hold off; 
% title('do the points match?');

%----check transform by plotting overlaid points----
idx=TS_filt(:,2)==img_no;
figure(1); 
subplot(2,1,1);
imagesc(data_top(:,:,img_no)); hold on; title(['top image no.' num2str(img_no) ' with overlaid TS-coords.']);
plot(TS_filt(idx,3)/(ux*1e9),TS_filt(idx,4)/(ux*1e9),'r.'); 
hold off;
subplot(2,1,2);
imagesc(data_bottom(:,:,img_no)); hold on; title(['bottom image no.' num2str(img_no) ' with transformed TS-coords.']);
plot(TS_t(idx,3)/(ux*1e9),TS_t(idx,4)/(ux*1e9),'r.'); 
hold off;

data_BG=zeros(size(data_top));
I_top=fun_thunderstorm_crop(TS_filt,data_top,data_BG,(nx-1)/2,ux); %crops molecule images in the top-image
I_bottom=fun_thunderstorm_crop(TS_t,data_bottom,data_BG,(nx-1)/2,ux); %crops molecule images in the bottom-image


%delting locs where I_top or I_bottom is zero and deleting them (can happen if a
%rotation is involved and the loc is close to the image border)
idx=find((squeeze(sum(sum(I_top,1),2))==0)|(squeeze(sum(sum(I_bottom,1),2))==0)); 
TS_filt(idx,:)=[];
TS_t(idx,:)=[];
I_top(:,:,idx)=[];
I_bottom(:,:,idx)=[];

%% <<<<<< joint MLE on the combined images >>>>>>>>>>>>
 
%correctly normalizing PSFs: each z-slice of the JOINT PSF is normalized
energies=repmat((trapz(trapz(PSF_tot,1),2)),[size(PSF_tot,1),size(PSF_tot,2)/2,1,nxi,nyi]); %energies of joint PSFs for each z-slice
PSF_top_norm=PSF_top;%./energies; 
PSF_bottom_norm=PSF_bottom;%./energies;

px=TS_filt(:,3); %approx. position of molecules in nm
py=TS_filt(:,4);

qual_thresh=inf; 

close all; 
clear x1_est y1_est z_est N_est BG1_est BG2_est resnorm_MLE loc_no frame_no
idx_disc=[];

m=0;
for mm=1:size(TS_filt,1)

    if mod(mm,10)==0 %show image every 20th frame or so
        showimage=1; 
    else
        showimage=1;
    end
    z_ini=[z_ini_info1; z_ini_info2];
    %z_ini=60;
    
    %---dual-channel evaluation----
    est=fun_MLE_2channels(I_top(:,:,mm), I_bottom(:,:,mm), PSF_top_norm, PSF_bottom_norm, interp_incr,x_data,z_ini,qual_thresh,ratio,showimage);
    %format: est=[x1,y1,x2,y2,z,Sig,BG1,BG2,resnorm1,resnorm2]
    
    %----single-channel eval-----------
    %tmp=fun_MLE(I_bottom(:,:,mm),PSF_bottom_norm,interp_incr,x_data,z_ini_info2,qual_thresh,showimage);
    %if isempty(tmp)==0
    %   est=[0 0 tmp(1:4) 0 tmp(5) 0 tmp(6)];
    %end
    
    if isempty(est) %if resnorm is too large, function returns est=[]
        idx_disc=[idx_disc mm]; %indidces of images to be discarded
    else
        m=m+1;
        
        x1_est(m)=(round(py(mm)/ux/1e9)+est(1))*ux*1e9; %note the round-operation! this is required because of the way the images are cropped in "fun_thunderstorm_crop"
        y1_est(m)=(round(px(mm)/ux/1e9)+est(2))*ux*1e9; %px measures along the column-direction
        z_est(m)=interp1(1:length(z_vec),z_vec,est(5))*1e9; %estimated z-position in nm
        N_est(m)=est(6)/gain*amp/QE; %N_est is the TOTAL signal in both images
        BG1_est(m)=est(7)/gain*amp/QE;
        BG2_est(m)=est(8)/gain*amp/QE;
        resnorm_MLE(m)=est(9)+est(10);
        loc_no(m)=TS_filt(mm,1); %number of localization (temporal order)
        frame_no(m)=TS_filt(mm,2); %number of image frame
        disp(m);
%         disp(x1_est(m));
%         disp(N_est(m));
%         disp(BG1_est(m));
%         disp(z_est(m)); 
       
    end
end
locdata_TS=[x1_est' y1_est' z_est' N_est' BG1_est' BG2_est' resnorm_MLE' loc_no' frame_no'];
locdata_TS(isnan(locdata_TS))=0;

disp('all data localized');
       
%% saving localization results

NameSave='biplane_alexa_sphere_2019-03-26';
save(['LOCDATA_' NameSave '.mat'],'locdata_TS','ux','NA','gain','amp','QE','z_vec','r_sphere','I_top','I_bottom');

%% -----data eval----- 

finaldata=locdata_TS; 
I_sel=I_bottom; %selected images

%----------deleting erroneous entries / filtering data -----------
    r_cutoff=15e-6; %cutoff distance from sphere center in m
%manually defining sphere center
    if exist('data_top')
        xc=[1 1]*size(data_top,1)*ux*1e9/2; %center of image (in nm)
    else
        xc=[1 1]*(max(finaldata(:,1))-min(finaldata(:,1)))/2; %center of image (in nm)
    end
    findcenter='y'; %automatic center-finding? set to 'y' or 'n'
%select temporal junks of data to see time trends
    no_c=1000*2; %central local. no.
    no_delta=1000;  %half width, set inf to get all data
%selecting a defined azimuthal section of the sphere 
    phi_c=200; %central angle of cake-slice (in deg.)
    dphi=inf; %angular width of "cake-slice"
%delting entries with too high fit-errors (residuals)
    qual_thresh=15;  %the lower the cutoff, the more strict the selection
%filtering data: photon threshold
    photon_loBND=0;
    photon_hiBND=inf;
%----------------------------------------

%format: finaldata=[x1_est' y1_est' z_est' N_est' BG1_est' BG2_est' resnorm_MLE' loc_no' frame_no'];
%delete molecules with zero signal photons
idx=finaldata(:,4)==0;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%delete entries that run into upper z-boundaries
idx=finaldata(:,3)>=max(z_vec)*1e9-1;
finaldata(idx,:)=[]; 
I_sel(:,:,idx)=[]; 

%delete entries that run into lower z-boundaries
idx=finaldata(:,3)<=min(z_vec)*1e9+1;
finaldata(idx,:)=[]; 
I_sel(:,:,idx)=[]; 

%selecting temporal sections
idx=abs(finaldata(:,8)-no_c)>no_delta;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%delete entries which are too far away from the sphere center
idx=((finaldata(:,1)-xc(1)).^2+(finaldata(:,2)-xc(2)).^2)>(r_cutoff*1e9)^2;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%selecting a defined azimuthal section of the sphere 
y_tmp=finaldata(:,2)-xc(2); 
x_tmp=finaldata(:,1)-xc(1); 
phi=atan2(y_tmp,x_tmp)/pi*180; %in deg
idx=abs(mod(phi-phi_c+180,360)-180)>dphi;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%delting entries with too high fit-errors (residuals)
idx=finaldata(:,7)>qual_thresh*min(finaldata(:,7));
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%filtering data: photon threshold
idx=(photon_hiBND<finaldata(:,4));
finaldata(idx,:)=[];
idx=(finaldata(:,4)<photon_loBND);
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

filtdata=finaldata;

%---automatically finding optimal sphere-center: 
r_coord=linspace(0,30,size(PSF_top,3))*1e-6; %radial coord. in m
z_theory=r_sphere-sqrt(r_sphere^2-(r_coord).^2);
z_error=@(x) std(interp1(r_coord,z_theory,sqrt((filtdata(:,1)-x(1)).^2+(filtdata(:,2)-x(2)).^2)*1e-9,'linear','extrap')*1e9-filtdata(:,3));

if strcmp(findcenter,'y')
    xc0=xc;
    if exist('data_top')
        xc=[1 1]*size(data_top,1)*ux*1e9/2; %center of image (in nm)
    else
        xc=[1 1]*(max(finaldata(:,1))-min(finaldata(:,1)))/2; %center of image (in nm)
    end
    xc=fminsearch(z_error,xc0);
end
    %z_target=interp1(r_coord,z_theory,r_mol)*1e9; %traget z-value in nm

delta_z=z_error(xc); %standard deviation in nm 
r_mol=sqrt((filtdata(:,1)-(xc(1))).^2+(filtdata(:,2)-(xc(2))).^2)*1e-9; %radial molecule coords. in meter



%----3D scatterplot: show all filtered data----
figure(1); 
%plot3(locdata_TS(:,1)/1e3,locdata_TS(:,2)/1e3,locdata_TS(:,3),'o'); grid on;
markersize=3; %(locdata_TS(:,4)+1)/1e3;
markercolor=filtdata(:,3); 
scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,filtdata(:,3),markersize,markercolor); grid on;
colormap jet; 
colorbar; 
title('local. data, thunderstorm');
zlabel('z / nm');
xlabel('x / µm');
ylabel('y / µm');


% showing estimated photon numbers
figure(4);
subplot(3,1,1); 
hist(filtdata(:,4),20);
grid on;
title(['signal, mean=' num2str(mean(filtdata(:,4)))]);
xlabel('photons');
ylabel('number of localizations');
subplot(3,1,2);
hist(filtdata(:,5),20);
grid on;
title(['background top image, mean=' num2str(mean(filtdata(:,5)),2)]);
xlabel('photons');
ylabel('number of localizations');
subplot(3,1,3);
hist(filtdata(:,6),20);
grid on;
title(['background bottom image, mean=' num2str(mean(filtdata(:,6)),2)]);
xlabel('photons');
ylabel('number of localizations');


%----2D scatterplot with color-coded z-value----
figure(6); 
markercolor=filtdata(:,3);% markercolor=uint8(markercolor);
%markercolor=[0 0 0];
scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,2,markercolor);
xlabel('x / µm');
ylabel('y / µm');
axis equal; axis tight; 
colormap jet; 
colorbar; 
title('loc. data, z-values color-coded');

% -----plotting radial distance versus estimated z-positions-----

mean_signal=mean(filtdata(:,4));
mean_bg=mean(filtdata(:,5));
[CRBx,CRBy,CRBz]=fun_CRB(PSF_tot./repmat(trapz(trapz(PSF_tot,1),2),[size(PSF_tot,1) size(PSF_tot,2) 1]),ux,uz,mean_signal,mean_bg,gain);

figure(3);
markercolor=(filtdata(:,4)); %fildata(:,6)=quality of fit; filtdata(:,7)=no. of fit; filtdata(:,4)...signal
markercolor=[0 0 0];
plot(r_mol*1e6,filtdata(:,3),'.',r_coord*1e6,z_theory*1e9,'r');
hold on; 
plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz)],'r.');
plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz)],'r.');
scatter(r_mol*1e6,filtdata(:,3),1,markercolor); grid on; colorbar; 
xlabel('radial coord / µm');
ylabel('z /  nm');
title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(1)/1000,3) '/' num2str(xc(2)/1000) ' µm, std=' num2str(delta_z,2) 'nm, qual-thresh=' num2str(qual_thresh)]);
grid on;
ylim([0 250]);
colormap jet; 
hold off;