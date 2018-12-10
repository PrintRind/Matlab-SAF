%evaluate experimental molecule Data --> 3D localisations
clc
clear all; 
close all

global TS pxSize dist_max data_BG PSF_top_norm PSF_bottom_norm interp_incr x_data z_vec ux
%% user parameters 

pxSize=115;
dist_max=400;

%camera parameters
gain=100; 
amp=9.9; %electrons per count 
QE=0.95; %quantum efficiency
r_sphere=2.0e-3/2; %radius of sphere
ratio=0.93; %fixed energy-ratio between both images ratio=energy_topimage/energy_bottomimage

z_ini=80; %initial z-estimate (frame-no of PSF-stack)

%% loading PSF model top image

load('./PSFs/PSF5D_tot(top)_0-2-250nm_RI=1,45_dz=0_aberrfree.mat'); PSF_top=PSF5D;

if ndims(PSF_top)==3 %if "standard" 3D PSF is loaded (no interpolation in MLE necessary -> faster)
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
disp('done');

%% loading PSF model bottom image

load('./PSFs/PSF5D_UAF(bottom)_0-2-250nm_RI=1,45_dz=0_aberrfree.mat'); PSF_bottom=PSF5D;%DONALD

if ndims(PSF_bottom)==3 %if "standard" 3D PSF is loaded (no interpolation in MLE necessary -> faster)
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
disp('done');

%% combine models and calculate CRLB

fw=3; %frame-width (width of border around the PSF-images, needed to allow for x-y-localization (shifts))
       
%create coordinate system 
x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
[X,Y]=ndgrid(x,y);
clear x_data;
x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx, ny]=size(X); %size of molecule image

%defining background mask for initial BG-estimate: 
%BGmask=zeros(nx,ny);
%BGmask(1,:)=1; BGmask(end,:)=1; BGmask(:,1)=1; BGmask(:,end)=1;

uz=z_vec(2)-z_vec(1); %z-increment provided with PSF-models
PSF_tot=[squeeze(PSF_top(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2))) squeeze(PSF_bottom(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)))];

[CRBx,CRBy,CRBz]=fun_CRB(PSF_tot./repmat(trapz(trapz(PSF_tot,1),2),[size(PSF_tot,1) size(PSF_tot,2) 1]),ux,uz,10000,500);

figure(1)
plot(z_vec(1:end-1)*1e9,sqrt(CRBz));
title(['mean CRB z: ' num2str(mean(sqrt(CRBz)))]);
xlabel('z-pos/nm');
ylabel('CRB /nm');
grid on;

%% reading in experimental data (movie of blinking molecules)

clear data_top data_bottom
[Name,Path,~]=uigetfile('*.tif','choose TOP and BOTTOM tif stacks','C:\Users\q004aj\Desktop\2018-11-28','MultiSelect','on');
    image_no=length(imfinfo([Path Name{1}]));
    for m=1:image_no
      data_bottom(:,:,m)=flipud(double(imread([Path Name{1}],m)));
      data_top(:,:,m)=double(imread([Path Name{2}],m))/ratio; %weighting image energy with "ratio"
    end
disp('done');
    
% % EVER - background reduction
% % find minima in data stack
% V_min=min(data,[],3);
% V_BG1=(V_min+27.8)/0.9513 ; 
% imagesc(V_BG1); title('background');
% 
% data_corr=data-repmat(V_BG1,1,1,size(data,3));
% data_BG=V_BG1;

%% creating coordinate transform between TOP and BOTTOM stacks 
%manual selection of point-pairs 

img_frame=1; %frame to compare

I1=data_top(:,:,img_frame);
I2=data_bottom(:,:,img_frame);

%>>>>>>>>>>>>>> execute next line if you want to map <<<<<<<<<<<<<<<
%cpselect(I1/max(I1(:)),I2/max(I2(:))); %save points as "movingPoints", "fixedPoints";

%saving points to workspace
%save('transformpoints_simu.mat','movingPoints','fixedPoints');
%load transformpoints_simu.mat

%% load Thunderstorm locations found in the top image and calculating corresponding locations in the bottom image
%>>>>>import under varname "TS" as matrix-type
%-> check if molecule images are accurately cropped from the image stacks

img_no=1;

moving_pts_adj= cpcorr(movingPoints, fixedPoints, I1, I2); %fine-tuning
%using cross-correlation

tform = fitgeotrans(moving_pts_adj,fixedPoints,'affine'); %finding transformation from selected points
%[tmp,RB]=imwarp(I1,tform);

TS_t=TS; %initializing Thunderstorm coordinates for transformed (second) image
offset=[0 0];%[-6 2.5]*ux*1e9; %visually determined offset for crop-coordinates
[xT,yT] = transformPointsForward(tform,TS(:,3),TS(:,4));
%[TS_t(:,3),TS_t(:,4)]=worldToIntrinsic(RB,xT,yT);

TS_t(:,3:4)=[xT,yT]-repmat(offset,size(TS,1),1);

data_BG=zeros(size(data_top));
I_top=fun_thunderstorm_crop(TS,data_top,data_BG,(nx-1)/2,ux);
I_bottom=fun_thunderstorm_crop(TS_t,data_bottom,data_BG,(nx-1)/2,ux);

 
figure(1);
subplot(1,2,1);
imagesc(I_top(:,:,img_no)); axis equal; axis tight; 
title('image top');
subplot(1,2,2);
imagesc(I_bottom(:,:,img_no)); axis equal; axis tight; 
title('image bottom');


%% joint MLE on the combined images
 
%correctly normalizing PSFs: each z-slice of the JOINT PSF is normalized
energies=repmat((trapz(trapz(PSF_tot,1),2)),[size(PSF_tot,1),size(PSF_tot,2)/2,1,nxi,nyi]); %energies of joint PSFs for each z-slice
PSF_top_norm=PSF_top;%./energies; 
PSF_bottom_norm=PSF_bottom;%./energies;

px=TS(:,3)/(ux*1e9); %approx. position of molecules 
py=TS(:,4)/(ux*1e9);

qual_thresh=5; 

close all; 
idx_disc=[];
for mm=1:size(TS,1)
    
    est=fun_MLE_2channels(I_top(:,:,mm), I_bottom(:,:,mm), PSF_top_norm, PSF_bottom_norm, interp_incr,x_data,z_ini,qual_thresh,ratio);
    %est=[x1,y1,x2,y2,z,Sig,BG1,BG2,resnorm1,resnorm2]
      
    if isempty(est) %if resnorm is too large, function returns est=[]
        idx_disc=[idx_disc mm]; %indidces of images to be discarded
    else
        m=m+1;
        x1_est(m)=(py(mm)+est(1))*ux*1e9;
        y1_est(m)=(px(mm)+est(2))*ux*1e9; %px measures along the column-direction
        z_est(m)=interp1(1:length(z_vec),z_vec,est(5))*1e9; %estimated z-position in nm
        N_est(m)=est(6)/gain*amp/QE;
        BG1_est(m)=est(7)/gain*amp/QE;
        BG2_est(m)=est(8)/gain*amp/QE;
        resnorm_MLE(m)=est(9)+est(10);
        loc_no(m)=TS(mm,1); %number of localization (temporal order)
        frame_no(m)=TS(mm,2); %number of image frame
        disp(m);
    end
end
locdata_TS=[x1_est' y1_est' z_est' N_est' BG1_est' BG2_est' resnorm_MLE' loc_no' frame_no'];
locdata_TS(isnan(locdata_TS))=0;

%% saving localization results

Name='donald_2018-11-19';
save(['LOCDATA_' Name '.mat'],'locdata_TS','ux','NA','gain','amp','QE','z_vec','r_sphere','I_top','I_bottom');

%% -----data eval----- finding optimal sphere center and plotting z-estimates against distance from sphere center

finaldata=locdata_TS; 
I_sel=I_top; %selected images

%----------deleting erroneous entries / filtering data -----------
    r_cutoff=25e-6; %cutoff distance from sphere center in m
%manually defining sphere center
    xc=[1 1]*size(data_top,1)*ux/2*1e9; %center of image
    findcenter='n'; %automatic center-finding? set to 'y' or 'n'
%select temporal junks of data to see time trends
    no_c=0; %central local. no.
    no_delta=inf;  %half width, set inf to get all data
%selecting a defined azimuthal section of the sphere 
    phi_c=200; %central angle of cake-slice (in deg.)
    dphi=inf; %angular width of "cake-slice"
%delting entries with too high fit-errors (residuals)
    qual_thresh=inf;  %the lower the cutoff, the more strict the selection
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
    xc0=[1 1]*size(data_top,1)*ux*1e9/2;
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
xlabel('x / �m');
ylabel('y / �m');


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
scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,2,markercolor);
xlabel('x / �m');
ylabel('y / �m');
axis equal; axis tight; 
colormap jet; 
colorbar; 
title('loc. data, z-values color-coded');



% -----plotting radial distance versus estimated z-positions-----

if gain==1 || strcmp(Name{1}(1:4),'simu')  %if the raw-data stems from a simulation, gain>1 has no influence on noise!
    mean_signal=mean(filtdata(:,4));
    mean_bg=mean(filtdata(:,5));
else
    mean_signal=mean(filtdata(:,4))/2; %considering excess noise when gain is switched on; this effectively reduces the signal number by factor 2
    mean_bg=mean(filtdata(:,5))/2;
end
[CRBx,CRBy,CRBz]=fun_CRB(PSF_tot./repmat(trapz(trapz(PSF_tot,1),2),[size(PSF_tot,1) size(PSF_tot,2) 1]),ux,uz,mean_signal,mean_bg);

figure(3);
markercolor=(filtdata(:,4)); %fildata(:,6)=quality of fit; filtdata(:,7)=no. of fit; filtdata(:,4)...signal
%markercolor=zeros(length(filtdata),1);
plot(r_mol*1e6,filtdata(:,3),'.',r_coord*1e6,z_theory*1e9,'r');
hold on; 
plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz)],'r.');
plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz)],'r.');
scatter(r_mol*1e6,filtdata(:,3),3,markercolor); grid on; colorbar; 
xlabel('radial coord / �m');
ylabel('z /  nm');
title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(1)/1000,3) '/' num2str(xc(2)/1000) ' �m, std=' num2str(delta_z,2) 'nm, qual-thresh=' num2str(qual_thresh)]);
grid on;
ylim([0 250]);
colormap jet; 
hold off;