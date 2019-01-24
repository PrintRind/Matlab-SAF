%evaluate experimental molecule Data --> 3D localisations

clear all; 

%% user parameters 

%camera parameters
gain=50;  %set gain to 1 if gain is switched off
amp=9.9; %9.9; %electrons per count 
QE=0.95; %quantum efficiency
ux=115e-9; %effective camera pixel size
r_sphere=2.0e-3/2; %radius of sphere

qual_thresh=inf; %quality thrshold for MLE-fits (bad fits are discarded from the data); %the smaller the value the stricter is the selection.

%% loading PSF model

PSFpath='C:\Users\q004aj\Desktop\PSFs\';

%3D PSFs:
%load([PSFpath 'PSF_0-3-250nm_RI=1,45_defoc=-400nm_aberrfree.mat']); PSF=PSF_tot; %defocused PSF-model; should give better z-estimates

%5D PSFs: (faster)
%load([PSFpath 'PSF5D_0-2-250nm_RI=1,45_dz=-400_aberrfree_os3.mat']); PSF=PSF5D; %defocused PSF-model; should give better z-estimates
%load([PSFpath 'PSF5D_17x17_0-2-250nm_RI=1.35_dz=-500_aberrfree_os3.mat']); PSF=PSF5D; %defocused PSF-model; should give better z-estimates
load([PSFpath 'PSF5D_13x13_0-2-250nm_RI=1,45_dz=-400_aberrfree_os3.mat']); PSF=PSF5D; %defocused PSF-model; should give better z-estimates

ux=ux*os; 

%"defocus" Stack (with varying defocus but constant SAF/UAF relations, i.e.
%no change of molecule distance from the coverslip)
%load('PSF5D_defocus_dz=-700 to -200nm_RI=1,45_aberrfree.mat'); PSF=PSF5D; %defocused PSF-model; should give better z-estimates


if length(z_vec)>1
    z=z_vec; %PSF model describes an "SAF-PSF" (fixed objetive defocus, but varying molecule distances from the coverslip
else
    z=dz_vec; %PSF model describes a "defocus" PSF (fixed mol. distance, but varying objective defoci
end

z_ini=round(size(PSF,3)/2); %initial z-estimate (frame-no of PSF-stack)


if ndims(PSF)==3 %if "standard" 3D PSF is loaded (no interpolation in MLE necessary -> faster)
    [nx0,ny0,nz0]=size(PSF); %size of model    clear PSF_tot;
    %PSF: normalization of energy in each z-slice
    energies=(trapz(trapz(PSF,1),2)); 
    PSF_norm=PSF./repmat(energies,[size(PSF,1),size(PSF,2),1]);
    interp_incr=1; 
elseif ndims(PSF)==5 
    [nx0,ny0,nz0,nxi,nyi]=size(PSF); %size of model
    clear PSF5D;
    %PSF: normalization of energy in each z-slice
    energies=(trapz(trapz(PSF(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)),1),2)); 
    PSF_norm=PSF./repmat(energies,[size(PSF,1),size(PSF,2),1,nxi,nyi]);
end
disp('done');

%dz_PSF=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
%PSF_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF

fw=0; %frame-width (width of border around the PSF-images, needed to allow for x-y-localization (shifts)); can be set to zero if a 5D PSF is used (becaues no interpolation is required in the function fun_MLE)
       
%create coordinate system 
x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
[X,Y]=ndgrid(x,y);
clear x_data;
x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx, ny]=size(X); %size of molecule image

%defining background mask for initial BG-estimate: 
BGmask=zeros(nx,ny);
BGmask(1,:)=1; BGmask(end,:)=1; BGmask(:,1)=1; BGmask(:,end)=1;

%% reading in experimental data (movie of blinking molecules)
use_EVER='y'; %background subtraction using EVER ('y') or no background subtraction ('n')

clear data
[Name,Path,~]=uigetfile('*.tif','choose multipage tif data','C:\Users\q004aj\Documents\PROJEKTE\SAF-TIRF_Gerhard Sch�tz - Immunolog. Synapse\Data\2018-10-18\');
    image_no=length(imfinfo([Path Name]));
    progressbar;
    for m=1:image_no
       data(:,:,m)=double(imread([Path Name],m));
       progressbar(m/image_no);
    end
    
% EVER - background reduction
%find minima in data stack

if strcmp(use_EVER,'y') %use EVER bg-subtraction
    V_min=min(data,[],3);
    V_BG1=(V_min+27.8)/0.9513 ; 
    imagesc(V_BG1); title('background');

    %more complex version: 
    % V_f=mean(data,3); %mean value of image set
    % V_vf=V_f./min(data,[],3); %variation factor
    % V_vr=mean(V_vf,[],3);
    % A=0.9433*exp(0.03336*V_vr)-186100*exp(-15.89*V_vr);
    % B=-600.6*exp(-4.131*V_vr)-22.9*exp(-.2063*V_vr); 

    data_corr=data-repmat(V_BG1,1,1,size(data,3));
    data_BG=V_BG1;
    disp('done; EVER BG-corr. used');

else
    data_corr=data; 
    data_BG=zeros(size(data_corr));
    disp('done; no BG subtracted!');
end



%% optional: saving EVER-corrected stack, e.g. to be used with thunderstorm
% 
if strcmp(use_EVER,'y')
    for m=1:image_no
        imwrite(uint16(data_corr(:,:,m)),[Name(1:end-4) '_EVER.tif'],'WriteMode','Append');
    end
end
disp('done');

%%
imagesc(data_corr(:,:,5));

%% METHOD A : click & evaluate (only use for 3D-PSFs! Ohterwise, you must click with 1 pixel precision!)

img_no=4;

figure(1);
imagesc(data_corr(:,:,img_no)); colormap gray; axis equal; axis tight; 
title('exp. image');
[px, py] = getpts(gcf); %click on several molecules and press enter 
no_mols=length(px); 

clear tmp;
%extracting small molecule-images from camera frame
for m=1:no_mols
    data_tmp=data_corr(round((py(m)-floor((nx0-1)/2)+fw):(py(m)+ceil((nx0-1)/2))-fw),round((px(m)-floor((nx0-1)/2)+fw):(px(m)+ceil((nx0-1)/2)-fw)),img_no);
    bg_tmp=data_BG(round((py(m)-floor((nx0-1)/2)+fw):(py(m)+ceil((nx0-1)/2))-fw),round((px(m)-floor((nx0-1)/2)+fw):(px(m)+ceil((nx0-1)/2)-fw)));
    tmp(:,:,m)=data_tmp+mean(bg_tmp(:)); %adding mean background of the molecule frame (important for MLE)
end

data2=double(tmp/gain*amp/QE); %converting into units of photons; 
imagesc(data2(:,:,1)); title('fist exp. image / photons'); axis equal; axis tight; colorbar; 

%--- Estimating parameters ---
    
clear x_est y_est z_est z_est_nm BG_est N_est

for m=1:no_mols
        
    I=data2(:,:,m); %pick one experimental image

    est=fun_MLE(I,PSF_norm,interp_incr,x_data,z_ini,qual_thresh,1); %estimate molecule position 
    
    if isempty(est)==0
        x_est(m)=(round(py(m))+est(1))*ux*1e9;
        y_est(m)=(round(px(m))+est(2))*ux*1e9; %px measures along the column-direction
        z_est(m)=interp1(1:length(z),z,est(3))*1e9; %estimated z-position in nm
        N_est(m)=est(4)/gain*amp/QE;
        BG_est(m)=est(5)/gain*amp/QE;

        disp(m);
        imagesc(I); axis equal; axis tight; title(['sig./ bg =' num2str(N_est(m)) ' / ' num2str(BG_est(m))]); pause(0);
    end
    
end

figure(2); 
plot3(x_est/1e3,y_est/1e3,z_est,'o'); 
grid on;
hold on; 
zlabel('z / nm');
xlabel('x / �m');
ylabel('y / �m');

locdata_struct{img_no}=[x_est' y_est' z_est' N_est' BG_est'];
%save('loc_data_-750nm-model_-400nm-data_2018-09-26.mat','locdata_struct','NA','RI','z','ux');

%converting final localization data from struct-array into array
locdata=[];
for m=1:length(locdata_struct)
    locdata=[locdata; locdata_struct{m}];
end

disp(est);

%% METHOD B: using thunderstorm coordinates
% (works for 3D and 5D PSF models)

TS_multi = csvread('TS.csv',1,0);

if size(TS_multi,2)<10 %adding extra columns if TS-matrix has too few
    TS_multi=[TS_multi zeros(size(TS_multi,1),10-size(TS_multi,2))];
end
        
%data pre-filtering options
thresh_dist=00; %minimum distance of a molecule to the next in nm (others are deleted in the function "RmNN")
thresh_sigma=1*max(TS_multi(:,5)); %sigma of Gaussian fit in nm
thresh_uncert=1*max(TS_multi(:,10)); %
TS = RmNN(TS_multi(end,2),TS_multi,thresh_dist,thresh_sigma,thresh_uncert); %picking only those molecule images with no neighbours in the vicinity

%use thunderstorm to find molecules. copy and paste resulting localization
%table into Matlab under the varname "TS" (import as "MATRIX")
%TS=[id, frame, x, y, sigma, intensity, offset,chi2]

%saving or loading Thunderstorm (TS) localization data: 
%save('TS_2018-10-10_0nm.mat','TS');
%load TS_2018-09-26_-400nm.mat
TS_full=TS;
%TS(10000:end,:)=[];

I=fun_thunderstorm_crop(TS,data_corr,data_BG,(nx-1)/2,ux);
px=TS(:,3);%/(ux*1e9); %approx. position of molecules in nm
py=TS(:,4);%/(ux*1e9);

clear x_est y_est z_est N_est BG_est resnorm_MLE loc_no frame_no chi2;

m=0;
clear locdata_TS resnorm_MLE loc_no x_est y_est z_est N_est BG_est;

idx_disc=[];
%progressbar; %initialize progressbar
loc_no_tot=length(TS_full);

tic
for mm=1:loc_no_tot
    
    if mm<100
        showimage=1;
    else
        showimage=0;
    end
        
    est=fun_MLE(I(:,:,mm),PSF_norm,interp_incr,x_data,z_ini,qual_thresh,showimage);
    %est=fun_MLE_fast(I(:,:,mm),PSF_norm,interp_incr,z_ini,qual_thresh,[relpos(mm,:) TS(mm,6:7)/amp],showimage)
    
    if isempty(est) %if resnorm is too large, function returns est=[]
        idx_disc=[idx_disc mm]; %indidces of images to be discarded
    else
        m=m+1;
        %x_est(m)=py(mm);%
        %y_est(m)=px(mm);% 
        x_est(m)=(round(py(mm)/ux/1e9)+est(1))*ux*1e9; %note the round-operation! this is required because of the way the images are cropped in "fun_thunderstorm_crop"
        y_est(m)=(round(px(mm)/ux/1e9)+est(2))*ux*1e9; %px measures along the row-direction (horizontal direction in image)
        
        z_est(m)=interp1(1:length(z),z,est(3))*1e9; %estimated z-position in nm
        N_est(m)=est(4)/gain*amp/QE;
        BG_est(m)=est(5)/gain*amp/QE;
        resnorm_MLE(m)=est(6);
        loc_no(m)=TS(mm,1); %number of localization (temporal order)
        frame_no(m)=TS(mm,2); %number of image frame
        chi2(m)=TS(mm,8);
        uncert(m)=TS(mm,9);
        %progressbar(mm/loc_no_tot);
        
        %test of x-y-localization
%         figure(1); 
%         imagesc(I(:,:,mm)); colormap gray; 
%         hold on; 
%         plot(est(2)+round(nx/2),est(1)+round(nx/2),'rx');
%         hold off; 
%         pause;
    
    end
end
toc

I(:,:,idx_disc)=[];
disp('done');

locdata_TS=[x_est' y_est' z_est' N_est' BG_est' resnorm_MLE' loc_no' frame_no' chi2'];
z_drift=0;

%% saving localization results

save(['LOCDATA_fulldata_13x13PSF_dz-400nm_' Name '.mat'],'locdata_TS','ux','NA','gain','amp','QE','z','r_sphere','I');

%% -- data eval, method B (thunderstorm)

finaldata=locdata_TS; 

I_sel=I; %selected images

% %DRIFT CORRECTION - interpolating frame-dependent drift values from the Thunderstorm-curves
% drift_x=interp1(drift(2:end,1),drift(2:end,2),finaldata(:,8));
% drift_y=interp1(drift(2:end,3),drift(2:end,4),finaldata(:,8));
% plot(drift_x); hold on; plot(drift_y); xlabel('frame no.'); ylabel('drift/pix'); hold off; 
% finaldata(:,1)=finaldata(:,1)-1*drift_y*ux*1e9;
% finaldata(:,2)=finaldata(:,2)-1*drift_x*ux*1e9;
% 
% finaldata=finaldata0; 
% -----data eval----- finding optimal sphere center and plotting z-estimates against distance from sphere center

%finaldata=locdata;   %manually selected molecules

%----------deleting erroneous entries / filtering data -----------
    r_cutoff=18e-6; %cutoff distance from sphere center in m; choose a smaller cutoff size(e.g. 18�m) if center of sphere has to be determinedd
%manually defining sphere center
    xc=[6 7.17]*1e3; %[x,y]
    %xc=[1 1]*size(data,1)/2*ux*1e9;
    findcenter='y'; %automatic center-finding? set to 'y' or 'n'
%select temporal junks of data to see time trends
    no_c=200000; %central local. no.
    no_delta=50000;  %half width, set inf to get all data
%selecting a defined azimuthal section of the sphere 
    phi_c=0; %central angle of cake-slice (in deg.)
    dphi=inf; %angular width of "cake-slice"
%delting entries with too high fit-errors (residuals)
    qual_thresh=10;  %the lower the cutoff, the more strict the selection
%filtering data: photon threshold
    photon_loBND=0;
    photon_hiBND=inf;
%----------------------------------------

%delete molecules with zero signal photons
idx=find(finaldata(:,4)==0);
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

% %delete entries that run into upper z-boundaries
% idx=finaldata(:,3)>=max(z)*1e9;
% finaldata(idx,:)=[]; 
% I_sel(:,:,idx)=[]; 
% 
% %delete entries that run into lower z-boundaries
% idx=finaldata(:,3)<=min(z)*1e9+1;
% finaldata(idx,:)=[]; 
% I_sel(:,:,idx)=[]; 

%selecting temporal sections
idx=abs(finaldata(:,7)-no_c)>no_delta;
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
idx=finaldata(:,6)>qual_thresh*min(finaldata(:,6));
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%filtering data: photon threshold
idx=(photon_hiBND<finaldata(:,4));
finaldata(idx,:)=[];
idx=(finaldata(:,4)<photon_loBND);
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

filtdata=finaldata;


%---calculating cramer-rao lower bounds---

if ndims(PSF)==5 
    PSF3D=PSF(:,:,:,ceil((size(PSF,4)+1)/2),ceil((size(PSF,4)+1)/2));
else
    PSF3D=PSF; 
end
    
if gain==1 || strcmp(Name(1:4),'simu')  %if the raw-data stems from a simulation, gain>1 has no influence on noise!
    mean_signal=mean(filtdata(:,4));
    mean_bg=mean(filtdata(:,5));
else
    mean_signal=mean(filtdata(:,4))/2; %considering excess noise when gain is switched on; this effectively reduces the signal number by factor 2
    mean_bg=mean(filtdata(:,5))/2;
end
[CRBx,CRBy,CRBz]=fun_CRB(PSF3D./repmat(sum(sum(PSF3D,1),2),[size(PSF3D,1) size(PSF3D,2) 1]),ux,z(2)-z(1),mean_signal,mean_bg);
disp('------------------');
disp('Cramer-Rao bounds:');
disp(['min. sqrt(CRLB-x): ' num2str(min(sqrt(CRBx(:))),3)]);
disp(['min. sqrt(CRLB-y): ' num2str(min(sqrt(CRBy(:))),3)]);
disp(['min. sqrt(CRLB-z): ' num2str(min(sqrt(CRBz(:))),3)]);

%%---automatically finding optimal sphere-center: 
r_coord=linspace(0,30,size(PSF,3))*1e-6; %radial coord. in m
z_theory=r_sphere-sqrt(r_sphere^2-(r_coord).^2);

clear r_mol
c=figure(3);
close(c);

z_error=@(x) std(interp1(r_coord,z_theory,sqrt((filtdata(:,1)-x(1)).^2+(filtdata(:,2)-x(2)).^2)*1e-9,'linear','extrap')*1e9-filtdata(:,3));

if strcmp(findcenter,'y')
    batchsize=1000; %batchsize for determining sphere-center
    N=ceil(length(filtdata)/batchsize);
    xc0=xc;

    for m=1:N %determining sphere center for each batch of data
       batch=filtdata((1+(m-1)*batchsize):min((1+m*batchsize),length(filtdata)),:);
       xc(m,:)=fun_find_sphere_center(xc0,r_coord,z_theory,batch);
       xc0=xc(m,:);
       r_mol{m}=sqrt((batch(:,1)-(xc(1))).^2+(batch(:,2)-(xc(2))).^2)*1e-9; %radial molecule coords. in meter
       delta_z(m)=z_error(xc(m,:)); %standard deviation in nm 
       z_batchmean(m)=mean(batch(:,3)); %mean z-value of batch; to identify z-drifts
       
       %----MAIN PLOT: plotting radial distance versus estimated z-positions-----
        figure(3);
        markercolor=(batch(:,8)); %fildata(:,6)=quality of fit; filtdata(:,7)=no. of fit; filtdata(:,4)...signal
        %markercolor=[0 0 0];
        plot(r_coord*1e6,z_theory*1e9,'r');
        hold on; 
        plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz)],'r.');
        plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz)],'r.');
        scatter(r_mol{m}*1e6,batch(:,3)-(z_batchmean(m)-z_batchmean(1)),3,markercolor); grid on; colorbar; 
        xlabel('radial coord / �m');
        ylabel('z /  nm');
        title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(m,1)/1000,3) '/' num2str(xc(m,2)/1000) ' �m, std=' num2str(delta_z(m),2) 'nm, qual=' num2str(qual_thresh) ', sig/bg=' num2str(mean(filtdata(:,4)),4) '/' num2str(mean(filtdata(:,5)),3)]);
        grid on;
        ylim([0 250]);
        colormap jet; 
        
    end
        figure(31);
        subplot(2,1,1);
        plot(xc,'.-');
        xlabel('batch no.');
        ylabel('sphere center');
        title('temporal development of sphere center');
        subplot(2,1,2);
        plot(z_batchmean); title('z-drift');
        xlabel('batch no.');
        ylabel('nm');
        
        %z_drift_fit=
        
        
else
        r_mol=sqrt((filtdata(:,1)-(xc(1))).^2+(filtdata(:,2)-(xc(2))).^2)*1e-9; %radial molecule coords. in meter
        delta_z=z_error(xc);
        
        figure(3);
        markercolor=(filtdata(:,8)); %fildata(:,6)=quality of fit; filtdata(:,7)=no. of fit; filtdata(:,4)...signal
        %markercolor=[0 0 0];
        plot(r_mol*1e6,filtdata(:,3),'.',r_coord*1e6,z_theory*1e9,'r');
        hold on; 
        plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz)],'r.');
        plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz)],'r.');
        scatter(r_mol*1e6,filtdata(:,3),3,markercolor); grid on; colorbar; 
        xlabel('radial coord / �m');
        ylabel('z /  nm');
        title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(1)/1000,3) '/' num2str(xc(2)/1000) ' �m, std=' num2str(delta_z,2) 'nm, qual=' num2str(qual_thresh) ', sig/bg=' num2str(mean(filtdata(:,4)),4) '/' num2str(mean(filtdata(:,5)),3)]);
        grid on;
        ylim([0 250]);
        colormap jet; 
        hold off;
end
figure(3);
hold off;

%z_target=interp1(r_coord,z_theory,r_mol)*1e9; %traget z-value in nm

%----3D scatterplot: show all filtered data----
figure(1); 
%plot3(locdata_TS(:,1)/1e3,locdata_TS(:,2)/1e3,locdata_TS(:,3),'o'); grid on;
markersize=1; %(locdata_TS(:,4)+1)/1e3;
markercolor=filtdata(:,3); 
scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,filtdata(:,3),markersize,markercolor); grid on;
colormap jet; 
colorbar; 
title('local. data, thunderstorm');
zlabel('z / nm');
xlabel('x / �m');
ylabel('y / �m');

% figure(2);
% z_dev=filtdata(:,3)-(r_sphere-sqrt(r_sphere.^2-r_mol.^2))*1e9; %deviation of est. z-pos. from the sphere surface
% markercolor=z_dev;
% scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,z_dev,markersize,markercolor); grid on;
% colormap jet; 
% colorbar; 
% title('deviation from sphere surface');
% zlabel('z / nm');
% xlabel('x / �m');
% ylabel('y / �m');

figure(4);
subplot(2,1,1); 
hist(filtdata(:,4),20);
grid on;
title(['signal, mean=' num2str(mean(filtdata(:,4)))]);
xlabel('photons');
ylabel('number of localizations');
subplot(2,1,2);
hist(filtdata(:,5),20);
grid on;
title(['background, mean=' num2str(mean(filtdata(:,5)),2)]);
xlabel('photons');
ylabel('number of localizations');

% ----2D scatterplot with color-coded z-value----
figure(6); 
markersize=1;
markercolor=filtdata(:,3);% markercolor=uint8(markercolor);
scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,markersize,markercolor);
xlabel('x / �m');
ylabel('y / �m');
axis equal; axis tight; 
colormap jet; 
colorbar; 
title('loc. data, z-values color-coded');
% %
% 
% %---2D-scatterplot; signal-photon-no. color-coded---
% figure(7);
% markercolor=filtdata(:,4); 
% %markercolor=uint8(markercolor/max(markercolor(:))*255*2);
% scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,2,markercolor);
% xlabel('x / �m');
% ylabel('y / �m');
% axis equal; axis tight; 
% colormap jet; 
% colorbar; 
% title('loc. data, photon-counts color-coded');

%% plotting only "heroes", i.e. molecules that are "on" for serveral consecutive frames


clear hero
minlength=5; 
bridgelength=5; 
hero=fun_findheroes(filtdata,minlength,bridgelength); %returns filtdata of heroes in a structure array 
no_heroes=length(hero);

for m=1:no_heroes;
        hero(m).mean=mean(hero(m).locdata);
        hero(m).std=std(hero(m).locdata);
end

hero_means=vertcat(hero.mean);
hero_sigmas=vertcat(hero.std);

%plotting 3D positions of heroes
figure(11); 
scatter3(hero_means(:,1),hero_means(:,2),hero_means(:,3),10);
title(['mean 3D positions of ' num2str(no_heroes) ' heros; minlength=' num2str(minlength) '; bridgelength=' num2str(bridgelength)]);

%finding best sphere center and plotting z-positions versus radial
%coordinate
xc_heroes=fun_find_sphere_center(xc0,r_coord,z_theory,hero_means); %finding sphere center based on hero mean positions



%calculating CRLBs for heroes
if gain==1 || strcmp(Name(1:4),'simu')  %if the raw-data stems from a simulation, gain>1 has no influence on noise!
    mean_signal_heroes=mean(hero_means(:,4));
    mean_bg_heroes=mean(hero_means(:,5));
else
    mean_signal_heroes=mean(hero_means(:,4))/2;
    mean_bg_heroes=mean(hero_means(:,5))/2;
end
[CRBx_h,CRBy_h,CRBz_h]=fun_CRB(PSF3D./repmat(sum(sum(PSF3D,1),2),[size(PSF3D,1) size(PSF3D,2) 1]),ux,z(2)-z(1),mean_signal_heroes,mean_bg_heroes);


figure(12); 
subplot(2,1,1);
r_heroes=sqrt((hero_means(:,1)-xc_heroes(1)).^2+(hero_means(:,2)-xc_heroes(2)).^2);
errorbar(r_heroes*1e-3,hero_means(:,3),hero_sigmas(:,3),'.');
hold on; 
plot(r_coord*1e6,z_theory*1e9,'r')
plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz_h)],'r.');
plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz_h)],'r.');
hold off; 
ylabel('z /nm');
xlabel('r /�m');
xlim([0 20]); 
grid on; 
title(['z-pos. of ' num2str(no_heroes) ' heros; minlength=' num2str(minlength) '; bridgelength=' num2str(bridgelength)]);
subplot(2,1,2);
plot(r_heroes*1e-3,hero_sigmas(:,3),'.');
hold on; 
plot(r_coord*1e6,[0 sqrt(CRBz_h)],'r.-');
hold off; 
xlabel('radial distance /�m');
ylabel('\sigma_z');
xlim([0 20]); 
  
%% show details to the CURSOR-selected molecule image

h=figure(3);
dcm_obj = datacursormode(h);
c_info = getCursorInfo(dcm_obj); 

XData=c_info.Target.XData;
YData=c_info.Target.XData;

selected_loc_no=find([XData' YData']==[c_info.Position(1:2)]); %c_info.Position(1))

figure(33)
imagesc(I_sel(:,:,selected_loc_no)); axis equal; axis tight; 
title(['z=' num2str(filtdata(selected_loc_no,3),4)]);

%% -----polynomial fit through data to provide better precision estimate-----
% %fitcoefs=polyfit(r_mol*1e6,filtdata(:,3),2);
% %fitcurve=polyval(fitcoefs,r_mol*1e6);
ft=fittype('a*x^6+b*x^4+c*x^2+d');
[fitobj,~,~,~]=fit(r_mol*1e6,filtdata(:,3),ft); %[fitobj, goodness, output, convmsg]=fit(x,y,ft)
coeff1=fitobj.a;
coeff2=fitobj.b;  
coeff3=fitobj.c;
coeff4=fitobj.d;
fitcurve=polyval([coeff1 0 coeff2 0 coeff3 0 coeff4],r_mol*1e6);                            %Calc residuals
figure(8);

plot(r_mol*1e6,filtdata(:,3),'.',r_mol*1e6,fitcurve,'.');
delta_z2=std(filtdata(:,3)-fitcurve)
title(['polynomial fit; std=' num2str(delta_z2,2) ' nm; qual-thresh=' num2str(qual_thresh)]);


%% --- click into 2D scatterplot and calculate precisions of closely neighboured localizations

figure(6);
[px, py] = getpts(gcf);
boxrad=50; %in nm; side length of squared selection "box" around selected point

for m=1:length(px)
    idx=find(logical(abs(filtdata(:,1)-px(m)*1e3)<=boxrad) & logical(abs(filtdata(:,2)-py(m)*1e3)<=boxrad)); %finding indices of localizations that are within the selected area
    
    figure(9);
    markersize=2;
    markercolor=filtdata(idx,3);% markercolor=uint8(markercolor);

    scatter3(filtdata(idx,1)/1e3,filtdata(idx,2)/1e3,filtdata(idx,3),markersize,markercolor);
    xlabel('x / �m');
    ylabel('y / �m');
    zlabel('z / nm');
    x_mean(m)=mean(filtdata(idx,1)); 
    x_std(m)=std(filtdata(idx,1));
    y_mean(m)=mean(filtdata(idx,2)); 
    y_std(m)=std(filtdata(idx,2));
    z_mean(m)=mean(filtdata(idx,3)); 
    z_std(m)=std(filtdata(idx,3));
    disp('----------------');
    disp(['x-pos (mean/std) = ' num2str(x_mean(m)/1e3,3) '�m /' num2str(x_std(m),3) 'nm']);
    disp(['y-pos (mean/std) = ' num2str(y_mean(m)/1e3,3) '�m /' num2str(y_std(m),3) 'nm']);
    disp(['z-pos (mean/std) = ' num2str(z_mean(m),3) 'nm /' num2str(z_std(m),3) 'nm']);

    figure(10);
    subplot(3,1,1);
    plot(filtdata(idx,8),filtdata(idx,1),'.');
    title('temporal trend - x');
    xlabel('frame no.'); grid on; 
    ylabel('nm');
    subplot(3,1,2);
    plot(filtdata(idx,8),filtdata(idx,2),'.');
    title('temporal trend - y');
    xlabel('frame no.');grid on; 
    ylabel('nm');
    subplot(3,1,3); 
    plot(filtdata(idx,8),filtdata(idx,3),'.');
    title('temporal trend - z');
    xlabel('loc. no.');grid on; 
    
    pause(0.1);
end
%save([Name '_eval.mat'],'x_mean','y_mean','z_mean','x_std','y_std','z_std','filtdata');

%% optional: spline fit to extract z-drift

drift_spline = fit(filtdata(idx,8),filtdata(idx,3),'smoothingspline','SmoothingParam',35e-2);
figure(100);
plot(drift_spline,filtdata(idx,8),filtdata(idx,3),'.'); xlabel('frame no.'); ylabel('z /nm');
tmp = ppval(drift_spline.p,filtdata(idx,8)); %z-correction curve
z_drift=interp1(filtdata(idx,8),tmp,filtdata(:,8),'linear','extrap');
z_drift=z_drift-mean(z_drift(:));


%% ---- z-drift correction using the blue-laser feedback signal recorded in Labview
%load drift-file (txt-file exported by Labview program "piezo control")
%drag & drop into workspace, import as "table" datatype
%call variable "d"

t0=datetime('16:01:20');  %recording start time of tif-stack (add 2(winter) or 3 hours! -> Andor software is 3hrs behind!) 
dt=0.05969; % 'Kinetic cycle time', (in seconds) take from tif-stack info                             
driftmag=+20; %drift factor in nm per drift signal unit; 


rectime=filtdata(:,8)*dt; %relative recording time of every location in seconds
t_end=t0+seconds(round(rectime(end))); 
drifttime=table2array(d(:,1));
driftsignal=table2array(d(:,2));

idx0=find(drifttime==t0, 1 ); %finding time point in dirft signal that matches the tif recordign start
%dt_drift=1/length(idx0); %time interval of a drift-value recordings
dt_drift=(seconds(drifttime(end)-drifttime(1)))/(find(drifttime==drifttime(end), 1 )-idx0); %time interval of a drift-value recordings

idx0(2:end)=[]; % first entry is the correct one
idx_end=round(rectime(end)/dt_drift);

z_drift=interp1(driftsignal,filtdata(:,8))*driftmag; %z-drift for a given camera frame

figure(11); 
plot(filtdata(:,8),z_drift,'.');
grid on; 
xlabel('frame no.');
ylabel('drift / nm');

%% calculate z-precision after drift correction:
z_dev_driftcorr=filtdata(:,3)-z_drift-(r_sphere-sqrt(r_sphere.^2-r_mol.^2))*1e9; %deviation of est. z-pos. from the sphere surface
delta_z3=std(z_dev_driftcorr);


%-----plotting drift-corrected z-data----- compare with fig3 
figure(33);
markercolor=(filtdata(:,7)); %fildata(:,6)=quality of fit; filtdata(:,7)=no. of fit
%markercolor=zeros(length(filtdata),1);
plot(r_mol*1e6,filtdata(:,3),'.',r_coord*1e6,z_theory*1e9,'r');
hold on; 
plot(r_coord*1e6,z_theory*1e9+[0 sqrt(CRBz)],'r.');
plot(r_coord*1e6,z_theory*1e9-[0 sqrt(CRBz)],'r.');
scatter(r_mol*1e6,filtdata(:,3)-z_drift,3,markercolor); grid on; colorbar; 
xlabel('radial coord / �m');
ylabel('z /  nm');
title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(1)/1000,3) '/' num2str(xc(2)/1000) ' �m, std=' num2str(delta_z3,2) 'nm, qual-th.=' num2str(qual_thresh) ', drift corr']);
grid on;
ylim([0 250]);
colormap jet; 
hold off;

%----3D scatterplot: show all filtered data----
figure(44); 
%plot3(locdata_TS(:,1)/1e3,locdata_TS(:,2)/1e3,locdata_TS(:,3),'o'); grid on;
markersize=3; %(locdata_TS(:,4)+1)/1e3;
markercolor=filtdata(:,7); 
scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,filtdata(:,3)-z_drift,markersize,markercolor); grid on;
colormap jet; 
colorbar; 
title('local. data, thunderstorm - drift corrected');
zlabel('z / nm');
xlabel('x / �m');
ylabel('y / �m');
zlim([0 250]);

