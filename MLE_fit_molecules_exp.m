%evaluate experimental molecule Data --> 3D localisations

clear all;  
close all; 

%% user parameters 

mode='single'; %set mode: 'biplane' or 'single'
                %NOTES for Biplane imaging: 
                %-) images of the 2nd channel will be flipped! (to ease
                %coord-transform)
preflip='LR';  %for biplane imaging: flips images of 2nd channel left-right LR or UD; set to '' if no flip is desired

%choose PSF file(s)
[PSF1_name,PSF1_path,~]=uigetfile('*.mat','choose PSF file for 1st channel.');
load([PSF1_path PSF1_name]); 
PSF1=PSF;

if strcmp(mode,'biplane')
    [PSF2_name,PSF2_path,~]=uigetfile('*.mat','choose PSF file of 2nd channel.');
    load([PSF2_path PSF2_name]); 
    PSF2=PSF;
end

%camera parameters
cam=camera('orca fusion')
obj=objective('olympus 1.49')
f_tube=200e-3; %focal length of the tube lens that was used in the actual experiment (200e-3 for Innsbruck setup)

ux=cam.pixsize/obj.M*obj.f_tube/f_tube; %effective camera pixel size
r_sphere=15e-6; %radius of sphere (only relevant for sphere-sample)
sample='nanoruler'; %type of sample, e.g. 'sphere', 'cos7', etc.
qual_thresh=inf; %quality thrshold for MLE-fits (bad fits are discarded from the data); %the smaller the value the stricter is the selection.

disp('PSF loaded');
if PSF.os*PSF.ux~=ux
    disp('WARNING! discrepancy of pixel sizes!')
end


%% reading in experimental data (movie of blinking molecules)

use_EVER='n'; %background subtraction using EVER ('y') or no background subtraction ('n')

clear data1 data2
[Name1,Path1,~]=uigetfile('*.tif','choose multipage tif data of channel 1:','C:\Users\q004aj\Documents\PROJEKTE\SAF-TIRF_Gerhard Schütz - Immunolog. Synapse\Data\2018-10-18\');
ULC1=[10 20]; %row and column indices of the upper left corner ot the image stack with respect to the entire camera chip (important for sCMOS)

if strcmp(mode,'biplane')
    [Name2,Path2,~]=uigetfile('*.tif','choose multipage tif data of channel 2:','C:\Users\q004aj\Documents\PROJEKTE\SAF-TIRF_Gerhard Schütz - Immunolog. Synapse\Data\2018-10-18\');
    ULC2=[10 20]; %row and column indices of the upper left corner ot the image stack with respect to the entire camera chip (important for sCMOS)
end

image_no=length(imfinfo([Path1 Name1]));

%------------load images of 1st channel-------------

progressbar;
for m=1:image_no
  data1(:,:,m)=double(imread([Path1 Name1],m));
  progressbar(m/image_no);
end

%data=bfopen([Path Name]);

% EVER - background reduction
%find minima in data stack

if strcmp(use_EVER,'y') %use EVER bg-subtraction
    [data1_corr, data_BG]=fun_EVER(data); 
else
    data1_corr=data1; 
    data1_BG=zeros(size(data1_corr));
end

clear data1;

%------------load images of 2nd channel-------------

if strcmp(mode,'biplane')
    
    progressbar;
    for m=1:image_no
      if strcmp(preflip,'LR')
        data2(:,:,m)=fliplr(double(imread([Path2 Name2],m)));
      elseif strcmp(preflip,'UD')
        data2(:,:,m)=flipud(double(imread([Path2 Name2],m)));
      else
        data2(:,:,m)=(double(imread([Path2 Name2],m))); 
      end
      progressbar(m/image_no);
    end
    
    if strcmp(use_EVER,'y') %use EVER bg-subtraction
        [data2_corr, data2_BG]=fun_EVER(data2); 
    else
        data2_corr=data2; 
        data2_BG=zeros(size(data2_corr));
    end
    
    clear data2;
end


disp('all images loaded');

%% optional: saving EVER-corrected stack, e.g. to be re-used with thunderstorm
% 
if strcmp(use_EVER,'y')
    for m=1:image_no
        imwrite(uint16(data1_corr(:,:,m)),[Name1(1:end-4) '_EVER.tif'],'WriteMode','Append');
    end

    if strcmp(mode,'biplane')
        for m=1:image_no
            imwrite(uint16(data2_corr(:,:,m)),[Name2(1:end-4) '_EVER_CH2.tif'],'WriteMode','Append');
        end
    end
end

disp('done');


%% for biplane imaging: creating coordinate transform between the two stacks 
%manual selection of point-pairs 

if strcmp(mode,'biplane')

    if isfile('transformpoints.mat') %if a coord. transform file is already contained in the current folder it is loaded
        load transformpoints.mat
        disp('transform matrix found and loaded.');
    else  %otherwise, create a new transform matrix by clicking on matching points in the top and bottom images

        prompt='no transform matrix found. Create one by clicking on matched point pairs. Do you have specific calib. images? (e.g. containing beads) (y/n)';
        answer = inputdlg(prompt);

        if strcmp(answer{1},'n')
            img_frame=1; %frame to compare points
            calib_image1=data1_corr(:,:,img_frame);
            calib_image2=data2_corr(:,:,img_frame);
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

end


%% ESTIMATE using thunderstorm coordinates

[Name_TS,Path_TS,~]=uigetfile('*.csv','choose Thunderstorm csv file:');
tmp=csvread([Path_TS Name_TS],1,0);

%TS=tmp(1:5000,:); %select parts of the entire data
TS=tmp;

if size(TS,2)<10 %adding extra columns if TS-matrix has too few
    TS=[TS zeros(size(TS,1),10-size(TS,2))];
end
        
%----data pre-filtering------ 
thresh_dist=1500; %minimum distance of a molecule to the next in nm (others are deleted in the function "RmNN")
thresh_sigma=1*max(TS(:,5)); %sigma of Gaussian fit in nm
thresh_uncert=1*max(TS(:,10)); %
disp('prefiltering data...');
TS_filt = RmNN(TS(end,2),TS,thresh_dist,thresh_sigma,thresh_uncert); %picking only those molecule images with no neighbours in the vicinity

px=TS_filt(:,3); %approx. position of molecules in nm
py=TS_filt(:,4);
loc_no_tot=length(TS_filt);
z_ini=0.5*PSF1.Nz*PSF1.uz; %initial z-estimate

fw=1; %frame width; size difference of cropped molecule images to PSF size; mol. images are smaller by [2*fw 2*fw]; should be set to 1 (if zero, "discretization artifacts" can occur in the estimated x,y,z coords)
nx=PSF1.Nx/PSF1.os-2*fw;
ny=PSF1.Ny/PSF1.os-2*fw;


if strcmp(mode,'biplane')

    %---------------calculating coords of locations in the bottom image---------
    img_no=1; %choose test image from stack
    tform = fitgeotrans(movingPoints,fixedPoints,'affine'); %finding transformation from selected points
    TS_t=TS_filt; %initializing Thunderstorm coordinates for transformed (second) image
    [xT,yT] = transformPointsForward(tform,TS_filt(:,3)/(ux*1e9),TS_filt(:,4)/(ux*1e9));

    %<<<<<<<fine-adjust centering of bottom image by varying offset-values
    offset=[0 0]*ux*1e9; %visually determined offset for crop-coordinates
    TS_t(:,3:4)=[xT,yT]*(ux*1e9)-repmat(offset,size(TS_filt,1),1);

    idx=TS_filt(:,2)==img_no;
    
    figure(2); 
    subplot(2,1,1);
    imagesc(data1_corr(:,:,img_no)); hold on; title(['CH-A image no.' num2str(img_no) ' with overlaid TS-coords.']);
    plot(TS_filt(idx,3)/(ux*1e9),TS_filt(idx,4)/(ux*1e9),'r.'); 
    hold off;
    subplot(2,1,2);
    imagesc(data2_corr(:,:,img_no)); hold on; title(['CH-B image no.' num2str(img_no) ' with transformed TS-coords.']);
    plot(TS_t(idx,3)/(ux*1e9),TS_t(idx,4)/(ux*1e9),'r.'); 
    hold off;

    data_BG=zeros(size(data1_corr));
    [I1, ULC1_stack]=fun_thunderstorm_crop(TS_filt,data1_corr,data_BG,(nx-1)/2,ux); %crops molecule images in the top-image
    [I2, ULC2_stack]=fun_thunderstorm_crop(TS_t,data2_corr,data_BG,(nx-1)/2,ux); %crops molecule images in the bottom-image

    %delting locs where I_top or I_bottom is zero and deleting them (can happen if a
    %rotation is involved and the loc is close to the image border)
    idx=find((squeeze(sum(sum(I1,1),2))==0)|(squeeze(sum(sum(I2,1),2))==0)); 
    TS_filt(idx,:)=[];
    TS_t(idx,:)=[];
    I1(:,:,idx)=[];
    I2(:,:,idx)=[];

    I{1}=I1;
    I{2}=I2; 
   
    tic
    est=fun_MLE3([PSF1 PSF2],I,z_ini);
    toc
    
    %adding coarse positional information of channel 2 from Thunderstorm
    px=TS_filt(:,3); %approx. position of molecules in nm
    py=TS_filt(:,4);
    x2_est=(round(py/ux/1e9)+est(:,7)/ux)*ux*1e9; %note the round-operation! this is required because of the way the images are cropped in "fun_thunderstorm_crop"
    y2_est=(round(px/ux/1e9)+est(:,8)/ux)*ux*1e9; %px measures along the row-direction (horizontal direction in image)

     
else %single plane mode

    %----cropping single molecule images----
    [I,ULC_stack]=fun_thunderstorm_crop(TS_filt,data1_corr,data1_BG,(nx-1)/2,ux);

    %data.images=I;
    %data.
    
    tic
    est=fun_MLE3(PSF,I,z_ini,cam);
    toc

end

%adding coarse positional information from Thunderstorm
px=TS_filt(:,3); %approx. position of molecules in nm
py=TS_filt(:,4);
loc_no_tot=length(TS_filt);
x_est=(round(py/ux/1e9)+est(:,1)/ux)*ux*1e9; %note the round-operation! this is required because of the way the images are cropped in "fun_thunderstorm_crop"
y_est=(round(px/ux/1e9)+est(:,2)/ux)*ux*1e9; %px measures along the row-direction (horizontal direction in image)

z_est=est(:,3)*1e9; %
N_est=est(:,4);%/cam.gain*cam.amp/cam.QE;
BG_est=est(:,5);%/cam.gain*cam.amp/cam.QE;
res_est=est(:,6); %resnorm

%adding data from the Thunderstorm file
loc_no=TS_filt(:,1); %number of localization (temporal order)
frame_no=TS_filt(:,2); %number of image frame
chi2=TS_filt(:,8);
uncert=TS_filt(:,9);

if strcmp(mode,'biplane')
    locdata=[x_est y_est z_est N_est BG_est loc_no frame_no chi2 uncert res_est x2_est y2_est];
else
    locdata=[x_est y_est z_est N_est BG_est loc_no frame_no chi2 uncert res_est];
end

z_drift=0;

disp('done');

%% saving localization results

if exist('Z_aberr2')~=1
    Z_aberr2=0;
end

if exist('Z_amp')~=1
    Z_amp=0;
end

save(['LOCDATA_-300nm_z-ini=' num2str(z_ini) '_mindist=' num2str(thresh_dist) '_' mode '.mat'],'locdata','PSF','lens','cam','I','Z_amp');
%save('filtdata_RI1,38_NA1,5_qt15.mat','filtdata');

%% -- data eval, method B (thunderstorm)

%NOTE: if you notice a kind of discretization artefact of x and y positions at
%higher values of x and y (e.g. to the upper right corner), this is
%probably due to a pixelsize mismatch between Thunderstorm and this MATLAB
%program (e.g. ux=115nm vs. ux=117nm) 

%NOTE: for execution of this section, a PSF model must be loaded too

if exist('sample')~=1 
    sample='';
end

finaldata=locdata; 
I_sel=I; %selected images

% %DRIFT CORRECTION - interpolating frame-dependent drift values from the Thunderstorm-curves
%drag "drift.csv" into workspace, slect "import as matrix" and select only
%the last 4 columns.
if exist('drift')
    drift_x=interp1(drift(2:end,1),drift(2:end,2),finaldata(:,7));
    drift_y=interp1(drift(2:end,3),drift(2:end,4),finaldata(:,7));
    plot(drift_x); hold on; plot(drift_y); xlabel('frame no.'); ylabel('drift/pix'); hold off; 
    finaldata(:,1)=finaldata(:,1)-1*drift_y*ux*1e9;
    finaldata(:,2)=finaldata(:,2)-1*drift_x*ux*1e9;
    disp(' ');
    disp('drift corrected');
end

% finaldata=finaldata0; 
% -----data eval----- finding optimal sphere center and plotting z-estimates against distance from sphere center

%finaldata=locdata;   %manually selected molecules

%----------deleting erroneous entries / filtering data -----------
    r_cutoff=3e-6; %cutoff distance from sphere center in m; choose a smaller cutoff size(e.g. 18µm) if center of sphere has to be determinedd
%manually defining sphere center
    xc=[4 4]*1e3; %[x,y]
    %xc=[1 1]*200/2*ux*1e9; %center of image (for simus)
    findcenter='n'; %automatic center-finding? set to 'y' or 'n'
%select temporal junks of data to see time trends
    no_c=2000; %central local. no.
    no_delta=inf;  %half width, set inf to get all data
%selecting a defined azimuthal section of the sphere 
    phi_c=0; %central angle of cake-slice (in deg.)
    dphi=inf; %angular width of "cake-slice"
%delting entries with too high fit-errors (residuals)
    qual_thresh=20;  %the lower the cutoff, the more strict the selection
%filtering data: photon threshold
    photon_loBND=500;
    photon_hiBND=inf;
%----------------------------------------

%delete molecules with zero signal photons
idx=find(finaldata(:,4)==0);
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%delete entries that run into upper z-boundaries
% idx=finaldata(:,3)>=0.99*PSF.uz*(PSF.Nz-1)*1e9;
% finaldata(idx,:)=[]; 
% I_sel(:,:,idx)=[]; 

%delete entries that run into lower z-boundaries
% idx=finaldata(:,3)<=min(z)*1e9+1;
% finaldata(idx,:)=[]; 
% I_sel(:,:,idx)=[]; 

%selecting temporal sections
idx=abs(finaldata(:,7)-no_c)>no_delta;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

% %delete entries which are too far away from the sphere center
% if strcmp(sample,'sphere')
%     idx=((finaldata(:,1)-xc(1)).^2+(finaldata(:,2)-xc(2)).^2)>(r_cutoff*1e9)^2;
%     finaldata(idx,:)=[];
%     I_sel(:,:,idx)=[]; 
% end

%selecting a defined azimuthal section of the sphere 
y_tmp=finaldata(:,2)-xc(2); 
x_tmp=finaldata(:,1)-xc(1); 
phi=atan2(y_tmp,x_tmp)/pi*180; %in deg
idx=abs(mod(phi-phi_c+180,360)-180)>dphi;
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%delting entries with too high fit-errors (residuals)
idx=finaldata(:,10)>qual_thresh;%
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

%filtering data: photon threshold
idx=(photon_hiBND<finaldata(:,4));
finaldata(idx,:)=[];
idx=(finaldata(:,4)<photon_loBND);
finaldata(idx,:)=[];
I_sel(:,:,idx)=[]; 

filtdata=finaldata;

% deleting NaNs from filtdata
tmp=isnan(filtdata);
[tmp_row,tmp_col]=find(tmp);
filtdata(tmp_row,:)=[];

mean_signal=mean(filtdata(:,4));
mean_bg=mean(filtdata(:,5));

%---calculating cramer-rao lower bounds---
if exist('PSF')
     
    [CRBx,CRBy,CRBz]=fun_CRB(PSF,mean_signal,mean_bg,cam.gain,0);
    disp('------------------');
    disp('Cramer-Rao bounds:');
    disp(['min. sqrt(CRLB-x): ' num2str(min(sqrt(CRBx(:))),3)]);
    disp(['min. sqrt(CRLB-y): ' num2str(min(sqrt(CRBy(:))),3)]);
    disp(['min. sqrt(CRLB-z): ' num2str(min(sqrt(CRBz(:))),3)]);
    
end

z_vec=(0:PSF.Nz-1)*PSF.uz;

%z_target=interp1(r_coord,z_theory,r_mol)*1e9; %traget z-value in nm

%----3D scatterplot: show all filtered data----
figure(1); 
%plot3(locdata_TS(:,1)/1e3,locdata_TS(:,2)/1e3,locdata_TS(:,3),'o'); grid on;
markersize=1; %(locdata_TS(:,4)+1)/1e3;
markercolor=filtdata(:,3); %(:,8)
scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,filtdata(:,3),markersize,markercolor); grid on;
colormap jet; 
colorbar; 
title('local. data, thunderstorm');
zlabel('z / nm');
xlabel('x / µm');
ylabel('y / µm');



h=figure(6); 
markersize=1;
markercolor=filtdata(:,3);% markercolor=uint8(markercolor);
markersize=2;
scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,markersize,markercolor);
xlabel('x / µm');
ylabel('y / µm');
axis equal; axis tight; 
colormap jet; 
colorbar; 
title('loc. data, z-values color-coded');
% %
set(gca,'Color','w');
colormap parula; 

% figure(2);
% z_dev=filtdata(:,3)-(r_sphere-sqrt(r_sphere.^2-r_mol.^2))*1e9; %deviation of est. z-pos. from the sphere surface
% markercolor=z_dev;
% scatter3(filtdata(:,1)/1e3,filtdata(:,2)/1e3,z_dev,markersize,markercolor); grid on;
% colormap jet; 
% colorbar; 
% title('deviation from sphere surface');
% zlabel('z / nm');
% xlabel('x / µm');
% ylabel('y / µm');

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

%% ----2D scatterplot with color-coded z-value----


% 
% %---2D-scatterplot; signal-photon-no. color-coded---
% figure(7);
% markercolor=filtdata(:,4); 
% %markercolor=uint8(markercolor/max(markercolor(:))*255*2);
% scatter(filtdata(:,1)/1e3,filtdata(:,2)/1e3,2,markercolor);
% xlabel('x / µm');
% ylabel('y / µm');
% axis equal; axis tight; 
% colormap jet; 
% colorbar; 
% title('loc. data, photon-counts color-coded');


%% show details to the CURSOR-selected molecule image
%slelect molecule in image first, then execute this block

h=figure(6);
dcm_obj = datacursormode(h);
c_info = getCursorInfo(dcm_obj); 

XData=c_info.Target.XData;
YData=c_info.Target.XData;

selected_loc_no=find([XData' YData']==[c_info.Position(1:2)]); %c_info.Position(1))

figure(33)
imagesc(I_sel(:,:,selected_loc_no)); axis equal; axis tight; 
title(['z=' num2str(filtdata(selected_loc_no,3),4)]);

%% plotting only "heroes", i.e. molecules that are "on" for serveral consecutive frames

%clear hero
minlength=6; %minimum frame-no a hero is on
bridgelength=2; %length of allowed "dark" frames in between on-frames
binsize=80; %side length of spatial bins in nm
  
batchsize=20e3; 
hero=[];
for mm=1:floor(length(filtdata)/batchsize) %for long data sets: evaluating in batches 
    
    disp(['processing batch no. ' num2str(mm)]); 
    idx0=(mm-1)*batchsize+1; 
    tmp=fun_findheroes(filtdata(idx0:(idx0+batchsize),:),minlength,bridgelength,binsize); %returns filtdata of heroes in a structure array 
    hero=[hero,tmp];
end
    
    no_heroes=length(hero);

    for m=1:no_heroes
            hero(m).mean=mean(hero(m).locdata);
            hero(m).std=std(hero(m).locdata);
            hero(m).lengths=length(hero(m).locdata);
    end

    hero_means=vertcat(hero.mean);
    hero_sigmas=vertcat(hero.std);
    hero_lengths=vertcat(hero.lengths);

    % plotting 3D positions of heroes
    figure(11); 
    colormap jet; 
    markercolor=hero_sigmas(:,3);
    scatter3(hero_means(:,1),hero_means(:,2),hero_means(:,3),10,markercolor);
    title(['mean 3D positions of ' num2str(no_heroes) ' heros; minlength=' num2str(minlength) '; bridgelength=' num2str(bridgelength)]);
    h = colorbar;
    ylabel(h, '\sigma_z / nm');

    %calculating CRLBs for heroes
    mean_signal_heroes=mean(hero_means(:,4));
    mean_bg_heroes=mean(hero_means(:,5));
    [CRBx_h,CRBy_h,CRBz_h]=fun_CRB(PSF3D,ux,z(2)-z(1),mean_signal_heroes,mean_bg_heroes,cam.gain,0);
    CRBz_h_extrapol=interp1(z_vec,[CRBz_h CRBz_h(end)],z_theory); 

%hero(13)=[];
    
%% --- plotting z-mean vs. z-sigmas

hero_sigmas_x_normed=hero_sigmas(:,1)./sqrt(hero_means(:,4)/mean_signal_heroes); %normalized standard deviations
hero_sigmas_y_normed=hero_sigmas(:,2)./sqrt(hero_means(:,4)/mean_signal_heroes); %normalized standard deviations
hero_sigmas_z_normed=hero_sigmas(:,3)./sqrt(hero_means(:,4)/mean_signal_heroes); %normalized standard deviations

if strcmp(sample,'sphere') % if sample is sphere
    %finding best sphere center and plotting z-positions versus radial
    %coordinate

    xc_heroes=fun_find_sphere_center(xc0,r_coord,z_theory,hero_means); %finding sphere center based on hero mean positions
    xc_heroes=[9.1 12.4]*1e3; %change sphere-center manually

    figure(12); 
    subplot(2,1,1);
    r_heroes=sqrt((hero_means(:,1)-xc_heroes(1)).^2+(hero_means(:,2)-xc_heroes(2)).^2);
    errorbar(r_heroes*1e-3,hero_means(:,3),hero_sigmas_z_normed,'.');
    hold on; 
    plot(r_coord*1e6,z_theory*1e9,'r')
    plot(r_coord*1e6,z_theory*1e9+[sqrt(CRBz_h_extrapol)],'r.');
    plot(r_coord*1e6,z_theory*1e9-[sqrt(CRBz_h_extrapol)],'r.');
    hold off; 
    ylabel('z /nm');
    xlabel('r /µm');
    xlim([0 20]); 
    grid on; 
    title(['z-pos. of ' num2str(no_heroes) ' heroes; minlength=' num2str(minlength) '; bridgelength=' num2str(bridgelength) ' xc=' num2str(xc_heroes/1000,4)]);
    subplot(2,1,2);
    plot(r_heroes*1e-3,hero_sigmas(:,3),'.');
    hold on; 
    plot(r_coord*1e6,[0 sqrt(CRBz_h)],'r.-');
    %hold off; 
    xlabel('radial distance /µm');
    ylabel('\sigma_z');
    xlim([0 20]); 
else
    figure(12);
    %color=hero_means(:,8);%-min(hero_means(:,8));
    color=hero_lengths; 
    color=[0 0 0]
    %----scatterplot----
    scatter(hero_means(:,3),hero_sigmas_z_normed,[],color,'.');
    ylim([0 30]); xlim([0 220]);
    %----with x-errorbars---
%     errorbar(hero_means(:,3),hero_sigmas_normed,hero_sigmas(:,3)./sqrt(hero_lengths),'horizontal','.'); 
%     ylim([0 50]); xlim([0 200]);
%     
    xlabel('mean z-value / nm');
    ylabel('\sigma_z / nm'); grid on;
    title(['normalized \sigma_z of heroes, binsze=' num2str(binsize) ' nm, minleng.=' num2str(minlength) '; bridgeleng.=' num2str(bridgelength)]);
    hold on; 
    plot(linspace(0,250,length(CRBz_h)),sqrt(CRBz_h)+0,'--');
    hold off; 
    colormap jet; 
    ylim([0 30]); xlim([0 220]);
    
    %plotting x-y-sigmas
    figure(13);

    %plot(hero_means(:,3),sqrt(hero_sigmas_x_normed.^2+hero_sigmas_y_normed.^2),'r.'); 
    
    plot(hero_means(:,3),hero_sigmas_x_normed,'r.'); 
    hold on; 
    plot(hero_means(:,3),hero_sigmas_y_normed,'b.'); 
    plot(linspace(0,250,length(CRBx_h)),sqrt(CRBx_h)+0,'--');
    hold off; 
    ylabel('\sigma_{x,y} / nm');
    xlabel('mean z-value / nm');
    grid on;   
    title(['normalized \sigma_xy of heroes, binsze=' num2str(binsize) ' nm, minleng.=' num2str(minlength) '; bridgeleng.=' num2str(bridgelength)]);
    ylim([0 55]); xlim([0 220]);
end


  
%% plot time traces of heroes 
%close fig(15)

m=5; 
meancorr=1;
%for m=1:no_heroes
    markersize=3; 
    figure(13);
    subplot(3,1,1);
    scatter(hero(m).locdata(:,8),hero(m).locdata(:,1)-meancorr*hero_means(m,1),markersize); ylabel('x /nm');
    %hold on;  
    title(['hero no.' num2str(m)]);
    grid on; 
    subplot(3,1,2);
    scatter(hero(m).locdata(:,8),hero(m).locdata(:,2)-meancorr*hero_means(m,2),markersize); ylabel('y /nm');
    %hold on;
    grid on; 
    subplot(3,1,3);
    scatter(hero(m).locdata(:,8),hero(m).locdata(:,3)-0*meancorr*hero_means(m,3),markersize); ylabel('z /nm');
    %hold on;
    xlabel('frame no.');
    grid on; 
   
    figure(14);
    %markercolor=hero(m).locata(:,8); %frame no is color
    plot3(hero(m).locdata(:,1),hero(m).locdata(:,2),hero(m).locdata(:,3),'.-');
    grid on; 
    xlabel('x/nm'); ylabel('y/nm'); zlabel('z/nm');
    axis equal; axis tight; 
    title(['hero no.' num2str(m)]);
   
%end

%---finding molecules images that belong to heroes----: 
loc_numbers=hero(m).locdata(:,7); %loc.numbers of molecule-images that belong to to hero no. m
herolength=length(loc_numbers);
for q=1:herolength
    idx=find(filtdata(:,7)==loc_numbers(q));
    hero(m).images(:,:,q)=I_sel(:,:,idx);
    %imagesc(hero(m).images(:,:,q)); pause; colormap gray; 
end

frames=hero(m).locdata(:,8);

figure(15); 
rowno=round(sqrt(herolength));
colno=ceil(herolength/rowno);
for q=1:herolength
    subplot(rowno,colno,q);
    imagesc(hero(m).images(:,:,q)); 
    title(['frame no. ' num2str(frames(q))]);
end

% plotting 3D positions of heroes
figure(11); 
scatter3(hero_means(:,1),hero_means(:,2),hero_means(:,3),10);
title(['mean 3D positions of ' num2str(no_heroes) ' heros; minlength=' num2str(minlength) '; bridgelength=' num2str(bridgelength)]);
hold on; 
scatter3(hero(m).mean(1),hero(m).mean(2),hero(m).mean(3),10,'r');
hold off; 
axis equal; axis tight; 

hero(m).std(3)


%% -----for sphere: polynomial fit through data to provide better precision estimate-----
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
boxrad=200; %in nm; side length of squared selection "box" around selected point

for m=1:length(px)
    idx=find(logical(abs(filtdata(:,1)-px(m)*1e3)<=boxrad) & logical(abs(filtdata(:,2)-py(m)*1e3)<=boxrad)); %finding indices of localizations that are within the selected area
    
    h=figure(9);
    markersize=3;
    markercolor=filtdata(idx,3);% markercolor=uint8(markercolor);

    scatter3(filtdata(idx,1),filtdata(idx,2),filtdata(idx,3),markersize,markercolor);
    xlabel('x / nm');
    ylabel('y / nm');
    zlabel('z / nm');
    x_mean(m)=mean(filtdata(idx,1)); 
    x_std(m)=std(filtdata(idx,1));
    y_mean(m)=mean(filtdata(idx,2)); 
    y_std(m)=std(filtdata(idx,2));
    z_mean(m)=mean(filtdata(idx,3)); 
    z_std(m)=std(filtdata(idx,3));
    disp('----------------');
    disp(['x-pos (mean/std) = ' num2str(x_mean(m)/1e3,3) 'µm /' num2str(x_std(m),3) 'nm']);
    disp(['y-pos (mean/std) = ' num2str(y_mean(m)/1e3,3) 'µm /' num2str(y_std(m),3) 'nm']);
    disp(['z-pos (mean/std) = ' num2str(z_mean(m),3) 'nm /' num2str(z_std(m),3) 'nm']);
    colormap parula;  
    
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

%%

figure(9); 
ax=gca; 
colormap parula; 
set(gca,'Color','k');
set(gca,'GridColor','w');

set(gca,'xcolor','w'); 
set(gca,'ycolor','w'); 
set(gca,'zcolor','w'); 
set(gca,'Linewidth',1); 
xlim([4500 5500]);
ylim([5000 5600]); 

view([0,0])

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
xlabel('radial coord / µm');
ylabel('z /  nm');
title([num2str(size(filtdata,1)) ' loc., spher-cent=' num2str(xc(1)/1000,3) '/' num2str(xc(2)/1000) ' µm, std=' num2str(delta_z3,2) 'nm, qual-th.=' num2str(qual_thresh) ', drift corr']);
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
xlabel('x / µm');
ylabel('y / µm');
zlim([0 250]);

%%
