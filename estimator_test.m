%program genereates test molecule images at different positions tests the
%estimator functions

clear all; 
close all; 

%% user parameters

two_ch='y'; %two imaging channels?

sig=5e3; %signal per molecule in the topimage in photons (refers to signal of brightest molecule)
BG1=100; %background level in photons in top image
BG2=100; %optional: background level in photons in bottom image

%camera parameters
gain=1; 
amp=9.9; %set to 9.9 for andor; electrons per count 
QE=0.95; %quantum efficiency

if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
    sig1=sig/2; 
    sig2=sig/2; 
else
    sig1=sig; 
end

%PSF-path
PSFpath='C:\Users\q004aj\Desktop\PSFs\';
%PSF for top image stack: 
name_PSF_top='PSF5D_13x13_0-2-250nm_RI=1.45_dz=0_aberrfree_os3.mat';
%PSF for bottom image stack: 
name_PSF_bottom='PSF5D_13x13_0-2-250nm_RI=1.45_dz=-500_aberrfree_os3.mat';


%% loading PSFs

load([PSFpath name_PSF_top]);
[nx0,ny0,nz0,nxi,nyi]=size(PSF5D);
PSF1=PSF5D; 
PSF1_s=PSF1(:,:,:,ceil((nxi+1)/2),ceil((nxi+1)/2));

ux0=ux*os;
uz=z_vec(2)-z_vec(1);

if exist('sigma')
    z_ini_info1=sigma/os; %if the loaded PSF-model contains a vector called "sigma"; sigma represents a gaussian-width-versus-z_ini-curve which provides a good initial z-estimate
else  %otherwise, take the center-value of the entire z-range as initial estimate for z
    z_ini_info1=round(nz0/2); %constant initial z-estimate (frame-no of PSF-stack)
end
disp('PSF1 loaded');  

if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
   load([PSFpath name_PSF_bottom]); %loading 5D PSF-model
   PSF2=PSF5D; 
   PSF2_s=PSF2(:,:,:,ceil((nxi+1)/2),ceil((nxi+1)/2)); %only central-pixel-PSF
   PSF_tot=[PSF1_s PSF2_s]; %concatenate PSFs

   if exist('sigma')
        z_ini_info2=sigma/os; %if the loaded PSF-model contains a vector called "sigma"; sigma represents a gaussian-width-versus-z_ini-curve which provides a good initial z-estimate
   else  %otherwise, take the center-value of the entire z-range as initial estimate for z
        z_ini_info2=round(nz0/2); %constant initial z-estimate (frame-no of PSF-stack)
   end
   disp('PSF2 loaded');  
end
   
clear PSF5D;

%create coordinate system 
x=(1:nx0)-ceil((nx0+1)/2);
y=x; 
[X,Y]=ndgrid(x,y);
clear x_data;
x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx, ny]=size(X); %size of molecule image

%% ----calculating Cramer-Rao lower bounds----

if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
    %calculate theoretical limits: 
    [CRBx_test,CRBy_test,CRBz_test,CRBsig_test,CRBbg_test,FI_test]=fun_CRB(PSF_tot./repmat(trapz(trapz(PSF_tot,1),2),[nx0 2*ny0 1]),ux0,uz,sig1+sig2,(BG1+BG2)/2,gain);
    
    %alternatively, adding the Fisher-information of both channels, then calculate the joint CRB: 
    [~,~,~,~,~,FI1]=fun_CRB(PSF1_s./repmat(trapz(trapz(PSF1_s,1),2),[nx0 ny0 1]),ux0,uz,sig1,BG1,gain);
    [~,~,~,~,~,FI2]=fun_CRB(PSF2_s./repmat(trapz(trapz(PSF2_s,1),2),[nx0 ny0 1]),ux0,uz,sig2,BG2,gain);
    FI=FI1+FI2; 
    for mm=1:nz0-1;
        FI_inv=inv(FI(:,:,mm));
        CRBx(mm)=FI_inv(1,1);
        CRBy(mm)=FI_inv(2,2);
        CRBz(mm)=FI_inv(3,3);
        CRBsig(mm)=FI_inv(4,4);
        CRBbg(mm)=FI_inv(5,5);
    end
else
    [CRBx,CRBy,CRBz,CRBsig,CRBbg,FI]=fun_CRB(PSF1_s./repmat(trapz(trapz(PSF1_s,1),2),[nx0 ny0 1]),ux0,uz,sig1,BG1,gain);
end

%------ plotting CRBs -------------
figure(1);
plot(sqrt(CRBx_test),'b'); hold on; 
plot(sqrt(CRBy_test),'b-'); 
plot(sqrt(CRBz_test),'b')
plot(sqrt(CRBx),'r'); hold on; 
plot(sqrt(CRBy),'r-'); 
plot(sqrt(CRBz),'r')

disp('done');

hold off; 

%% creating test_images with defined signal and background-levels

no_images=50; 
x_pos=0*ones(1,no_images); %integer number
y_pos=0*ones(1,no_images);  %integer number
z_pos=100*ones(1,no_images); %integer number

%-------calculating images---------------
I1_nf=zeros(nx0,ny0,no_images);
I1=I1_nf; 
for m=1:no_images %creating images
    I1_nf(:,:,m)=(PSF1(:,:,z_pos(m),x_pos(m)+ceil((nxi+1)/2),y_pos(m)+ceil((nxi+1)/2))); %noisefree
    I1(:,:,m)=poissrnd(I1_nf(:,:,m)*sig1+BG1); %adding background and noise
end

if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
    I2_nf=zeros(nx0,ny0,no_images);
    I2=I2_nf; 
    for m=1:no_images %creating images
        I2_nf(:,:,m)=(PSF2(:,:,z_pos(m),x_pos(m)+ceil((nxi+1)/2),y_pos(m)+ceil((nxi+1)/2))); %noisefree
        I2(:,:,m)=poissrnd(I2_nf(:,:,m)*sig1+BG1); %adding background and noise
    end    
end
%----------------------------------------
no_images=size(I1,3);
disp('test images created.');

%% run estimator

resnorm_thresh=inf; 
showimage=1;
clear est x_est y_est z_est N_est BG_est

figure(2); 
for m=1:no_images
    if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
        ratio=1; %brightness ratio between images of both channels
        est(m,:)=fun_MLE_2channels(I1(:,:,m),I2(:,:,m),PSF1,PSF2,interp_incr,x_data,[z_ini_info1; z_ini_info2],resnorm_thresh,ratio,showimage);
        %format: est=[x1_est y1_est x2_est y2_est z_est N_est BG1_est BG2_est resnorm1 resnorm2];
        x_est(m)=est(m,1)*ux0*1e9; 
        y_est(m)=est(m,2)*ux0*1e9; 
        z_est(m)=est(m,5)*uz*1e9; 
        N_est(m)=est(m,6)/gain*amp/QE;
        BG1_est(m)=est(m,7)/gain*amp/QE;
        BG2_est(m)=est(m,8)/gain*amp/QE;
    else
        est(m,:)=fun_MLE(I1(:,:,m),PSF1,interp_incr,x_data,z_ini_info1,resnorm_thresh,showimage);
        x_est(m)=est(m,1)*ux0*1e9; 
        y_est(m)=est(m,2)*ux0*1e9; 
        z_est(m)=est(m,3)*uz*1e9; 
        N_est(m)=est(m,4)/gain*amp/QE;
        BG_est(m)=est(m,5)/gain*amp/QE;
    end
      
end

disp('all data estimated.');

%% ------------------------ show results--------------------;

figure(1);

subplot(3,1,1);
plot(x_est,'.'); xlabel('image no.'); ylabel('x_{est} /nm'); 
hold on; 
plot(mean(x_est)+sqrt(CRBx(z_pos)),'r--'); plot(mean(x_est)-sqrt(CRBx(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_x=' num2str(std(x_est),3) ', CRBx=' num2str(mean(sqrt(CRBx)),3)]); 

subplot(3,1,2);
plot(y_est,'.'); xlabel('image no.'); ylabel('y_{est} /nm');
hold on; 
plot(mean(y_est)+sqrt(CRBy(z_pos)),'r--'); plot(mean(y_est)-sqrt(CRBy(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_y=' num2str(std(y_est),3) ', CRBy=' num2str(mean(sqrt(CRBy)),3)]); 

subplot(3,1,3);
plot(z_est,'.'); xlabel('image no.'); ylabel('z_{est} /nm');
hold on; 
plot(mean(z_est)+sqrt(CRBz(z_pos)),'r--'); plot(mean(z_est)-sqrt(CRBz(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_z=' num2str(std(z_vec(z_pos)*1e9-z_est),3) ', CRBz=' num2str(mean(sqrt(CRBz)),3)]); 



