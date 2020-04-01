%program genereates test molecule images at different positions and tests the
%estimator functions

clear all; 
close all; 

%% user parameters

two_ch='y'; %two imaging channels?

sig=2e3; %signal per molecule in the topimage in photons (refers to signal of brightest molecule)
BG1=100; %background level in photons in top image
BG2=100; %optional: background level in photons in bottom image

%camera parameters
gain=1; 
amp=1; %set to 9.9 for andor; electrons per count 
QE=1; %quantum efficiency

if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
    sig1=sig/2; 
    sig2=sig/2; 
else
    sig1=sig; 
end

%PSF-path
PSFpath='C:\Users\q004aj\Desktop\PSFs\NA1,7\';
%PSF for top image stack: 
name_PSF_top='PSF5D_NA1.67_13x13_0-2-250nm_RI=1.33_dz=0_aberrfree_os3.mat';
%PSF for bottom image stack: 
name_PSF_bottom='PSF5D_UAF_NA1.67_13x13_0-2-250nm_RI=1.33_dz=0_aberrfree_os3.mat';

disp('done');


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

    %CRB of first channel: 
    [CRBx_1,CRBy_1,CRBz_1,CRBsig_1,CRBbg_1,FI1]=fun_CRB(PSF1_s,ux0,uz,sig1,BG1/2,gain,0);
    %CRB of second channel: 
    [CRBx_2,CRBy_2,CRBz_2,CRBsig_2,CRBbg_2,FI2]=fun_CRB(PSF2_s,ux0,uz,sig2,BG2/2,gain,0);

    %joint CRB for both channels: 
    %[CRBx,CRBy,CRBz,CRBsig,CRBbg,FI]=fun_CRB(PSF_tot,ux0,uz,sig1+sig2,(BG1+BG2)/2,gain,0);
    
    %alternatively, adding the Fisher-information of both channels, then calculate the joint CRB: 
    %(this seems to provide a worse match to the simulated precisions)
    FI_test=FI1+FI2; 
    for mm=1:nz0-1
        FI_inv=inv(FI_test(:,:,mm));
        CRBx(mm)=FI_inv(1,1);
        CRBy(mm)=FI_inv(2,2);
        CRBz(mm)=FI_inv(3,3);
        CRBsig(mm)=FI_inv(4,4);
        CRBbg(mm)=FI_inv(5,5);
    end
else
    [CRBx,CRBy,CRBz,CRBsig,CRBbg,FI]=fun_CRB(PSF1_s,ux0,uz,sig1,BG1,gain,0);

end

%------ plotting CRBs -------------
% figure(1);
% plot(sqrt(CRBx_1),'g--'); hold on; 
% plot(sqrt(CRBy_1),'g--'); 
% plot(sqrt(CRBz_1),'g')

% plot(sqrt(CRBx_2),'b--'); hold on; 
% plot(sqrt(CRBy_2),'b--'); 
% plot(sqrt(CRBz_2),'b')

% plot(sqrt(CRBx_test),'r--'); hold on; 
% plot(sqrt(CRBy_test),'r--'); 
%plot(sqrt(CRBz_test),'r'); hold on; 
%plot(sqrt(CRBsig_test),'r');
%plot(sqrt(CRBbg_test),'r');

plot(sqrt(CRBx),'cy--'); hold on; 
plot(sqrt(CRBy),'cy--'); 
plot(sqrt(CRBz),'cy');
%plot(sqrt(CRBsig),'cy');
%plot(sqrt(CRBbg),'cy');
ylabel('\sigma /nm'); 
xlabel('z'); 
hold off; 

disp('done');

%% creating test_images with defined signal and background-levels

no_zsteps=10; 
no_avg=50; %no of images per z-step
no_images=no_zsteps*no_avg; 
%---------------
z_vals=2:10:100; 
%z_vals=80; 
x_vals=zeros(length(z_vals));
y_vals=x_vals; 

tmp=repmat(z_vals,[no_avg, 1]);
z_pos=reshape(tmp,1,numel(tmp));

x_pos=0*ones(1,length(z_pos));  %integer number
y_pos=0*ones(1,length(z_pos));  %integer number

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
clear est x_est y_est z_est N_est BG1_est BG2_est

figure(2); 
for m=1:no_images
    if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
        
        %----VARIANT A: JOINT TWO-CHANNEL ESTIMATION---------
%         ratio=1; %brightness ratio between images of both channels
%         %definition of outputs: est=[x1_est y1_est x2_est y2_est z_est N_est BG1_est BG2_est resnorm1 resnorm2];
%         tmp=fun_MLE_2channels(I1(:,:,m),I2(:,:,m),PSF1,PSF2,interp_incr,x_data,[z_ini_info1; z_ini_info2],resnorm_thresh,ratio,showimage);
%         if isempty(tmp)
%             est(m,:)=0;
%         else
%             est(m,:)=tmp; 
%         end
%         %format: est=[x1_est y1_est x2_est y2_est z_est N_est BG1_est BG2_est resnorm1 resnorm2];
%         x1_est(m)=est(m,1)*ux0*1e9; 
%         y1_est(m)=est(m,2)*ux0*1e9;
%         x2_est(m)=est(m,3)*ux0*1e9; 
%         y2_est(m)=est(m,4)*ux0*1e9;
%         z_est(m)=est(m,5)*uz*1e9; 
%         N_est(m)=est(m,6)/gain*amp/QE;
%         BG1_est(m)=est(m,7)/gain*amp/QE;
%         BG2_est(m)=est(m,8)/gain*amp/QE;
% 
%         %combine x-y estimates of both channels to a single, joint estimate of higher precision: 
%         idx=max(min(round(est(m,5))+1,nz0-1),0); 
%         x_est(m)=(x1_est(m)/CRBx_1(idx)+x2_est(m)/CRBx_2(idx))/(1/CRBx_1(idx)+1/CRBx_2(idx));
%         y_est(m)=(y1_est(m)/CRBy_1(idx)+y2_est(m)/CRBy_2(idx))/(1/CRBy_1(idx)+1/CRBy_2(idx));

        
        %----VARIANT B: COMBINING TWO SINGLE-CHANNEL ESTIMATES----------
        
        est1(m,:)=fun_MLE(I1(:,:,m),PSF1,interp_incr,x_data,z_ini_info1,resnorm_thresh,showimage);
        x1_est(m)=est1(m,1)*ux0*1e9; 
        y1_est(m)=est1(m,2)*ux0*1e9; 
        z1_est(m)=est1(m,3)*uz*1e9; 
        N1_est(m)=est1(m,4)/gain*amp/QE;
        BG1_est(m)=est1(m,5)/gain*amp/QE;

        est2(m,:)=fun_MLE(I2(:,:,m),PSF2,interp_incr,x_data,z_ini_info2,resnorm_thresh,showimage);
        x2_est(m)=est2(m,1)*ux0*1e9; 
        y2_est(m)=est2(m,2)*ux0*1e9; 
        z2_est(m)=est2(m,3)*uz*1e9; 
        N2_est(m)=est2(m,4)/gain*amp/QE;
        BG2_est(m)=est2(m,5)/gain*amp/QE;
        
        %combine estimates of both channels to a single, joint estimate of higher precision: 
        idx=z_pos(1); %we use for now the ground truth
        
        x_est(m)=(x1_est(m)/CRBx_1(idx)+x2_est(m)/CRBx_2(idx))/(1/CRBx_1(idx)+1/CRBx_2(idx));
        y_est(m)=(y1_est(m)/CRBy_1(idx)+y2_est(m)/CRBy_2(idx))/(1/CRBy_1(idx)+1/CRBy_2(idx));       
        z_est(m)=(z1_est(m)/CRBz_1(idx)+z2_est(m)/CRBz_2(idx))/(1/CRBz_1(idx)+1/CRBz_2(idx));
        N_est(m)=(1+1/ratio)*(N1_est(m)/CRBsig_1(idx)+N2_est(m)/CRBsig_2(idx))/(1/CRBsig_1(idx)+1/CRBsig_2(idx));
                
    else
       
        %----ONE-CHANNEL ESTIMATION---------
        est(m,:)=fun_MLE(I1(:,:,m),PSF1,interp_incr,x_data,z_ini_info1,resnorm_thresh,showimage);
        x_est(m)=est(m,1)*ux0*1e9; 
        y_est(m)=est(m,2)*ux0*1e9; 
        z_est(m)=est(m,3)*uz*1e9; 
        N_est(m)=est(m,4)/gain*amp/QE;
        BG1_est(m)=est(m,5)/gain*amp/QE;
    end
      
end

disp('all data estimated.');

% ------------------------ show results--------------------;

figure(1);

subplot(3,2,1);
plot(x_est,'.'); xlabel('image no.'); ylabel('x_{est} /nm'); 
hold on; 
plot(mean(x_est)+sqrt(CRBx(z_pos)),'r--'); plot(mean(x_est)-sqrt(CRBx(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_x=' num2str(std(x_est),3) ', CRBx=' num2str(sqrt(CRBx(z_pos(1))),3)]); 

subplot(3,2,2);
plot(y_est,'.'); xlabel('image no.'); ylabel('y_{est} /nm');
hold on; 
plot(mean(y_est)+sqrt(CRBy(z_pos)),'r--'); plot(mean(y_est)-sqrt(CRBy(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_y=' num2str(std(y_est),3) ', CRBy=' num2str(sqrt(CRBy(z_pos(1))),3)]); 

subplot(3,2,3);
plot(z_est,'.'); xlabel('image no.'); ylabel('z_{est} /nm');
hold on; 
plot((z_pos)*uz*1e9+sqrt(CRBz(z_pos)),'r--'); plot(z_pos*uz*1e9-sqrt(CRBz(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_z=' num2str(std(z_vec(z_pos)*1e9-z_est),3) ', CRBz=' num2str(sqrt(CRBz(z_pos(1))),3)]); 

subplot(3,2,4);
plot(N_est,'.'); xlabel('image no.'); ylabel('N_{est} /nm');
hold on; 
plot(sig+sqrt(CRBsig(z_pos)),'r--'); plot(sig-sqrt(CRBsig(z_pos)),'r--');
hold off; 
xlim([1,no_images]);
title(['\sigma_N=' num2str(std(N_est),3) ', CRBsig=' num2str(sqrt(CRBsig(z_pos(1))),3)]); 


    if strcmp(two_ch,'y') %flag, indicating that two channels are simulated
        subplot(3,2,5);
        plot(BG1_est,'.'); xlabel('image no.'); ylabel('BG1_{est} /nm');
        hold on; 
        plot(BG1+sqrt(CRBbg_1(z_pos)),'r--'); plot(BG1-sqrt(CRBbg_1(z_pos)),'r--');
        hold off; 
        xlim([1,no_images]);
        title(['\sigma_{BG1}=' num2str(std(BG1_est),3) ', CRBbg=' num2str(sqrt(CRBbg_1(z_pos(1))),3)]); 
        subplot(3,2,6);
        
        plot(BG2_est,'.'); xlabel('image no.'); ylabel('BG2_{est} /nm');
        hold on; 
        plot(BG2+sqrt(CRBbg_2(z_pos)),'r--'); plot(BG2-sqrt(CRBbg_2(z_pos)),'r--');
        hold off; 
        xlim([1,no_images]);
        title(['\sigma_{BG2}=' num2str(std(BG2_est),3) ', CRBbg=' num2str(sqrt(CRBbg_2(z_pos(1))),3)]); 
    end
    
    %%

ratio=N1_est./N2_est.*squeeze(sum(sum(PSF1_s,1),2))./squeeze(sum(sum(PSF2_s,1),2)); 
plot(ratio)