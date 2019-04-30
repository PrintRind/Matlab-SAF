%MLE fit to molecule images based on a numerical model of the PSF
%estimated are: x,y,z-position  / no. photons / background

clear all; 
close all; 
clc;

%loading 3D PSF model: [x,y, z-distance from coverslip]
%load('PSF_SAF_2nm_3,52rad.mat') %data exported from "optimal_SAF_PSF.m"; variables: PSF_tot, PSF_UAF, ux, z_vec, RI, NA
%load('PSF_SAF_0-2nm-300nm.mat') %data exported from "optimal_SAF_PSF.m"; variables: PSF_tot, PSF_UAF, ux, z_vec, RI, NA
%load('PSF_SAF_0-2nm-300nm_dz=-430nm.mat') %defocused PSF-model; should give better z-estimates
%load('PSF_SAF_0-2nm-200nm_dz=-430nm.mat') %defocused PSF-model; should give better z-estimates
%load('PSF_SAF_0-2nm-200nm_dz=0nm.mat') %defocused PSF-model; should give better z-estimates
%load('PSF_SAF_0-2nm-300nm_phi=PI.mat') %data exported from "optimal_SAF_PSF.m"; variables: PSF_tot, PSF_UAF, ux, z_vec, RI, NA
%load('PSF_SAF_NA1,49_0-2nm-200nm_dz=0nm.mat') %defocused PSF-model; should give better z-estimates
%load('PSF_SAF_NA1,49_0-2nm-200nm_dz=-430nm.mat') %defocused PSF-model; should give better z-estimates

%vectashield-PSFs:
PSFpath='C:\Users\q004aj\Desktop\PSFs\';

load([PSFpath 'PSF_0-3-250nm_RI=1,45_defoc=-400nm_aberrfree.mat']) %defocused PSF-model; should give better z-estimates
%load([PSFpath 'PSF_0-3-250nm_RI=1,45_defoc=0nm_aberrfree.mat']) %defocused PSF-model; should give better z-estimates

%load('PSF_SAF_0-2nm-200nm_RI=1,45_dz=-750nm_2018-09-26.mat') %defocused PSF-model; should give better z-estimates

%use second image channel?
ch='n';  %choose 'y' or 'n'; if 'y', separate UAF and SAF images are considered. The photon no. is then the total photon number for both images

%<<<<<<<<<<<<
PSF=PSF_tot; %choose PSF: PSF_tot or PSF_UAF
PSF2=PSF_UAF; 
%<<<<<<<<<<<<

%Note: if estimate is "combined", i.e. PSF_tot and PSF_UAF are simultaneously
%recorded and used for the estimate: 
%BG and photon no. 

uz=z_vec(2)-z_vec(1); %z-increment in PSF-stack

%PSF: normalization of energy in each z-slice
energies=(trapz(trapz(PSF,1),2)); 
PSF_norm=PSF;%./repmat(energies,[size(PSF,1),size(PSF,2),1]);
dz_PSF=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
PSF_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF

if strcmp(ch,'y')
    energies2=(trapz(trapz(PSF2,1),2)); 
    PSF2_norm=PSF2;%./repmat(energies2,[size(PSF2,1),size(PSF2,2),1]);
    dz_PSF2=diff(PSF2_norm,1,3); %calculating derivative along z (required for ML estimation)
    PSF2_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF
    ratio_gt=squeeze((energies-energies2)./energies2); %"ground-truth" ratio-curve
end

fw=2; %frame-width; assumed molecule image must be smaller than the PSF-stack images! Ohterwise, NaNs appear when interpolating at different x-y positions
      % frame-width must be larger than any possible x-y-shift; 

      
%create coordinate system 
[nx0, ny0, nz]=size(PSF_norm);
x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
[X,Y]=ndgrid(x,y);

x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx, ny]=size(X); %size of molecule image


%% taking out an arbitraty x-y slice of the PSF-stack: this is our "measurement"

photon_no=500000; %total photon number
gain=1; 
BG=0; %background level in photons

%calculating CRLBs for single-channel imaging (all photons in one channel)
[CRBx,CRBy,CRBz]=fun_CRB(PSF./repmat(sum(sum(PSF,1),2),[size(PSF,1) size(PSF,2) 1]),ux,uz,photon_no,BG,gain);

fun_estimate=@(v) v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3)); %calculates molecule image from estimated parameters x
fun_estimate2=@(v) v(5)+v(4)*interpn(PSF2_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3)); %calculates molecule image from estimated parameters x

BGmask1=zeros(nx,ny);
BGmask2=zeros(nx,ny);

z_truth=1:10:nz; %z-slide no. of PSF_stack which is assumed to be the ground truth
clear z_mean z_sigma x_mean x_sigma y_mean y_sigma N_mean N_sigma BG_mean BG_sigma z_prob N BG_est z_est x_est y_est N_est dx dy var0 N0 N0_2 ratio_mean ratio_sigma
clear z2_mean z2_sigma x2_mean x2_sigma y2_mean y2_sigma N2_mean N2_sigma BG2_mean BG2_sigma

for mm=1:length(z_truth) %stepping through all z-slices and simulating several measurements for each slice
    
    %---creating a simulated measurement: I ---
    dx(mm)=fw*(rand-0.5); %assumed lateral shifts
    dy(mm)=fw*(rand-0.5);
    %I=PSF_norm(1+fw-dx:end-fw-dx,1+fw-dy:end-fw-dy,z_truth(mm)); %simulated image, w/o noise
    I_nf=interpn(PSF_norm,(fw+1:nx+fw)-dx(mm),(fw+1:ny+fw)'-dy(mm),z_truth(mm)); %imagesc(I); %nf...noise-free
    
    %generate second ground-truth image (UAF):
    if strcmp(ch,'y')
        I2_nf=interpn(PSF2_norm,(fw+1:nx+fw)-dx(mm),(fw+1:ny+fw)'-dy(mm),z_truth(mm)); %imagesc(I); %nf...noise-free
        frac_tot=0.5; %half of the energy goes into the tot-image
        frac_UAF=0.5*squeeze(energies2(z_truth)./energies(z_truth)); %a bit less goes into the UAF-channel due to the SAF-block
    else
        frac_tot=1;
    end
    
    
    for m=1:10 %repetitive measurements loop
                
        if strcmp(ch,'y') %for 2-image-channels
            I=poissrnd(I_nf/sum(I_nf(:))*photon_no*frac_tot+BG*frac_tot); %adding noise
            I2=poissrnd(I2_nf/sum(I2_nf(:))*photon_no*frac_UAF(mm)+BG*frac_UAF(mm)); %adding noise
        else
            I=poissrnd(I_nf/sum(I_nf(:))*photon_no+BG); %adding noise
        end
        
        %neg. log-likelihood for a single Pixel i: -LLH=[PSF_i-I_i*log(PSF_i)]
        fun_LLH=@(v) sum(sum(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))-I.*log(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))),1),2);
        fun_LLH2=@(v) sum(sum(v(5)+v(4)*interpn(PSF2_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))-I2.*log(v(5)+v(4)*interpn(PSF2_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))),1),2);
        fun_LLH_joint=@(v) fun_LLH(v(1:5))+fun_LLH2([v(1:3) v(6:7)]); %joint log-likelihood function

        %initial BG estimate for I1 (for Gaussfit)
        BGmask1(1,:)=1; BGmask1(end,:)=1; BGmask1(:,1)=1; BGmask1(:,end)=1;
        BGmask2(2,2:end-1)=1; BGmask2(end-1,2:end-1)=1; BGmask2(2:end-1,2)=1; BGmask2(2:end-1,end-1)=1;
        BG2=sum(BGmask2(:).*I(:))/sum(BGmask2(:)); %mean signal in inner frame
        BG1=sum(BGmask1(:).*I(:))/sum(BGmask1(:)); %mean signal in outer frame
        %BG0=abs(BG1-(BG2-BG1)*2);
        BG0=BG1;
        %BG0=BG; %take true BG as initial estimate: this should be the best
       
        %----Gaussfit for initial estimates of x-y-position,BG,N0-----
        param_ini=[BG0 max(I(:))-BG0 0 0 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
        lb=[0 0.5*(max(I(:))-BG0) -nx/3 -ny/3 0.5]; %lower bounds for fit-parameters
        ub=[2*BG0 1.5*(max(I(:))-BG0) nx/3 ny/3 3]; %upper bounds
        [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_gauss_and_offset_test,param_ini,x_data,I,lb,ub);
        Gaussfit=fun_gauss_and_offset_test(Param,x_data);
        x0=Param(3);
        y0=Param(4);
        BG0=Param(1); %this serves as initial estimate for the MLE fit
        %N0=Param(2)*2*pi*Param(5)^2; %energy contained (see e.g. formula in Thunderstorm script)
        N0=sum(I(:))-BG0*nx*ny; %initial guess for number of photons - seems to work better 
                
        %----MLE estimation----------        
        gauss_est(m,:)=[x0 y0 round(nz/2) N0 BG0]; %initial estimates (e.g. from Gaussfit)
        LB=[-fw/2 -fw/2 1 0.5*N0 0.5*BG0]; %lower bounds
        UB=[+fw/2 +fw/2 nz 2*N0 2*BG0]; %upper bounds
        %tmp=fmincon(fun_LLH,x0,[],[],[],[],[1 0.75*N_ini 0.75*BG_ini],[nz 1.5*N_ini 1*BG_ini]);
        tmp=fminsearchbnd(fun_LLH,gauss_est(m,:),LB,UB);
        x_est(m)=tmp(1);
        y_est(m)=tmp(2);
        z_est(m)=tmp(3)
        N_est(m)=tmp(4);
        BG_est(m)=tmp(5);
        
        if strcmp(ch,'y');        %if second image channel is used
            %initial BG estimate for second image:
            BG2=sum(BGmask2(:).*I2(:))/sum(BGmask2(:)); %mean signal in inner frame
            BG1=sum(BGmask1(:).*I2(:))/sum(BGmask1(:)); %mean signal in outer frame
            %BG0=abs(BG1-(BG2-BG1)*2);
            BG0=BG1;
            
            %----2nd image: Gaussfit for initial estimates of x-y-position,BG,N0-----
            param2_ini=[BG0 max(I2(:))-BG0 0 0 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
            lb2=[0 0.5*(max(I2(:))-BG0) -nx/3 -ny/3 0.5]; %lower bounds for fit-parameters
            ub2=[2*BG0 1.5*(max(I2(:))-BG0) nx/3 ny/3 3]; %upper bounds
            [Param2,resnorm2,residual2,exitflag2]=lsqcurvefit(@fun_gauss_and_offset_test,param2_ini,x_data,I2,lb2,ub2);
            Gaussfit2=fun_gauss_and_offset_test(Param2,x_data);
            x0_2=Param2(3);
            y0_2=Param2(4);
            BG0_2=Param2(1); %this serves as initial estimate for the MLE fit
            N0_2=sum(I2(:))-BG0_2*nx*ny; %initial guess for number of photons
            %N0_2=Param2(2)*2*pi*Param2(5)^2; %energy contained (see e.g. formula in Thunderstorm script)
            
%             N0_2=photon_no*frac_UAF(mm); %the exact value - just to test
%             BG0_2=BG/nx/ny*frac_UAF(mm);
            
            %----2nd image: MLE estimation----------   
            gauss_est2(m,:)=[x0_2 y0_2 round(nz/2) N0_2 BG0_2];
            LB2=[-fw/2 -fw/2 1 0.5*N0_2 0.5*BG0_2]; %lower bounds
            UB2=[+fw/2 +fw/2 nz 2*N0_2 2*BG0_2]; %upper bounds
            %tmp=fmincon(fun_LLH,x0,[],[],[],[],[1 0.75*N_ini 0.75*BG_ini],[nz 1.5*N_ini 1*BG_ini]);
            tmp=fminsearchbnd(fun_LLH2,gauss_est2(m,:),LB2,UB2);
            x2_est(m)=tmp(1);
            y2_est(m)=tmp(2);
            z2_est(m)=tmp(3)
            N2_est(m)=tmp(4);
            BG2_est(m)=tmp(5);
            
            
            %---joint MLE-estimation based on both images----
            gauss_est_j(m,:)=[(x0+x0_2)/2 (y0+y0_2)/2 round(nz/2) N0 BG0 N0_2 BG0_2]; %seven parameters for the joint est.
            LB_j=[-fw/2 -fw/2 1 0.5*N0 0.5*BG0 0.5*N0_2 0.5*BG0_2]; %lower bounds
            UB_j=[+fw/2 +fw/2 nz 2*N0 2*BG0 2*N0_2 2*BG0_2]; %upper bounds    
            tmp=fminsearchbnd(fun_LLH_joint,gauss_est_j(m,:),LB_j,UB_j);
            xj_est(m)=tmp(1);
            yj_est(m)=tmp(2);
            zj_est(m)=tmp(3)
            Nj_est(m)=tmp(4);
            BGj_est(m)=tmp(5);
            
        end
  
        
    end
    
                   
    %--MLE estimates from a single image
    I_est=fun_estimate(tmp); %calculate estimated molecule image
    
    z_mean(mm)=interp1(1:length(z_vec),z_vec,mean(z_est))*1e9
    z_sigma(mm)=std(z_est)*uz*1e9
    z_err(mm)=mean(abs(z_est-z_truth(mm))); %mean error
    
    x_mean(mm)=mean(x_est)*ux*1e9;
    x_sigma(mm)=std(x_est)*ux*1e9;
    x_err(mm)=mean(abs((x_est-dx(mm))))*ux*1e9; %mean error
    
    y_mean(mm)=mean(y_est)*ux*1e9;
    y_sigma(mm)=std(y_est)*ux*1e9;
    y_err(mm)=mean(abs((y_est-dy(mm))))*ux*1e9; %mean error
    
    N_mean(mm)=mean(N_est);
    N_sigma(mm)=std(N_est);
    N_err(mm)=mean(abs(N_est-photon_no*frac_tot)); %mean error
    
    BG_mean(mm)=mean(BG_est);
    BG_sigma(mm)=std(BG_est);
    BG_err(mm)=mean(abs(BG_est-BG*frac_tot)); %mean error

    
        if strcmp(ch,'y'); %if a second image-channel is used
            
            %--results from Gauss-fits--
            N_tot=gauss_est(:,4);
            N_UAF=gauss_est2(:,4);
            ratio_mean_gauss(mm)=(mean(N_tot)-mean(N_UAF))/mean(N_UAF); %mean SAF/UAF ratio
            ratio_sigma_gauss(mm)=sqrt(std(N_UAF).^2/mean(N_UAF).^2+std(N_tot).^2/mean(N_UAF).^2).*ratio_mean_gauss(mm); %error propagation

            %--MLE estimates from the second image
            I2_est=fun_estimate(tmp); %calculate estimated molecule image

            z2_mean(mm)=interp1(1:length(z_vec),z_vec,mean(z2_est))*1e9
            z2_sigma(mm)=std(z2_est)*uz*1e9
            z2_err(mm)=mean(abs(z2_est-z_truth(mm))); %mean error

            x2_mean(mm)=mean(x2_est)*ux*1e9;
            x2_sigma(mm)=std(x2_est)*ux*1e9;
            x2_err(mm)=mean(abs((x2_est-dx(mm))))*ux*1e9; %mean error

            y2_mean(mm)=mean(y2_est)*ux*1e9;
            y2_sigma(mm)=std(y2_est)*ux*1e9;
            y2_err(mm)=mean(abs((y2_est-dy(mm))))*ux*1e9; %mean error

            N2_mean(mm)=mean(N2_est);
            N2_sigma(mm)=std(N2_est);
            N2_err(mm)=mean(abs(N2_est-photon_no*frac_UAF(mm))); %mean error

            BG2_mean(mm)=mean(BG2_est);
            BG2_sigma(mm)=std(BG2_est);
            BG2_err(mm)=mean(abs(BG2_est-BG*frac_UAF(mm))); %mean error
                       
            %--results from MLE-estimates
            ratio_mean_MLE(mm)=(N_mean(mm)-N2_mean(mm))/N2_mean(mm); %mean SAF/UAF ratio
            ratio_sigma_MLE(mm)=sqrt(N_sigma(mm)^2/N_mean(mm)^2+N2_sigma(mm)^2/N2_mean(mm)^2).*ratio_mean_MLE(mm); %error propagation
            
            
            %--MLE estimates of the joint images---
            zj_mean(mm)=interp1(1:length(z_vec),z_vec,mean(zj_est))*1e9
            zj_sigma(mm)=std(zj_est)*uz*1e9
            zj_err(mm)=mean(abs(zj_est-z_truth(mm))); %mean error

            xj_mean(mm)=mean(xj_est)*ux*1e9;
            xj_sigma(mm)=std(xj_est)*ux*1e9;
            xj_err(mm)=mean(abs((xj_est-dx(mm))))*ux*1e9; %mean error

            yj_mean(mm)=mean(yj_est)*ux*1e9;
            yj_sigma(mm)=std(yj_est)*ux*1e9;
            yj_err(mm)=mean(abs((yj_est-dy(mm))))*ux*1e9; %mean error
            
            
        end   
           
    disp('');
    disp(['true z=' num2str(z_vec(z_truth)*1e9)])
    
    %visualization
    figure(1); 
    subplot(3,1,1);
    imagesc(I,[0 max(I(:))]); 
    colormap gray;
    axis equal; axis tight; title('measurement'); 
    subplot(3,1,2);
    imagesc(I_est, [0 max(I(:))]);
    colormap gray;
    axis equal; axis tight; title('estimate'); 
    subplot(3,1,3);
    errorbar(z_vec(z_truth(1:mm))*1e9,z_mean,z_sigma,'x');  axis equal, axis tight;
    grid on;
    hold on; 
    plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'r-');
    plot(z_vec*1e9,z_vec*1e9+sqrt([CRBz 0]),'r.');
    plot(z_vec*1e9,z_vec*1e9-sqrt([CRBz 0]),'r.');
    hold off;
    xlabel('z_{truth} / nm');
    ylabel('z_{estimate} / nm');
    pause(0);    

       
end


%% ---visualization of final results

figure(3);

%x-y-estimate
subplot(2,2,1);
errorbar(z_vec(z_truth)*1e9,x_mean,x_sigma,'rx'); 
hold on;
errorbar(z_vec(z_truth)*1e9,y_mean,y_sigma,'gx');  
plot(z_vec(z_truth)*1e9,dx*1e9*ux,'ro');
plot(z_vec(z_truth)*1e9,dy*1e9*ux,'go');
grid on;
hold off;
title(['1st image: MLE of x-y-pos.(x=red); cts=' num2str(photon_no*frac_tot)]);
xlabel('z_{truth} / nm');
ylabel('x/y_{estimate} / nm');

%z-estimate
subplot(2,2,2);
errorbar(z_vec(z_truth)*1e9,z_mean,z_sigma,'x'); 
grid on;
hold on; 
plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'-');  
hold off;
title(['1st image: MLE of z-dist. cts=' num2str(photon_no*frac_tot)]);
xlabel('z_{truth} / nm');
ylabel('z_{estimate} / nm');
xlim([z_vec(z_truth(1)), z_vec(z_truth(end))]*1e9);
ylim([0 z_vec(end)]*1e9);

%photon-number estimate
subplot(2,2,3);
errorbar(z_vec(z_truth)*1e9,N_mean,N_sigma,'-x'); ylabel('est. photon no.'); 
line([z_vec(z_truth(1)),z_vec(z_truth(end))]*1e9,[photon_no photon_no]*frac_tot,'Color','red');
title(['1st image: MLE of photon-no.; cts=' num2str(photon_no*frac_tot)]);
xlabel('z_{truth} / nm');
ylabel('signal photons');

%Background estimate
subplot(2,2,4);
errorbar(z_vec(z_truth)*1e9,BG_mean,BG_sigma,'-x'); ylabel('est. BG');
line([z_vec(z_truth(1)),z_vec(z_truth(end))]*1e9,[BG BG]*frac_tot,'Color','red');
title(['1st image: MLE of BG; cts=' num2str(photon_no*frac_tot)]);
xlabel('z_{truth} / nm');
ylabel('BG photons');

disp('--- mean values over all measurements ----');
disp(' ');
disp(['x_sigma / x_err: ' num2str(mean(x_sigma)) ' / ' num2str(mean(x_err)) ' nm']);
disp(['y_sigma / y_err: ' num2str(mean(y_sigma)) ' / ' num2str(mean(y_err)) ' nm']);
disp(['z_sigma / z_err: ' num2str(mean(z_sigma)) ' / ' num2str(mean(z_err)) ' nm']);
disp(['N_sigma / N_err: ' num2str(mean(N_sigma)) ' / ' num2str(mean(N_err)) ]);
disp(['BG_sigma / BG_err: ' num2str(mean(BG_sigma)) ' / ' num2str(mean(BG_err)) ]);


if strcmp(ch,'y'); 

    figure(4);

    %x-y-estimate
    subplot(2,2,1);
    errorbar(z_vec(z_truth)*1e9,x2_mean,x2_sigma,'rx'); 
    hold on;
    errorbar(z_vec(z_truth)*1e9,y2_mean,y2_sigma,'gx');  
    plot(z_vec(z_truth)*1e9,dx*1e9*ux,'ro');
    plot(z_vec(z_truth)*1e9,dy*1e9*ux,'go');
    grid on;
    hold off;
    title(['2nd image: MLE of x-y-pos.(x=red)']);
    xlabel('z_{truth} / nm');
    ylabel('x/y_{estimate} / nm');

    %z-estimate
    subplot(2,2,2);
    errorbar(z_vec(z_truth)*1e9,z2_mean,z2_sigma,'x'); 
    grid on;
    hold on; 
    plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'-');  
    hold off;
    title(['2nd image: MLE of z-distance']);
    xlabel('z_{truth} / nm');
    ylabel('z_{estimate} / nm');
    xlim([z_vec(z_truth(1)), z_vec(z_truth(end))]*1e9);
    ylim([0 z_vec(end)]*1e9);

    %photon-number estimate
    subplot(2,2,3);
  
    errorbar(z_vec(z_truth)*1e9,N2_mean,N2_sigma,'-x'); ylabel('est. photon no.'); 
    hold on;
    plot(z_vec(z_truth)*1e9,frac_UAF*photon_no,'Color','red');
    hold off;
    title(['2nd image: MLE of photon-no.']);
    xlabel('z_{truth} / nm');
    ylabel('signal photons');

    %Background estimate
    subplot(2,2,4);
    errorbar(z_vec(z_truth)*1e9,BG2_mean*nx*ny,BG2_sigma,'-x'); ylabel('est. BG');
    hold on;
    plot(z_vec(z_truth)*1e9,frac_UAF*photon_no,'Color','red');
    hold off;
    title(['2nd image: MLE of background']);
    xlabel('z_{truth} / nm');
    ylabel('BG photons');

    disp(' ');
    disp('--- 2nd image: mean values over all measurements ----');
    disp(' ');
    disp(['x_sigma / x_err: ' num2str(mean(x2_sigma)) ' / ' num2str(mean(x2_err)) ' nm']);
    disp(['y_sigma / y_err: ' num2str(mean(y2_sigma)) ' / ' num2str(mean(y2_err)) ' nm']);
    disp(['z_sigma / z_err: ' num2str(mean(z2_sigma)) ' / ' num2str(mean(z2_err)) ' nm']);
    disp(['N_sigma / N_err: ' num2str(mean(N2_sigma)) ' / ' num2str(mean(N2_err)) ]);
    disp(['BG_sigma / BG_err: ' num2str(mean(BG2_sigma)) ' / ' num2str(mean(BG2_err)) ]);
    
    
    %----joint estimates---
    
    figure(5);
    %x-y-estimate
    subplot(1,2,1);
    errorbar(z_vec(z_truth)*1e9,xj_mean,xj_sigma,'rx'); 
    hold on;
    errorbar(z_vec(z_truth)*1e9,yj_mean,yj_sigma,'gx');  
    plot(z_vec(z_truth)*1e9,dx*1e9*ux,'ro');
    plot(z_vec(z_truth)*1e9,dy*1e9*ux,'go');
    grid on;
    hold off;
    title(['joint MLE of x-y-pos.(x=red)']);
    xlabel('z_{truth} / nm');
    ylabel('x/y_{estimate} / nm');

    %z-estimate
    subplot(1,2,2);
    errorbar(z_vec(z_truth)*1e9,zj_mean,zj_sigma,'x'); 
    grid on;
    hold on; 
    plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'-');  
    hold off;
    title(['MLE of z-distance']);
    xlabel('z_{truth} / nm');
    ylabel('z_{estimate} / nm');
    xlim([z_vec(z_truth(1)), z_vec(z_truth(end))]*1e9);
    ylim([0 z_vec(end)]*1e9);
    
    disp(' ');
    disp('--- joint estimation: mean values over all measurements ----');
    disp(' ');
    disp(['x_sigma / x_err: ' num2str(mean(xj_sigma)) ' / ' num2str(mean(xj_err)) ' nm']);
    disp(['y_sigma / y_err: ' num2str(mean(yj_sigma)) ' / ' num2str(mean(yj_err)) ' nm']);
    disp(['z_sigma / z_err: ' num2str(mean(zj_sigma)) ' / ' num2str(mean(zj_err)) ' nm']);

    
end

% evaluating dual-channel data via ratio calculation

if strcmp(ch,'y') %evaluating dual-channel data (intenstiy ratio-calculation --> z-position)
            
        figure(6);
        subplot(2,1,1);
        %Gaussfit-result:
        errorbar(z_vec(z_truth)*1e9,ratio_mean_gauss,ratio_sigma_gauss,'-x'); ylabel('est. ratio from Gaussfits');
        hold on;
        plot(z_vec(z_truth)*1e9,ratio_gt(z_truth),'.-');
        hold off;
        grid on;
        xlim([z_vec(z_truth(1)) z_vec(z_truth(end))]*1e9);
        title(['SAF/UAF-ratio from Gaussfits, photon-no.='  num2str(photon_no)]);
        xlabel('z_{truth} / nm');
        
        %calculating z-position errors from ratio errors
        z_est_gauss=interp1(ratio_gt,z_vec,ratio_mean_gauss);
        z_sigma_plus=interp1(ratio_gt,z_vec,ratio_mean_gauss+ratio_sigma_gauss);
        z_sigma_minus=interp1(ratio_gt,z_vec,ratio_mean_gauss-ratio_sigma_gauss);

        subplot(2,1,2);
        %plot(z_vec(z_truth)*1e9,z_est_gauss*1e9,'x'); axis equal; axis tight;
        errorbar(z_vec(z_truth)*1e9,z_est_gauss*1e9,abs(z_sigma_plus-z_est_gauss)*1e9,abs(z_est_gauss-z_sigma_minus)*1e9,'x')
        grid on;
        hold on;
        plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'-'); 
        hold off;
        xlabel('z_{truth} / nm');
        ylabel('z_{estimate} / nm');
        xlim([z_vec(z_truth(1)), z_vec(z_truth(end))]*1e9);
        ylim([0 z_vec(end)]*1e9);
        title('z-estimates from ratios');

        %---from MLE-------
        
        figure(7);
        subplot(2,1,1);
        plot(z_vec(z_truth)*1e9,ratio_gt(z_truth),'.-');
        hold on;
        %plot MLE-result: 
        errorbar(z_vec(z_truth)*1e9,ratio_mean_MLE,ratio_sigma_MLE,'-x'); ylabel('est. ratio from Gaussfits');
        hold off;
        grid on;
        xlim([z_vec(z_truth(1)) z_vec(z_truth(end))]*1e9);
        title(['SAF/UAF-ratio from MLE, photon-no.='  num2str(photon_no)]);
        xlabel('z_{truth} / nm');
        
        %calculating z-position errors from ratio errors
        z_est_MLE=interp1(ratio_gt,z_vec,ratio_mean_MLE);
        z_sigma_min_MLE=interp1(ratio_gt,z_vec,ratio_mean_MLE+ratio_sigma_MLE);
        z_sigma_plu_MLE=interp1(ratio_gt,z_vec,ratio_mean_MLE-ratio_sigma_MLE);

                
        subplot(2,1,2);
        %plot(z_vec(z_truth)*1e9,z_est_gauss*1e9,'x'); axis equal; axis tight;
        errorbar(z_vec(z_truth)*1e9,z_est_MLE*1e9,(z_sigma_plu_MLE-z_est_MLE)*1e9,(z_est_MLE-z_sigma_min_MLE)*1e9,'x')
        grid on;
        hold on;
        plot(z_vec(z_truth)*1e9,z_vec(z_truth)*1e9,'-'); 
        hold off;
        xlabel('z_{truth} / nm');
        ylabel('z_{estimate} / nm');
        xlim([z_vec(z_truth(1)), z_vec(z_truth(end))]*1e9);
        ylim([0 z_vec(end)]*1e9);
        title('z-estimates from ratios');

    disp(' ');
    disp('--- 2nd image: mean values over all measurements ----');
    disp(' ');
    disp(['x_sigma / x_err: ' num2str(mean(x2_sigma)) ' / ' num2str(mean(x2_err)) ' nm']);
    disp(['y_sigma / y_err: ' num2str(mean(y2_sigma)) ' / ' num2str(mean(y2_err)) ' nm']);
    disp(['N_sigma / N_err: ' num2str(mean(N2_sigma)) ' / ' num2str(mean(N2_err)) ]);
    disp(['BG_sigma / BG_err: ' num2str(mean(BG2_sigma)) ' / ' num2str(mean(BG2_err)) ]);
    disp(' ');
    disp('z-position from MLE-estimated ratios:');
    tmp=(z_sigma_plu_MLE-z_sigma_min_MLE)/2;
    disp(['z_sigma / z_err: ' num2str(mean(tmp(not(isnan(tmp))))*1e9) ' / ' num2str(mean(z2_err)) ' nm']);

end

%%     
