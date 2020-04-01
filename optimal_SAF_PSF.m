%calculating CRB for SAF-detection and finding optimal PSF for SAF-detection by minimizing the CRB:
%-------------------------------------------------------------
%program calculates BFP-fields for various distances (z) of the dipole from the coverslip
%from these data the CRB for x-y-z expectation values are calculated


clear all;
%close all;

global Z PSF uk Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz n_photon bg mask E_tot z_vec

%% calculating BFP fields for SAF

%---user parameters----

N=128; %pupil diameter in pixels

cam=camera('orca fusion');
lens=objective('olympus 1.49');
f_tube=200e-3; %focal lenght of tube lens that was used in the actual experiment

%--defining properties of the PSF to be created
PSF=psf; %initialize PSF class
PSF.lambda=670e-9; %vacuum wavelength
PSF.os=3; %oversampling (set to 3 when calculating stack for performing estimates later)
PSF.uz=5e-9;
PSF.Nx=13*PSF.os; %desired simulated field size in pixel
PSF.Ny=PSF.Nx;
 
T=lens.transmission(N); 
[~,~,R,pupil]=create_coord(N,1,'FFT');


fits='n'; %perform fits to the PSFs? This will slow down the loop; choose 'y' or 'n'
noise='n';  %set to 'y' or 'n'
n_photon=2000; %number of camera counts in the brightest dipole-image
bg=100; %mean background-counts level

z_vec=(0e-9:PSF.uz:250e-9); %simulated dipole distances above layer 2

PSF.Nz=length(z_vec);
RI_buffer=1.33;
%intermediate layer
RI_layer=RI_buffer; %refractive index of an intermediate layer on top of the coverlip (e.g. lipid bilayer)
d_layer=0e-9; %thickness of intermediate layer (layer 2)

%dz_vec=(-0.9:.005:0.1)*1e-6; %vector of defocus values - to simulate 3D PSFs
dz_vec=-500e-9; %only a single defocus 

%----------------------------------------------------------------------

mu_x=1e-10; %magnitude of x-dipole (arbitrary value)
mu_y=1e-10; %magnitude of y-dipole (arbitrary value)
mu_z=1e-10; %magnitude of z-dipole (arbitrary value)
mu=[mu_x mu_y mu_z]; %vector of dipole strengths

PSF.ux=cam.pixsize/lens.M*lens.f_tube/f_tube/PSF.os; %effective pixel size in focal plane
RI=[RI_buffer RI_layer lens.RI]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]

%----loading exp. aberrations-----
%Z_aberr2=zeros(1,35);
%load('coeff_2020-02-07_NA1,47_err0,25.mat'); Z_aberr2=Z_phase2;
%PSF.Zernike=[modes; Z_phase2];
if exist('Z_amp')
    %additional apodization in the objective pupil
    Apo_mode(:,:,1)=1.*pupil;
    Apo_mode(:,:,2)=-pupil.*(R/R(end/2,end)).^2+1;
    Apo_mode(:,:,3)=-pupil.*(R/R(end/2,end)).^4+1;
    Apo_mode(:,:,4)=-pupil.*(R/R(end/2,end)).^6+1;
    
    tmp=zeros(1,1,size(Apo_mode,3)); 
    tmp(1,1,:)=Z_amp;
    Apo=sum(Apo_mode.*repmat(tmp,[N,N,1]),3);
    Apo=Apo/max(Apo(:));
else
    Apo=1;
end

%Z_aberr2=zeros(1,length(Z_aberr2)); Z_aberr2(3)=0; Z_aberr2(4)=0;  
%Z_aberr2=Z_aberr2;
%for aberrations; %deactivate line above to simulate aberration-free case
%---------------------------------

for mm=1:length(dz_vec)
   dz=dz_vec(mm);
   
%dz=-300e-9;  %defocus of objective lens (negative values mean moving the focus into the fluid)

    [~,Defocus,~] =fun_SA_RImismatch(N,RI(3),RI(3),lens.NA,PSF.lambda,1); %Defocus function refers to refractive index n2
    [SA_out,~,~] =fun_SA_RImismatch(N,RI(3),RI(1),lens.NA,PSF.lambda,1); %Defocus function refers to refractive index n2

    uk=4*pi/PSF.lambda*lens.NA/N; %unit in pupil space (k-space)
    

    pupil_UAF=circshift(R<=((N/2)*(RI(1)/lens.NA))*1,[0 0]); %pupil containing UAF light

    %-----considering aberrations-----
    if exist('Z_aberr2') %if aberration Zernikes are loaded
        Z_aberr2(1:2)=0; %delete tip, tilt to center the PSF
        bar(Z_aberr2); 
        aberr=sum(ZernikeCalc([2:(1+length(Z_aberr2))],Z_aberr2',pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded 
    else
        aberr=0; 
        %aberr=sum(ZernikeCalc([4 11 22],[2.1603    1.1705    1.1048]',pupil,'NOLL'),3);
        disp('assuming NO aberrations');
    end

%calculating BFP-fields for all dipole orientations
for m=1:length(z_vec)
    dipole=[0,0]; %[theta, phi], e.g. [0,0] for z-dipole, [pi/2,0] for x-dipole
    [Ex_Pz(:,:,m),Ey_Pz(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,RI,[0,0],d_layer,z_vec(m),lens.f,mu_z,T); %z-dipole
    [Ex_Px(:,:,m),Ey_Px(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,RI,[pi/2, 0],d_layer,z_vec(m),lens.f,mu_x,T); %x-dipole
    [Ex_Py(:,:,m),Ey_Py(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,RI,[pi/2, pi/2],d_layer,z_vec(m),lens.f,mu_y,T); %y-dipole
end
disp('done');

% display phase in BFP
%imagesc(mod(angle(Ex_Px(:,:,m)),2*pi)); colormap hsv; colorbar; axis equal; axis tight; title('angle(Ey); Dipole Py is emitter')

%% calculating BFP images

clear I_BFP ratio I_BFP

for m=1:length(z_vec) %BFP image number

        %user-defined additional pupil mask (incorporates objective transmission profile):
        %phase=real(-SA_out).*(-dz_vec).*pupil_UAF; %correcting the RI-mismatch 
        phase=0; 
        mask=Apo.*pupil.*exp(1i*phase+1i*aberr+1i*dz*Defocus);

        figure(1);
        %I_BFP(m,:,:)=(abs(Ex_Px(:,:,m).*mask).^2+abs(Ex_Py(:,:,m).*mask).^2+abs(Ey_Px(:,:,m).*mask).^2+abs(Ey_Py(:,:,m).*mask).^2+abs(Ex_Pz(:,:,m).*mask).^2+abs(Ey_Pz(:,:,m).*mask).^2);
        I_BFP(m,:,:)=abs(mask).^2.*(abs(Ex_Px(:,:,m)).^2+abs(Ex_Py(:,:,m)).^2+abs(Ey_Px(:,:,m)).^2+abs(Ey_Py(:,:,m)).^2+abs(Ex_Pz(:,:,m)).^2+abs(Ey_Pz(:,:,m)).^2);

        E_tot(m)=sum(sum(I_BFP(m,:,:)));
        E_UAF(m)=sum(sum(squeeze(I_BFP(m,:,:)).*pupil_UAF));
        
        imagesc(squeeze(I_BFP(m,:,:))); %axis equal; axis tight; colorbar; colormap gray;
        axis equal; axis tight; colormap gray;
        title(['BFP intens. for z=' num2str(z_vec(m))]);
        pause(0);
        %plot(I_BFP(end/2,:));
        %imwrite(squeeze(I_BFP(m,:,:))/max(I_BFP(:)),['IBFP_' num2str(m) '_z=' num2str(z_vec(m)) '_' num2str([RI(1) RI(3)]) '.tiff']);
end

ratio=(E_tot-E_UAF)./E_UAF; %SAF/UAF ratio

figure(2);
plot(z_vec*1e9,ratio); xlabel('z / nm'); ylabel('SAF/UAF ratio'); grid on;
title(['NA=' num2str(lens.NA) ', RIs=' num2str(RI) ', \lambda_0=' num2str(PSF.lambda*1e9)]);
hold on;


%% -----calculating PSF as seen on the camera for different emitter z-positions-----
% calc of CCD images and CRBs for all z-values contained in z_vec

clear PSF_tot PSF_SAF PSF_UAF Gfit_UAF Gfit_SAF

if strcmp(fits,'y')
%preparations for Gaussfit:
    x_vec=linspace(-PSF.Nx/2+1,PSF.Nx/2,PSF.Nx); %required for Gaussfits
    x_data=zeros(length(x_vec), length(x_vec),2); %for lsqcurvefit Gaussfit
    [tmpX, tmpY]=meshgrid(x_vec,x_vec);
    x_data(:,:,1)=tmpX;
    x_data(:,:,2)=tmpY;
end

figure(1);
boundary=zeros(PSF.Nx,PSF.Nx); %border-mask to estimate background from PSF-images
tmp=1/(4*(PSF.Nx-1)); %normalization such that sum(boundary(:)) equals 1
boundary(1,:)=tmp; boundary(end,:)=tmp; boundary(:,1)=tmp; boundary(:,end)=tmp; 


for m=1:length(z_vec)
    %I_CCD_x=abs(fftshift(fft2(ifftshift(embed(pupil.*(E_x(:,:,m)),N_pad,0))))).^2;
        
    %-----calculating total (SAF+UAF) images-----
    
    I_xx=abs(czt2(Ex_Px(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    PSF_tot(:,:,m)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m); %normalization to total energy in BFP; we thus assume that the SAME energy comes through the objective lens for every z-distance (note that the energy would otherwise be different in the simulation for varying z-distances)
      
    %if m==1; C_norm=sum(sum(PSF_tot(:,:,1))); end %normalization to total intensity in first image (m=1)
    %C_norm=sum(sum(PSF_tot(:,:,m)));  %normalization if info is contained in shape-changes
    
    if strcmp(noise,'y') %if noise is selected
        tmp=PSF_tot(:,:,m)*n_photon+bg; 
        PSF_tot(:,:,m)=poissrnd(tmp,size(tmp,1),size(tmp,2)); %normalization of PSF   
    end
    
    est_offset=sum(sum(boundary.*PSF_tot(:,:,m))); %returns the mean value in the boundary-region
    energy_tot(m)=sum(sum(PSF_tot(:,:,m)-est_offset));
    
    %-----calculating UAF-images-----
    I_xx=abs(czt2(Ex_Px(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:,m).*mask.*pupil_UAF,uk,PSF.ux,PSF.Nx)).^2;
    PSF_UAF(:,:,m)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m);    
    
    if strcmp(noise,'y') %if noise is selected
        tmp=PSF_UAF(:,:,m)*n_photon+bg; 
        PSF_UAF(:,:,m)=poissrnd(tmp,size(tmp,1),size(tmp,2)); %normalization of PSF   
    end
    est_offset=sum(sum(boundary.*PSF_tot(:,:,m))); %returns the mean value in the boundary-region
    energy_UAF(m)=sum(sum(PSF_UAF(:,:,m)-est_offset));

    
    if strcmp(fits,'y')

    %----- performing Gaussfits/double-Gaussfits to PSF_UAF and PSF_tot-----
        
        %----fit to UAF-PSFs----:
        %offset=3;
        tmp=squeeze(PSF_UAF(:,:,m));%+offset;

        %FIT METHOD A: one 2D-Gauss with offset: [offset amplitude x-shift y-shift width]
        param_ini=[0 max(tmp(:)) 1 1 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
        lb=[0 0 -length(x_vec)/2 -length(x_vec)/2 0.5]; %lower bounds for fit-parameters
        ub=[max(tmp(:)) max(tmp(:)) length(x_vec)/2 length(x_vec)/2 length(x_vec)/2]; %upper bounds
        [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_gauss_and_offset,param_ini,x_data,tmp,lb,ub);
        Gaussfit=fun_gauss_and_offset(Param,x_data);
        Gfit_UAF(m)=Param(2)*2*pi*Param(5)^2; %energy contained (see e.g. formula in Thunderstorm script)

        %FIT METHOD B: two Gaussians with offset: [offset amp1 x-shift y-shift width1 amp2 width2]
%         param_ini=[50 max(tmp(:)) 1 1 1 max(tmp(:))/8 3]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
%         lb=[0 0 -length(x_vec)/2 -length(x_vec) 0 0 0]; %lower bounds for fit-parameters
%         ub=[max(tmp(:)) max(tmp(:)) length(x_vec)/2 length(x_vec) length(x_vec) max(tmp(:))/3 6]; %upper bounds
%         [Param,resnorm_UAF,residual,exitflag]=lsqcurvefit(@fun_2gauss_and_offset,param_ini,x_data,tmp,lb,ub);
%         Param_UAF(m,:)=Param;
%         Gaussfit=fun_2gauss_and_offset(Param,x_data);
%         Gfit_UAF(m)=(Param(2)*Param(5)^2+Param(6)*Param(7)^2)*2*pi; %energy contained (see e.g. formula in Thunderstorm script)

        subplot(2,1,1); 
        plot(x_vec,sum(tmp,1),x_vec,sum(Gaussfit,1)); 
        ylabel('UAF');

    %----fit to TOT-PSFs----:
        tmp=squeeze(PSF_tot(:,:,m));%+offset;
    
        %FIT METHOD A: single Gauss + offset
        param_ini=[0 max(tmp(:)) 1 1 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
        lb=[0 0 -length(x_vec)/2 -length(x_vec)/2 0.5]; %lower bounds for fit-parameters
        ub=[max(tmp(:)) max(tmp(:)) length(x_vec)/2 length(x_vec)/2 length(x_vec)/2]; %upper bounds
        [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_gauss_and_offset,param_ini,x_data,tmp,lb,ub);
        Gaussfit=fun_gauss_and_offset(Param,x_data);
        Gfit_tot(m)=Param(2)*2*pi*Param(5)^2; %energy contained (see e.g. formula in Thunderstorm script)
        sigma(m)=Param(5);
        
        
        % 
        %FIT METHOD B: two Gaussians with offset: [offset amp1 x-shift y-shift width1 amp2 width2]
%         param_ini=[50 max(tmp(:)) 1 1 1 max(tmp(:))/5 3]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
%         lb=[0 0 -length(x_vec)/2 -length(x_vec) 0 0 0]; %lower bounds for fit-parameters
%         ub=[max(tmp(:)) max(tmp(:)) length(x_vec)/2 length(x_vec) length(x_vec) max(tmp(:))/3 6]; %upper bounds
%         [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_2gauss_and_offset,param_ini,x_data,tmp,lb,ub);
%         Param_tot(m,:)=Param;
%         Gaussfit=fun_2gauss_and_offset(Param,x_data);
%         Gfit_tot(m)=(Param(2)*Param(5)^2+Param(6)*Param(7)^2)*2*pi; %energy contained (see e.g. formula in Thunderstorm script)

        %imagesc(Gaussfit); pause(0);
        subplot(2,1,2); 
        plot(x_vec,sum(tmp,1),x_vec,sum(Gaussfit,1)); pause(0.1);
        ylabel('TOT');
    end
end

%calculating and plotting SAF/UAF ratio from PSFs (=bead images)
figure(2);
plot(z_vec*1e9,(energy_tot-energy_UAF)./energy_UAF,'r.-'); grid on; %ratios determined from integration over PSFs
hold on; 
plot(z_vec*1e9,ratio,'g.-'); grid on; %ratio determined directly in BFP
title('ratios obtained from BFP (green) and simple integration (red)');
if strcmp(fits,'y')
    plot(z_vec*1e9,(Gfit_tot-Gfit_UAF)./Gfit_UAF,'b.-'); %ratios determined from Gaussfits
end

hold off;
ylabel('SAF/UAF-ratio'); title('SAF/UAF-ratios: red=integr. over PSFs; green=measured in BFP; blue=from Gaussfits');
xlabel('z / nm');
     
    %plotting total PSFs for different z-values
    z_plot=round([1, 2, 2*length(z_vec)/3, length(z_vec)]);
    figure(3);
    subplot(2,2,1);
    imagesc(PSF_tot(:,:,z_plot(1))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(z_vec(z_plot(1)))]);
    subplot(2,2,2);
    if length(z_vec)>=4
        imagesc(PSF_tot(:,:,z_plot(2))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(z_vec(z_plot(2)))]);
        subplot(2,2,3);
        imagesc(PSF_tot(:,:,z_plot(3))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(z_vec(z_plot(3)))]);
        subplot(2,2,4);
        imagesc(PSF_tot(:,:,z_plot(4))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(z_vec(z_plot(4)))]);
    end
    
%%%saving tot- and UAF-PSF images to tif-stacks
%     for m=1:length(z_vec);
%         imwrite((PSF_tot(:,:,m)/max(PSF_tot(:))),['PSF_tot_'  ' PSF.uz=' num2str(PSF.uz) '_focus=' num2str(dz) '.tif'],'WriteMode','append');
%         imwrite((PSF_UAF(:,:,m)/max(PSF_tot(:))),['PSF_UAF_'  ' PSF.uz=' num2str(PSF.uz) '_focus=' num2str(dz) '.tif'],'WriteMode','append');
%     end
if length(dz_vec)>1 
    PSF_defocus(:,:,mm)=PSF_tot; %if mm-loop is activated; %here you can choose between PSF_tot or PSF_UAF
end

end  %end of mm-index loop 


if exist('Z_aberr2')~=1
    Z_aberr2=0; 
end
if exist('Z_amp')~=1
    Z_amp=0; 
end

PSFpath='c:/users/q004aj/desktop/PSFs/NA1,49/';
PSFname=[PSFpath 'PSF_NA' num2str(lens.NA) '_' num2str(PSF.Nx/PSF.os) 'x' num2str(PSF.Nx/PSF.os) '_' num2str(z_vec(1)*1e9) '-' num2str(PSF.uz*1e9) '-' num2str(z_vec(end)*1e9) 'nm_RI=' num2str(RI(1),3) '_dz=' num2str(dz*1e9) '_aberrfree.mat'];
PSF.data=PSF_tot;

%save(PSFname,'PSF','lens','RI','cam');

%% optional: preparing 5D-PSF-model for use with "MLE_fit_molecules_exp" or "MLE_fit_molecules_2channels_exp"
% PSF-model is expanded to 5D (interpolated along x-y directions)
% this avoids the use of interpolation in the log-likelihood function ->
% faster

if length(dz_vec)>1 %if multiple defocus-values are defined
    Nz=length(dz_vec);
    PSF_tmp=PSF_defocus;
    disp('creating 5D "defocus" PSF');
elseif length(z_vec)>1   %if multiple z-distances are defined (and only one defocus-value!)
    Nz=length(z_vec);
    PSF_tmp=PSF_tot;%here you can choose between PSF_tot or PSF_UAF
    disp('creating 5D SAF-PSF');
end

Ns=101; %interpolation-steps in x-y
w=1; %(half) range of molecule x-y-offsets (in camera pixels) that should be covered with the interpolation
sx=linspace(-w*PSF.os,w*PSF.os,Ns);
sy=sx; 
interp_incr=(sx(2)-sx(1))/PSF.os; %interpolation increment in units of pixels

PSF5D=zeros(round(PSF.Nx/PSF.os)-2*w,round(PSF.Nx/PSF.os)-2*w,Nz,Ns,Ns);
for mz=1:Nz
    for mx=1:Ns
        for my=1:Ns 
            tmp=PSF.os^2*imresize(interp2(PSF_tmp(:,:,mz),(1:PSF.Nx)-sx(mx),(1:PSF.Nx)'-sy(my)),1/PSF.os,'box');
            PSF5D(:,:,mz,mx,my)=tmp(1+w:end-w,1+w:end-w); %cut the outer pixel frame because of potential sampling artifacts
        end
    end
    disp(Nz-mz);
end
disp('done');
PSF5D(isnan(PSF5D))=0;

PSFpath='c:/users/q004aj/desktop/PSFs/NA1,50/';
name=[PSFpath 'PSF5D_NA' num2str(lens.NA) '_' num2str(size(PSF5D,1)) 'x' num2str(size(PSF5D,1)) '_' num2str(z_vec(1)*1e9) '-' num2str(PSF.uz*1e9) '-' num2str(z_vec(end)*1e9) 'nm_RI=' num2str(RI(1),3) '_dz=' num2str(dz*1e9) '_aberrfree.mat'];



%save(name,'PSF5D','z_vec','PSF.ux','lens','RI','interp_incr','PSF.os','Z_aberr2','Z_amp','mu'); 
    
%% ---calculating CRBs----
%for calculation of CRBs, the PSFs are individually normalized to contain the same
%energy -- this "deletes" information contained in differences of total energy
%between individual PSFs, which are also not available in experiment. 

%optional: investigating a combination of two molecule images (two-channel imaging)
    %load('PSF_0-3-250nm_RI=1,45_defoc=-500nm_aberrfree.mat'); PSF_A=PSF_UAF; PSF_B=PSF_UAF;
    %load('PSF_0-3-250nm_RI=1,45_defoc=0nm_aberrfree.mat'); PSF_B=PSF_tot;
    %PSF_tot=[PSF_A PSF_B];
    %PSF_tot=PSF_B; %single channel 
    %PSF_tot=PSF_defocus; %if a "defocus"-stack has been calculated
    
    
nor=0; %should each z-slice of the PSF be normalized individually to contain n_photon photons?
[CRBx,CRBy,CRBz,~,~,FI]=fun_CRB(PSF,n_photon,bg,cam.gain,nor); %note that if PSF is oversampled the background needs to be adapted correspondingly
figure(4);
plot((1:length(CRBz))*PSF.uz*1e9,sqrt(CRBz),'b.-'); xlabel('z-pos in nm'); ylabel('nm');
title(['sqrt(CRBz), cts=' num2str(n_photon) ', bg=' num2str(bg)]); grid on;
%ylim([0 70]);

figure(5);
plot((1:length(CRBz))*PSF.uz*1e9,sqrt(CRBx)); xlabel('z-pos in nm'); ylabel('nm');
hold on; 
plot((1:length(CRBz))*PSF.uz*1e9,sqrt(CRBy),'r.-'); xlabel('z-pos in nm'); ylabel('nm');
title(['sqrt(CRBx), sqrt(CRBy) (red), cts='  num2str(n_photon) ', bg=' num2str(bg)]);
hold off;
grid on;
%ylim([0 15]);

metric1=mean(sqrt((CRBx.*CRBy.*CRBz)).^(1/3))  %"localization volume"
metric2=mean(sqrt(CRBz))

%plotting FIsher infomration 

figure(6); 
plot((1:length(CRBz))*PSF.uz*1e9,squeeze(FI(1,1,:)));
hold on; 
plot((1:length(CRBz))*PSF.uz*1e9,squeeze(FI(3,3,:)));
title('Fisher information x and z');
xlabel('z-pos in nm'); ylabel('1/nm²');
hold off; 

%% A) for ZENIKE OPTIMIZATION: calculating Zernikes; the Zernike table is stored in Z and used by the
%sub-routine "fun_optPSFforSAF.m"
tmp=[4,6]; %Zernike modes (Noll scheme) to be included in the optimization
Z=ZernikeCalc(tmp,ones(length(tmp),1),pupil,'Noll');
no_modes=size(Z,3); 

a=rand(1,size(Z,3));
a2=zeros(1,1,length(a));
a2(1,1,:)=a;
phase=sum(repmat(a2,[size(Z,1),size(Z,2),1]).*Z,3); %pupil phase
imagesc(phase); title('initial pupil phase'); pause(0.01);
a_ini=zeros(no_modes,1); %initial vector of Zernikes, beginning with astigmatism (no. 5)

%% B) for stepped phase mask optimization
%clear Z
a=[1 0]; %phase-step between UAF and SAF
no_modes=2;
Z(:,:,1)=pupil_UAF;
%Z(:,:,1)=a(1)*zernike_modes(N,N,[ones(1,3) 1]);
%Z(:,:,2)=a(1)*zernike_modes(N,N,[ones(1,10) 1]);
%Z(:,:,3)=a(1)*zernike_modes(N,N,[ones(1,21) 1]);
%Z=ZernikeCalc([4,11,22],[1 1 1]',pupil,'NOLL');

Z(:,:,1)=Defocus*1e-6; 
%Z(:,:,2)=[];
a_ini=0;

%solution: Z=[2.1603    1.1705    1.1048]

%% minimizing CRB using fmincon

%a_opt=fmincon(@fun_optPSFforSAF,a,eye(no_modes),ones(no_modes,1))
options=optimset('Display','iter','Tolfun',1e-5);
a_opt=fminsearch(@fun_optPSFforSAF,a_ini,options);
disp(a_opt);

%% minimizing CRB for biplane imaging

Z=Defocus*1e-6;
split_ratio_ini=0.5; 
a_ini=[-0.25 -0.6 split_ratio_ini]; %initial defocus values in micron
options=optimset('Display','iter','Tolfun',1e-5);

a_opt=fminsearch(@fun_optPSFforSAF_biplane,a_ini,options);
disp(a_opt);
