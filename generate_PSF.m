%calculation of 3D PSFs

clear all;
close all;

%% ---user parameters----

N=128; %pupil diameter in pixels
z_max=250e-9; %maximal z-value (working range)

cam=camera('orca fusion'); %e.g. 'orca fusion'
lens=objective('olympus 1.49');

f_tube=200e-3; %focal length of tube lens that was used in the actual experiment

%--defining properties of the PSF to be created
PSF=psf; %initialize PSF class
PSF.lambda=670e-9; %vacuum wavelength

%z_max=PSF.lambda/2; 

PSF.os=2; %oversampling (set to 3 when calculating stack for performing estimates later)
PSF.uz=20e-9; %z-increment in meter
PSF.Nx=15*PSF.os; %desired simulated field size in pixel
PSF.Ny=PSF.Nx;
z_vec=(0e-9:PSF.uz:z_max); %z-range, i.e., simulated dipole distances above layer 2 
PSF.Nz=length(z_vec);
PSF.RI=1.33; %refractive index of buffer
PSF.defocus=-70e-9; %negative values = moving objective towards sample

T=lens.transmission(N); 
[X,Y,R,pupil]=create_coord(N,1,'FFT');

UseAberr='y'; %use aberrations?
    %if yes, choose Zernike vector (1st row: modes; 2nd row: magnitudes)and
    %apodization pupil image or load aberration file;
    if strcmp(UseAberr,'y')
        %PSF.Zernike=[5 ; -1];Apo=1; %define modes and magnitudes here
        aberr_filename='aberr_coeff_cyl.mat'; 
    end
    
%define signal and BG for CRLB calculations
sig=2000; %number of camera counts in the brightest dipole-image
bg=20; %mean background level in photons per pixel

%intermediate layer
RI_layer=PSF.RI; %refractive index of an intermediate layer on top of the coverlip (e.g. lipid bilayer)
d_layer=0e-9; %thickness of intermediate layer (layer 2)

%----------------------------------------------------------------------

mu_x=1e-10; %magnitude of x-dipole (arbitrary value)
mu_y=1e-10; %magnitude of y-dipole (arbitrary value)
mu_z=1e-10; %magnitude of z-dipole (arbitrary value)
mu=[mu_x mu_y mu_z]; %vector of dipole strengths

PSF.ux=cam.pixsize/lens.M*lens.f_tube/f_tube/PSF.os; %effective pixel size in focal plane
RI=[PSF.RI RI_layer lens.RI]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]

%----loading experimentally determined aberrations-----
if strcmp(UseAberr,'y')
    
    if isempty(PSF.Zernike)
        load(aberr_filename); 
        Z_aberr=Z_phase2;
        PSF.Zernike=[modes; Z_aberr];
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
    end
else
            Apo=1;
end

dz_vec=PSF.defocus;

for mm=1:length(dz_vec) %loop if a vector of different defocus values has been defined
    dz=dz_vec(mm);

    [~,Defocus,~] =fun_SA_RImismatch(N,lens.RI,lens.RI,lens.NA,PSF.lambda,1); %Defocus function refers to refractive index n2
    [SA_out,~,~] =fun_SA_RImismatch(N,lens.RI,PSF.RI,lens.NA,PSF.lambda,1); %Defocus function refers to refractive index n2

    uk=4*pi/PSF.lambda*lens.NA/N; %unit in pupil space (k-space)
    pupil_UAF=circshift(R<=((N/2)*(PSF.RI/lens.NA))*1,[0 0]); %pupil containing UAF light

    %-----considering aberrations-----
    if strcmp(UseAberr,'y') %if aberration Zernikes are loaded
        
        %delete tip, tilt (if present) to center the PSF
        idx1=find(PSF.Zernike(1,:)==2); 
        idx2=find(PSF.Zernike(1,:)==3); 
        PSF.Zernike(2,idx1)=0;
        PSF.Zernike(2,idx2)=0;
         
        bar(PSF.Zernike(1,:),PSF.Zernike(2,:));
        xlabel('Zernike modes (Noll scheme)');
        ylabel('rad RMS');
        aberr=sum(ZernikeCalc(PSF.Zernike(1,:),PSF.Zernike(2,:)',pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded 
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

%% calculating BFP images

clear I_BFP ratio I_BFP

for m=1:length(z_vec) %BFP image number

        %user-defined additional pupil mask (incorporates additional objective transmission profile as defined in the aberrations file):
        %phase=real(-SA_out).*(-dz_vec).*pupil_UAF; %correcting the RI-mismatch 
        phase=0; 
        
        %for simulation of a UAF-PSF, replace "pupil" by "pupil_UAF" in the
        %line below. 
        mask=Apo.*pupil.*exp(1i*phase+1i*aberr+1i*dz*Defocus);
       
        I_BFP(m,:,:)=abs(Apo.*pupil).^2.*(abs(Ex_Px(:,:,m)).^2+abs(Ex_Py(:,:,m)).^2+abs(Ey_Px(:,:,m)).^2+abs(Ey_Py(:,:,m)).^2+abs(Ex_Pz(:,:,m)).^2+abs(Ey_Pz(:,:,m)).^2);
        E_tot(m)=sum(sum(I_BFP(m,:,:)));
        E_UAF(m)=sum(sum(squeeze(I_BFP(m,:,:)).*pupil_UAF));
        
        %show BFP intensity
        figure(1);
        imagesc(squeeze(I_BFP(m,:,:))); %axis equal; axis tight; colorbar; colormap gray;
        axis equal; axis tight; colormap gray;
        title(['BFP intens. for z=' num2str(z_vec(m)*1e9) ' nm']);
        pause(0);
end

ratio=(E_tot-E_UAF)./E_UAF; %SAF/UAF ratio

figure(2);
plot(z_vec*1e9,ratio); xlabel('z / nm'); ylabel('SAF/UAF ratio'); grid on;
title(['NA=' num2str(lens.NA) ', RIs=' num2str(RI) ', \lambda_0=' num2str(PSF.lambda*1e9)]);
hold on;

%% -----calculating PSF as seen on the camera for different emitter z-positions-----
% calc of CCD images and CRBs for all z-values contained in z_vec

clear PSF_tot PSF_SAF PSF_UAF Gfit_UAF Gfit_SAF

for m=1:length(z_vec)
        
    I_xx=abs(czt2(Ex_Px(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:,m).*mask,uk,PSF.ux,PSF.Nx)).^2;
    PSF_tot(:,:,m)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m); %normalization to total energy in BFP; we thus assume that the SAME energy comes through the objective lens for every z-distance (note that the energy would otherwise be different in the simulation for varying z-distances)
    
end
     
    %plotting total PSFs for different z-values
    z_plot=round([1, 2, 2*length(z_vec)/3, length(z_vec)]);
    figure(3);
    subplot(2,2,1);
    imagesc(PSF_tot(:,:,z_plot(1))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(z_vec(z_plot(1)))]);
    subplot(2,2,2);
    if length(z_vec)>=4
        imagesc(PSF_tot(:,:,z_plot(2))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(1e9*z_vec(z_plot(2))) ' nm']);
        subplot(2,2,3);
        imagesc(PSF_tot(:,:,z_plot(3))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(1e9*z_vec(z_plot(3))) ' nm']);
        subplot(2,2,4);
        imagesc(PSF_tot(:,:,z_plot(4))); colorbar; axis equal; axis tight; title(['tot-PSF, z=' num2str(1e9*z_vec(z_plot(4))) ' nm']);
    end
    
if length(dz_vec)>1 
    PSF_defocus(:,:,mm)=PSF_tot; %if mm-loop is activated; %here you can choose between PSF_tot or PSF_UAF
end

end  %end of mm-index loop 


if exist('Z_aberr2')~=1
    Z_aberr=0; 
end
if exist('Z_amp')~=1
    Z_amp=0; 
end

%PSFpath='c:/users/q004aj/desktop/PSFs/NA1,49/';
PSFpath='';
PSFname=[PSFpath 'PSF_Cyl_NA' num2str(lens.NA) '_' num2str(PSF.Nx/PSF.os) 'x' num2str(PSF.Nx/PSF.os) '_' num2str(z_vec(1)*1e9) '-' num2str(PSF.uz*1e9) '-' num2str(z_vec(end)*1e9) 'nm_RI=' num2str(RI(1),3) '_dz=' num2str(dz*1e9) 'nm_aberr=' UseAberr '.mat'];
PSF.data=PSF_tot;
%save(PSFname,'PSF','lens','RI','cam');
  
%--------calculating CRBs------------------------------

%nor=0; %should each z-slice of the PSF be normalized individually to contain n_photon photons?
eta=1; %detection efficiency of channel
[CRBx,CRBy,CRBz,CRBsig,~,FI]=PSF.CRLB(sig,bg,cam,eta);

figure(4);
plot(z_vec*1e9,sqrt(CRBz),'r'); xlabel('z-pos in nm'); ylabel('nm');
title(['sqrt(CRBz), cts=' num2str(sig) ', bg=' num2str(bg)]); grid on;
hold on; 
plot(z_vec*1e9,sqrt(CRBx),'b'); xlabel('z-pos in nm'); ylabel('nm');
plot(z_vec*1e9,sqrt(CRBy),'g'); xlabel('z-pos in nm'); ylabel('nm');
title(['x,y,z precisions(z=red), cts='  num2str(sig) ', bg=' num2str(bg)]);
hold off;
grid on;

%plotting Fisher information 
% figure(6); 
% plot((1:length(CRBz))*PSF.uz*1e9,squeeze(FI(1,1,:)));
% hold on; 
% plot((1:length(CRBz))*PSF.uz*1e9,squeeze(FI(3,3,:)));
% title('Fisher information x and z');
% xlabel('z-pos in nm'); ylabel('1/nm²');
% hold off; 


%% Optimizing PSF

clear Z;

PSFtype=1; %1...Defocus
           %2...Zernikes (z.B. Astigm.) + Defocus
           %3...Biplane 
           %4...phase ramp
           %5...Donald

%----for optimal defocus in single-channel imaging------
if PSFtype==1

    Z=Defocus*1e-6; 
    a_ini=-0.5;
    method='single ch., defocus';

elseif PSFtype==2

    %----allowing for many zernike modes
    modes=[6]; 
    Z(:,:,1)=Defocus*1e-6;
    Z(:,:,2:length(modes)+1)=ZernikeCalc(modes,ones(length(modes),1),pupil,'NOLL');
    a_ini=ones(1,size(Z,3));
    method=['single ch.: ' num2str(modes)]; 

elseif PSFtype==3

    %----for bi-plane imaging
    method='biplane'; 
    Z(:,:,1)=Defocus*1e-6; 
    Z(:,:,2)=ZernikeCalc(6,1,pupil,'NOLL');%additional Zernike term for channel 1
    Z(:,:,3)=ZernikeCalc(6,1,pupil,'NOLL');%additional Zernike term for channel 2
    a_ini=[0 -0.5]; %[defocus x, defocus y, split-ratio, Z1, Z2]

elseif PSFtype==4 
    
    %phase ramp (Baddeley et al, Nano Res. 2011)
    Z(:,:,1)=Defocus*1e-6; 
    tmp=ZernikeCalc(3,1,pupil,'NOLL');
    Z(:,:,2)=(X<0).*tmp+(X>=0).*(-tmp); %phase ramp
    a_ini=[0 0]; 
    method=['phase ramp']; 

elseif PSFtype==5

    %DONALD
    method='donald'; 
    Z(:,:,1)=Defocus*1e-6; 
    %Z(:,:,2)=Z(:,:,1); 
    a_ini=[-0.5]; %one defocus value for both channels
       
end
%----------------------

E_tensor=zeros(N,N,PSF.Nz,6); 
E_tensor(:,:,:,1)=Ex_Px; 
E_tensor(:,:,:,2)=Ex_Py; 
E_tensor(:,:,:,3)=Ex_Pz; 
E_tensor(:,:,:,4)=Ey_Px; 
E_tensor(:,:,:,5)=Ey_Py; 
E_tensor(:,:,:,6)=Ey_Pz; 

pupilmask=exp(1i*aberr); %additional, fixed, complex pupil mask
metric=2; %1...x-y-z; 2...optimal z
a_opt=fun_optimize_PSF(PSF,lens,cam,E_tensor,pupilmask,Z,a_ini,sig,bg,metric,method,pupil_UAF);
disp('optimal parameters: ');
disp(a_opt);

grid on; ylim([0 50]); 


%% for paper - generating CRLB-maps for defined parameter ranges

if PSFtype==1 %defocus
    %a_opt=-0.54; 
    a{1}=a_opt+[-1:0.1:1]; a{2}=1;  
elseif PSFtype==2  %ast
    a{1}=a_opt(1)+[-0.6:0.1:0.6]; %defocus
    a{2}=a_opt(2)+[-1:0.1:1]; %Z6
elseif PSFtype==3  %biplane
    a{1}=a_opt(1)+[-0.5:0.05:0.5]; %defocus 1
    a{2}=a_opt(2)+[-0.5:0.05:0.5]; %defocus 2
elseif PSFtype==5 %Donald
    a{1}=a_opt(1)+[-0.5:0.05:0.5]; %defocus 1   
    %a{2}=a_opt(1)+[-0.5:0.05:0.5]; %defocus 1    
end

map=fun_parameter_maps(PSF,lens,cam,E_tensor,pupilmask,Z,a,sig,bg,method);

%rearrange dimensions if only one parameter was varied: 
if ndims(map)==2 
    tmp=map; 
    map=ones(3,1,length(tmp)); 
    map(:,1,:)=tmp(:,:); 
end    

%% showing parameters maps

figure(1);

subplot(3,1,1);
imagesc(a{1},a{2},(squeeze(map(1,:,:)))'); 
colormap jet; 
colorbar; 
axis equal; axis tight; 

subplot(3,1,2);
imagesc(a{1},a{2},(squeeze(map(2,:,:)))'); 
colorbar; 
axis equal; axis tight; 

subplot(3,1,3);
imagesc(a{1},a{2},(squeeze(map(3,:,:)))'); 
axis equal; axis tight; 
colorbar; 

%% saving stability map

%save([cd '\stab maps\stab_map_salm.mat'],'map','a')