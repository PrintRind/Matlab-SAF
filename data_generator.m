%Data generator for "MLE_fit_molecules.m" and
%"MLE_fit_molecules_2channels.m"

%simulated "raw data" (i.e. molecule images) are generated from one (or two) PSF-model(s) 

clear all; 
close all; 

%%
no_images=200; %number of simulated CCD images (=length of stack)
n_img=400;  %pixel size length of each simulated camera image
rho=0.001; %average molecule density 1/µm^2
ux=115e-9;  %pixel size in m
sig=3000; %signal per molecule in photons (refers to signal of brightest molecule)
BG_top=500; %background level in photons in top image
BG_bottom=500; %background level in photons in bottom image

%camera parameters
gain=100; 
amp=9.9; %electrons per count 
QE=0.95; %quantum efficiency

r_sphere=2.0e-3/2; %radius of calibration sphere (Alexa-coated ball lens)

%distortion of second image (NOTE: bottom image is also flipped-up/down in
%addition to the settings below)
theta=5;  %rotation angle in deg
ratio=0.93; %ratio=energy_topimage/engergy_bottomimage: power ratio between the two imaging channels(representing nonideal beamsplitter)

%PSF-path
path='C:\Users\q004aj\Desktop\PSFs\';

%% loading "top-image" model (upper image on camera)

%load('./PSFs/PSF5D_tot(top)_0-2-250nm_RI=1,45_dz=0_aberrfree.mat'); %loading 3D or 5D PSF-model
load([path 'PSF5D_0-2-250nm_RI=1,45_dz=-400_aberrfree.mat']); %loading 3D or 5D PSF-model

two_ch='n'; %flag, indicating that two channels are simulated

if ndims(PSF5D)==3 %if "standard" 3D PSF is loaded 
    PSF_top=PSF_tot; 
    [nx0,ny0,nz0]=size(PSF_top); %size of model    
elseif ndims(PSF5D)==5 
    [nx0,ny0,nz0,nxi,nyi]=size(PSF5D);
    PSF_top=PSF5D(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)); %remove the unnecessary extra-dimensions
end

E_top=trapz(trapz(PSF_top,1),2);
clear PSF5D;
disp('done');


%% optional: loading "bottom-image" model
%if a second image model is loaded, two image stacks are created rather than just
%one

two_ch='y'; %flag, indicating that two channels are simulated

load([path 'PSF5D_UAF(bottom)_0-2-250nm_RI=1,45_dz=0_aberrfree.mat']); %loading 3D or 5D PSF-model

if ndims(PSF5D)==3 %if "standard" 3D PSF is loaded 
    PSF_bottom=PSF_tot;
elseif ndims(PSF5D)==5 
    PSF_bottom=PSF5D(:,:,:,ceil((nxi+1)/2),ceil((nyi+1)/2)); %remove the unnecessary extra-dimensions
end
E_bottom=trapz(trapz(PSF_bottom,1),2);

clear PSF5D;
disp('done');

%% creating coordinate systems

%creating real space coordinate system
x=linspace(0,(n_img-1)*ux,n_img);
xc=x(end/2); 
yc=xc; %center coords.
[X,Y]=meshgrid(x,x); %coordinate system on camera

%creating k-space coordinate system
uk=2*pi/ux/n_img;
Kx=X/ux*uk;
Ky=Y/ux*uk;

%% calculating molecule positions and creating image stacks
clear CCD_top CCD_bottom

no_mols=round(rho*(n_img*ux*1e6)^2); %total number of molecules in a single CCD image

for mm=1:no_images
    disp(no_images-mm)
    
    x_mol=rand(1,no_mols)*(n_img-1)*ux; 
    y_mol=rand(1,no_mols)*(n_img-1)*ux; 
    z_mol=min(r_sphere-sqrt(r_sphere.^2-(x_mol-xc).^2-(y_mol-yc).^2),250e-9);

%     figure(1);
%     scatter3(x_mol*1e6,y_mol*1e6,z_mol*1e9); %show molecule positions
%     xlabel('x/µm'); ylabel('y/µm');
%     zlabel('z/nm');
%     title('molecule positions - ground truth');

    if strcmp(two_ch,'y')
        %molecule positions in rotated & distorted second image: 
        tform = affine2d([cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1]);
        [xT_mol,yT_mol] = transformPointsForward(tform,x_mol,y_mol);
    end
    %generating molecule images from the PSF and moving them to the correct positions in the CCD image

    dz=z_vec(2)-z_vec(1); %z-increment
    I_top=zeros(nx0,nx0,no_mols);
    if strcmp(two_ch,'y')
        I_bottom=I_top;
    end
    
    for m=1:no_mols %generating images according ot z-positions
        I_top(:,:,m)=interpn(PSF_top,1:nx0,1:ny0,1+z_mol(m)/dz,'linear')*sig;
        
        if strcmp(two_ch,'y')
            I_bottom(:,:,m)=interpn(PSF_bottom,1:nx0,1:ny0,1+z_mol(m)/dz,'linear')*sig;
        end
        
    end
    disp('done');

    
    CCD_ideal_top=zeros(n_img,n_img);
    if strcmp(two_ch,'y')
        CCD_ideal_bottom=zeros(n_img,n_img);
    end
    
    for m=1:no_mols %moving molecules laterally
        tmp=embed(I_top(:,:,m),[n_img,n_img],0); %embedding molecules image in empty CCD image
        
        %A)Fourier-based shifting
        %F_tmp=fftshift(fft2((tmp))); 
        %phaseslope=exp(-1i*Kx*x_mol(m)-1i*Ky*y_mol(m)); %phase-factor for lateral shift
        %CCD_ideal_top=CCD_ideal_top+abs(fftshift(ifft2(ifftshift(F_tmp.*phaseslope))));

        %B) interp-based shifting
        CCD_ideal_top=CCD_ideal_top+interp2(X,Y,tmp,X-max(X(:))/2+x_mol(m),Y+max(Y(:))/2-y_mol(m),'linear',0);
                
        if strcmp(two_ch,'y')
            tmp2=embed(I_bottom(:,:,m),[n_img,n_img],0); %embedding molecules image in empty CCD image
            
            %A)Fourier-based shifting
            %F2_tmp=fftshift(fft2((tmp2))); 
            %phaseslope2=exp(-1i*Kx*xT_mol(m)-1i*Ky*yT_mol(m)); %phase-factor for lateral shift
            %CCD_ideal_bottom=CCD_ideal_bottom+flipud(abs(fftshift(ifft2(ifftshift(F2_tmp.*phaseslope2)))));
     
            %B) interp-based shifting
            CCD_ideal_bottom=CCD_ideal_bottom+flipud(interp2(X,Y,tmp2,X-max(X(:))/2+xT_mol(m),Y+max(Y(:))/2-yT_mol(m),'linear',0));
        end
        
    end

    % adding background and noise
    CCD_top(:,:,mm)=uint16(poissrnd((CCD_ideal_top+BG_top)*ratio)*gain/amp*QE);
    
    if strcmp(two_ch,'y')
        CCD_bottom(:,:,mm)=uint16(poissrnd(CCD_ideal_bottom+BG_bottom)*gain/amp*QE);
    end
    
    
    %showing images 
    if strcmp(two_ch,'y')
        
        subplot(2,1,1);
        imagesc(CCD_top(:,:,mm)); axis equal; axis tight; colormap gray; colorbar;
        title('top image');
        subplot(2,1,2);
        imagesc(CCD_bottom(:,:,mm)); axis equal; axis tight; colormap gray; colorbar;
        title('bottom image');
        pause(0);
    else
        figure(2);
        imagesc(CCD_top(:,:,mm)); 
        axis equal; axis tight; colormap gray; colorbar;
        title('top image');
        pause(0);
    end
    
end

%% saving data

filename=['sig-bg=' num2str(sig,4) '-' num2str(BG_top,3) '_gain=' num2str(gain)  '_rho=' num2str(rho)];
for mm=1:no_images
    imwrite(CCD_top(:,:,mm),['simu_top_(tot)_' filename '.tif'],'WriteMode','Append');
    
    if strcmp(two_ch,'y')
        imwrite(CCD_bottom(:,:,mm),['simu_bottom_(UAF)_' filename '.tif'],'WriteMode','Append');
    end
end
disp('done');

