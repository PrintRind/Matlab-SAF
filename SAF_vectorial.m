%Super-critical angle fluorescence
%Calculation of PSFs for emitting dipoles near a high-refractive index
%surface (near-field coupling and polarization issues are included)
%based on Axelrod J. of Microscopy 2012, vol. 247

clear all; 
close all;
%run('C:\Program Files\DIPimage 2.4.1\dipstart.m'); %run dipimage

%% user parameters

lambda_0=532e-9; %wavelength
n1=1.33; %refractive index of solution
n2=1.42; %of layer on glass (e.g. lipid bilayer)
n3=1.518; %of glass & immersion oil
f=1.8e-3; %focal length of microscope objective
NA=1.49; %NA of objective
h=0.0*lambda_0; %thickness of intermediate layer
scalemax=.3e37; %max scalebar  

%coordinates in the objective pupil
N=64;
[X,Y,R,mask]=create_coord(N,2/N,'exact');
PHI3=atan2(Y,X);
mask_obj=R<(NA/n3); %aperture of the objective
mask_UAF=R<(n1/n3);

%----------SLM phase masks------------------

phi_rel=0; %relative phase shift between UAF and SAF light (this value applied to SAF-light)
phase_all=zernike_modes(N,N,[0 0 0 0 0]);%corr; %phasemask applies to total BFP (UAF+SAF)
%phase_all=0;    
%corr=-angle(E_BFP_x);

%----------define dipole position and orientation
z=(0:5:25)*1e-9; 
 
theta_mu=0;  %theta-angle of dipole (0=parallel to z-axis, pi/2...in x-y-plane)
phi_mu=0;    %phi-angle of dipole
mu=1; %strength of dipole (field)

%----------polarized detection (set px and py to 1 for unpolarized detection)
px=1;
py=1;

%----------------------------------------

%calculations
k0=2*pi/lambda_0;
k1=k0*n1;
k2=k0*n2;
k3=k0*n3;

THETA3=asin(R).*mask; %angle in medium 3  (R=-1...1; THETA3=-90...90°)
THETA2=acos(sqrt(1-(n3/n2*R).^2)); %angle in medium2
THETA1=acos(sqrt(1-(n3/n1*R).^2)); %angle in medium1

%calculations according to paper of Axelrod, 2012
   
%Fresnel-coefs  %Eq. 3 of paper
    tp12=2*n1*cos(THETA1)./(n1*cos(THETA2)+n2*cos(THETA1));
    tp23=2*n2*cos(THETA2)./(n2*cos(THETA3)+n3*cos(THETA2));
    ts12=2*n1*cos(THETA1)./(n1*cos(THETA1)+n2*cos(THETA2));
    ts23=2*n2*cos(THETA2)./(n2*cos(THETA2)+n3*cos(THETA3));
    rp12=(n2*cos(THETA1)-n1*cos(THETA2))./(n1*cos(THETA2)+n2*cos(THETA1));
    rp23=(n3*cos(THETA2)-n2*cos(THETA3))./(n2*cos(THETA3)+n3*cos(THETA2));
    rs12=(n1*cos(THETA1)-n2*cos(THETA2))./(n1*cos(THETA1)+n2*cos(THETA2));
    rs23=(n2*cos(THETA2)-n3*cos(THETA3))./(n2*cos(THETA2)+n3*cos(THETA3));

    %from them we calculate Fresnel coefs for three-layer system
    tp=(tp12.*tp23.*exp(1i*k2*h*cos(THETA2)))./(1+rp12.*rp23.*exp(2i*k2*h*cos(THETA2)));
    ts=(ts12.*ts23.*exp(1i*k2*h*cos(THETA2)))./(1+rs12.*rs23.*exp(2i*k2*h*cos(THETA2)));

    %dipole projections onto directions s, p and z:
    mu_p=mu*sin(theta_mu).*cos(phi_mu-PHI3);
    mu_s=mu*sin(theta_mu).*sin(phi_mu-PHI3);
    mu_z=mu*cos(theta_mu);
   
for m=1:length(z); %defocus-steps

    C=(k3^2*exp(1i*k3*f).*cos(THETA3))/f/n1.*exp(-1i*k3*h*cos(THETA3)).*exp(1i*k1.*cos(THETA1).*z(m)); %eq.11 in paper

    %field magnitudes in layer 3 (pre-objective zone), along the s,p and z-axes
    E3p=C.*tp.*cos(THETA3).*(mu_p./n3+mu_z.*sin(THETA3)./cos(THETA1));
    E3s=C.*ts.*(mu_s/n3./cos(THETA1));
    E3z=C.*tp.*sin(THETA3).*(mu_p./n3+mu_z.*sin(THETA3)./cos(THETA1));

    %influence of objective: rotation of rays by their angle theta3, such that they are all parallel to
    %the optical axis: 
    Apodization=1./sqrt(cos(THETA3)).*mask_obj; %Apodization of objective lens
    E_BFP_p=(E3p.*cos(THETA3)+E3z.*sin(THETA3)).*Apodization;
    E_BFP_s=E3s.*Apodization;  %s-polarization remains unchanged by this rotation

    %coordinate transform into x-and y-polarization
    E_BFP_x=cos(PHI3).*E_BFP_p-sin(PHI3).*E_BFP_s;
    E_BFP_y=sin(PHI3).*E_BFP_p+cos(PHI3).*E_BFP_s;
    E_x(:,:,m)=E_BFP_x; %for saving as mat-file-stack
    E_y(:,:,m)=E_BFP_y; %for saving as mat-file-stack

    I_BFP_x=abs(E_BFP_x).^2;
    I_BFP_y=abs(E_BFP_y).^2;

        
    % display BFP intensities and phases
    figure(1);
    colormap gray;
    subplot(3,2,1);
    imagesc(abs(E_BFP_x)); axis equal; axis tight; title('|E_x|'); %colorbar;
    subplot(3,2,2);
    imagesc(angle(E_BFP_x)+phase_all); axis equal; axis tight; title('angle(Ex)+phase_{all}'); %colorbar;
    subplot(3,2,3);
    imagesc(abs(E_BFP_y)); axis equal; axis tight; title('|E_y|'); %colorbar;
    subplot(3,2,4);
    imagesc(angle(E_BFP_y)+phase_all); axis equal; axis tight; title('angle(Ey)+phase_{all}'); %colorbar;
    subplot(3,2,5);
    imagesc(abs(E_BFP_y)+abs(E_BFP_x)); axis equal; axis tight; title('|E_x|+|E_y|'); %colorbar;
    %subplot(3,2,6);
    
    figure(6);  
    v=6;
    quiver(X(1:v:end,1:v:end),Y(1:v:end,1:v:end),E_BFP_x(1:v:end,1:v:end),E_BFP_y(1:v:end,1:v:end));
    axis equal; 
    title('Electric field in BFP');

    
    % measure UAF and SAF intensities
%     SAF_image=(I_BFP_x+I_BFP_y).*not(mask_UAF);
%     SAF(m)=sum(SAF_image(:)); %sum of SAF intensity;
%     
%     UAF_image=(I_BFP_x+I_BFP_y).*mask_UAF;
%     UAF(m)=sum(UAF_image(:)); %sum of UAF intensity;
%     UAFSAF(m)=UAF(m)+SAF(m);
%     ratio(m)=SAF(m)/UAF(m);
    
    phasemask=((mask_obj+not(mask_UAF))-1).*phi_rel;

    % -----------------------------------------------------------------------
    % PSF calculations

    N_pad=2*N; %padding array to increase resolution in x-space
    uk=2*k3/N; %k-space unit
    ux=2*pi/uk/N_pad; %x-space unit
    
    fov=10; %field of view in x-space
    
    %UAF and UAF+SAF PSFs
        E_x_SAF=fftshift(fft2(ifftshift(embed(E_BFP_x.*not(mask_UAF).*exp(1i*phase_all),N_pad,0))))/N_pad;
        E_y_SAF=fftshift(fft2(ifftshift(embed(E_BFP_y.*not(mask_UAF).*exp(1i*phase_all),N_pad,0))))/N_pad;
        I_x_SAF=abs(E_x_SAF(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2; %intensities on camera
        I_y_SAF=abs(E_y_SAF(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2;
        I_tot_SAF=px*I_x_SAF+py*I_y_SAF;
        
        E_x_UAF=fftshift(fft2(ifftshift(embed(E_BFP_x.*mask_UAF.*exp(1i*phase_all),N_pad,0))))/N_pad;
        E_y_UAF=fftshift(fft2(ifftshift(embed(E_BFP_y.*mask_UAF.*exp(1i*phase_all),N_pad,0))))/N_pad;
        I_x_UAF=abs(E_x_UAF(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2; %intensities on camera
        I_y_UAF=abs(E_y_UAF(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2;
        I_tot_UAF=px*I_x_UAF+py*I_y_UAF;

        E_x_tot=fftshift(fft2(ifftshift(embed(E_BFP_x.*exp(1i*phase_all),N_pad,0))))/N_pad;
        E_y_tot=fftshift(fft2(ifftshift(embed(E_BFP_y.*exp(1i*phase_all),N_pad,0))))/N_pad;
        I_x_tot=abs(E_x_tot(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2; %intensities on camera
        I_y_tot=abs(E_y_tot(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2;
        I_tot=px*I_x_tot+py*I_y_tot;
    
    %PSF resulting from phaseshift between UAF/SAF light

        E_x_phi=fftshift(fft2(ifftshift(embed(E_BFP_x.*exp(1i*phasemask+1i*phase_all),N_pad,0))))/N_pad;
        E_y_phi=fftshift(fft2(ifftshift(embed(E_BFP_y.*exp(1i*phasemask+1i*phase_all),N_pad,0))))/N_pad;
        I_x_phi=abs(E_x_phi(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2; %intensities on camera
        I_y_phi=abs(E_y_phi(end/2-fov:end/2+fov,end/2-fov:end/2+fov)).^2;
        I_tot_phi=px*I_x_phi+py*I_y_phi;

        %UAF and total light
        figure(2);
        subplot(3,3,1);
        imagesc(I_x_tot,[0 scalemax]); axis equal; axis tight;  title('I_{x, tot}'); %colorbar;
        subplot(3,3,4);
        imagesc(I_y_tot,[0 scalemax]); axis equal; axis tight; title('I_{y, tot}'); %colorbar;
        subplot(3,3,7);
        imagesc(I_tot,[0 scalemax]);
        colormap gray; axis equal; axis tight; title('total intensity'); %colorbar;
        subplot(3,3,2);
        imagesc(I_x_UAF,[0 scalemax]); axis equal; axis tight;  title('I_{x, UAF}'); %colorbar;
        subplot(3,3,5);
        imagesc(I_y_UAF,[0 scalemax]); axis equal; axis tight; title('I_{y, UAF}'); %colorbar;
        subplot(3,3,8);
        imagesc(I_tot_UAF,[0 scalemax]); axis equal; axis tight; title('total UAF intensity'); %colorbar;
        subplot(3,3,3);
        imagesc(I_x_phi,[0 scalemax]); axis equal; axis tight;  title('I_{x, phi}'); %colorbar;
        subplot(3,3,6);
        imagesc(I_y_phi,[0 scalemax]); axis equal; axis tight; title('I_{y, phi}'); %colorbar;
        subplot(3,3,9);
        imagesc(I_tot_phi,[0 scalemax]); axis equal; axis tight; title('total intensity with phi'); %colorbar;

        %SAF-light
        figure(3); 
        colormap gray;
        subplot(1,3,1)
        imagesc(I_x_SAF,[0 scalemax]); axis equal; axis tight;  title('I_{x, SAF}'); %colorbar;
        subplot(1,3,2);
        imagesc(I_y_SAF,[0 scalemax]); axis equal; axis tight; title('I_{y, SAF}'); %colorbar;
        subplot(1,3,3);
        imagesc(I_tot_SAF,[0 scalemax]); axis equal; axis tight; title('total SAF intensity'); %colorbar;        
   
   %---measuring total intensities in the focal plane----
   
   %defining summation ROI
   cx=12; cy=12; %center coords.         
   HW=2; %half widths
   
   UAF(m)=sum(sum(I_tot_UAF(cy-HW:cy+HW,cx-HW:cx+HW))); %total PSF intensity for UAF-light
   TOT(m)=sum(sum(I_tot(cy-HW:cy+HW,cx-HW:cx+HW)));     %total PSF intensity for UAF+SAF-light
   PHI(m)=sum(sum(I_tot_phi(cy-HW:cy+HW,cx-HW:cx+HW))); %total PSF intensity when phaseshift is applied to SAF light
   ratio(m)=UAF(m)/TOT(m);
   ratio_phi(m)=PHI(m)/TOT(m);
   
   UAF_x(m)=sum(sum(I_x_UAF(cy-HW:cy+HW,cx-HW:cx+HW)));
   TOT_x(m)=sum(sum(I_x_tot(cy-HW:cy+HW,cx-HW:cx+HW)));
   PHI_x(m)=sum(sum(I_x_phi(cy-HW:cy+HW,cx-HW:cx+HW)));

   UAF_y(m)=sum(sum(I_y_UAF(cy-HW:cy+HW,cx-HW:cx+HW)));
   TOT_y(m)=sum(sum(I_y_tot(cy-HW:cy+HW,cx-HW:cx+HW)));
   PHI_y(m)=sum(sum(I_y_phi(cy-HW:cy+HW,cx-HW:cx+HW)));

end

uz=z(2)-z(1);
%saving data for evaluation with "optimal_SAF_PSF.m"
save(['SAF_theta=' num2str(theta_mu) '_phi=' num2str(phi_mu) '.mat'],'E_x','E_y','ux','uk','uz','NA','lambda_0','n1','n2','n3','theta_mu','phi_mu','mask_obj','mask_UAF','z');

figure(4);
subplot(1,2,1);
    plot(z*1e9,UAF,'b-','LineWidth',2); xlabel('z / nm (distance dipole to surface)'); ylabel(['integr. Int. over ' num2str(2*HW+1) 'x' num2str(2*HW+1) 'pix']); 
    title(['dipole-theta=' num2str(theta_mu/pi) 'pi; UAF+SAF(red), UAF(blue), phi(green)']);
    hold on;
    plot(z*1e9,PHI,'g-','LineWidth',2); 
    plot(z*1e9,TOT,'r-','LineWidth',2); 
    %plot(z,ratio,'bo--');
    hold off;
    xlim([0 100])
    ylim([0.5e37 5e37]) 
subplot(1,2,2);
    plot(z,UAF_x,'bo-'); xlabel('z (distance dipole from surface)'); ylabel(['integrated Int. over ' num2str(2*HW+1) 'x' num2str(2*HW+1) 'pixel']); title(['X-POL and Y-POL (dashed) separately']);
    hold on;
    plot(z,PHI_x,'go-'); 
    plot(z,TOT_x,'ro-'); 
    plot(z,UAF_y,'bo--'); 
    plot(z,PHI_y,'go--'); 
    plot(z,TOT_y,'ro--'); 
    hold off;

figure(5);
plot(z,ratio,'bo-');
hold on ;
plot(z,ratio_phi,'go-');
title('UAF/TOT ratios for the normal case (blue) and the interference case (green)');
xlabel('z');
ylabel('ratio');
hold off;
