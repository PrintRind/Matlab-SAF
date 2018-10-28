
%simulation of supercritical-angle-fluorescence detection
%aim: z-localization by measuring SAF/UAF ratio

close all;
clear all;

%%

dia_pupil=128;
lambda_0=532e-9;
ux_desired=40e-9;
NA=1.49;
n1=1.35; 
n2=1.52;
NA_crit=n1; 

[Kx,Ky,Kr,pupil,defocus,N_pad,uk,ux,v]=fun_usualsuspects(dia_pupil,lambda_0,ux_desired,NA,n2);

SAF_mask=(Kr>NA_crit/NA*Kr(end/2,end)).*pupil;
UAF_mask=not(SAF_mask).*pupil;



M=20;
I_UAF=linspace(0.5,1,M);
I_SAF=1-I_UAF;

figure(1);

for m=1:M;

%I_UAF=0.8;

Int=UAF_mask*I_UAF(m)+SAF_mask*I_SAF(m); %imagesc(Int);

phase_SAF=pi;
phase_UAF=0; %atan2(Ky,Kx);%atan2(Ky,Kx);
holo=UAF_mask.*exp(1i*phase_UAF)+SAF_mask.*exp(1i*phase_SAF);
PSF=abs(fftshift(fft2(ifftshift(embed(sqrt(Int).*holo,N_pad,0))))).^2;

subplot(1,2,1);
imagesc(Int,[0, 1]); colormap gray;

subplot(1,2,2);
imagesc(PSF(v,v)); pause(0.1);

w(m)=fun_FWHM(PSF(round(end/2),round(end/2)-10:round(end/2)+10),0.5)*ux;
I_max(m)=max(PSF(:));
%plot(PSF(round(end/2),round(end/2)-5:round(end/2)+5)); pause(0.01);

end
figure(1);
plot(I_UAF,w); xlabel('I_{UAF} (I_{SAF}=1-I_{UAF})');

figure(2);
plot(I_UAF,I_max); xlabel('I_{UAF} (I_{SAF}=1-I_{UAF})'); ylabel('max PSF intenstiy');










%% radiation model described by Enderlein in Chem.PhysLett 308 (1999).

n1=1.35;
n2=1.52;
k_0=2*pi/lambda_0;
c=3e8;

%THETA2=80/180*pi;
%PHI=pi/4;

%function ddS=emission(k_0,n1,n2,THETA2,K_p2,dip,K_s,z_0)
%k_0...vacuum wave vector
%THETA2...angle of emission direction (into glass)
%n1...refractive index of solution
%n2...refractive index of glass
%dip...emission dipole (vector)
%z_0...distance of dipole from interface
%K_p2, K_s....unit vectors; dependent on PHI and THETA


%coordinates
N=128;
[X,Y,R,mask]=create_coord(N,2/N,'exact');
PHI=atan2(Y,X);
THETA2=(asin(R)).*mask;


z_0=0;   %distance of molecule from interface in units of wavelengths
dipole=[0;1;0];  %emission dipole orientation


    w=sqrt(n1^2-n2^2*sin(THETA2).^2); %regions where w is imaginary correspond to the "normal" fluorescence
    
    K_p2=[w.*cos(PHI), w.*sin(PHI), n2*sin(THETA2)]/n1;   
        K_p2_res=ones(N,N,3);   %rearrangeing the matrix to a tensor 3 x N x N
        K_p2_res(:,:,1)=K_p2(:,1:N);
        K_p2_res(:,:,2)=K_p2(:,(1:N)+N);
        K_p2_res(:,:,3)=K_p2(:,(1:N)+2*N);
        
    K_s=[-sin(PHI), cos(PHI), zeros(N,N)];
        K_s_res=ones(N,N,3);   %rearrangeing the matrix to a tensor 3 x N x N
        K_s_res(:,:,1)=K_s(:,1:N);
        K_s_res(:,:,2)=K_s(:,(1:N)+N);
        K_s_res(:,:,3)=K_s(:,(1:N)+2*N);

    
    T_p=(2*n1*sqrt(n1^2-n2^2.*sin(THETA2).^2))./(n2*sqrt(n1^2-n2^2*sin(THETA2).^2)+n1^2*cos(THETA2)); %Fresnel transmission
    T_s=(2*sqrt(n1^2-n2^2*sin(THETA2).^2))./(sqrt(n1^2-n2^2*sin(THETA2).^2)+n2*cos(THETA2));

    T_p_stack=repmat(T_p,[1,1,3]);  
    T_s_stack=repmat(T_s,[1,1,3]); 
    
    dipole_stack=permute(repmat(dipole,[1,N,N]),[2,3,1]);
    
    ddS=c*k_0^4*n2/8/pi*abs(cos(THETA2)./w).^2.*(abs(sum(T_p_stack.*K_p2_res.*dipole_stack,3)).^2+abs(sum(T_s_stack.*K_s_res.*dipole_stack,3)).^2).*exp(-2*imag(w).*z_0);
  
    figure(1);
    imagesc(ddS); axis equal; colorbar; xlabel('x'); ylabel('y');
    figure(2);
    plot(X(end/2,:),ddS(end/2,:));
    
    