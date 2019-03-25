function metric=fun_optPSFforSAF(a)
%a...vector of Zernikes

global Z ux uk uz n_photon bg Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz Nx mask 

%------------------------------

a2=zeros(1,1,length(a));
a2(1,1,:)=a;
phase=sum(repmat(a2,[size(Z,1),size(Z,2),1]).*Z,3); %pupil phase calculation
imagesc(phase); title('pupil phase'); pause(0.1);

for m=1:size(Ex_Px,3)
    I_xx=abs(czt2(Ex_Px(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    PSF(:,:,m)=I_xx+I_yx+I_xy+I_yy+I_xz+I_yz;
    PSF(:,:,m)=PSF(:,:,m)/sum(sum(PSF(:,:,m))); %normalization of PSF
end
figure(2);
imagesc(PSF(:,:,1)); 

figure(3);
[CRBx,CRBy,CRBz]=fun_CRB(PSF,ux,uz,n_photon,bg,1);
plot(sqrt(CRBz)); pause(0);

metric1=mean(sqrt((CRBx.*CRBy.*CRBz)).^(1/3));  %"localization volume"
metric2=mean(sqrt(CRBz)); 

metric=metric2;