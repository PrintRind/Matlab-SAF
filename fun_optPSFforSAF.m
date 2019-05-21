function metric=fun_optPSFforSAF(a)
%a...vector of Zernikes


global Z ux uk uz n_photon bg Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz Nx mask E_tot

%------------------------------

a2=zeros(1,1,length(a));
a2(1,1,:)=a;
%phase=sum(repmat(a2,[size(Z,1),size(Z,2),1]).*Z,3); %pupil phase calculation
phase=sum(Z.*a2,3); 
figure(2); 
subplot(2,1,1); 
imagesc(phase); title('pupil phase'); 

for m=1:size(Ex_Px,3)
    I_xx=abs(czt2(Ex_Px(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:,m).*exp(1i*phase).*mask,uk,ux,Nx)).^2;
    PSF(:,:,m)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m);
end

figure(2);
subplot(2,1,2); 
nor=0; %no normalization of PSFs
[CRBx,CRBy,CRBz]=fun_CRB(PSF,ux,uz,n_photon,bg,1,nor);
plot(sqrt(CRBz(2:end))); pause(0);

%choose metric for optimization

metric1=mean(sqrt((CRBx(2:end).*CRBy(2:end).*CRBz(2:end))).^(1/3));  %"localization volume"
metric2=mean(sqrt(CRBz(2:end))); 

metric=metric2;

