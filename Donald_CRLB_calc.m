%CRLB calcutations for paper 

clear all; 

%loading PSFs
tmp=load('PSF_NA1.49_15x15_0-10-250nm_RI=1.33_dz=0nm_aberr=n.mat')
PSF1=tmp.PSF;
cam=tmp.cam;
tmp=load('PSF_UAF_NA1.49_15x15_0-10-250nm_RI=1.33_dz=0nm_aberr=n.mat')
PSF2=tmp.PSF; 

%calculating precisions
sig=2000;
bg=100; 

%% adding Fisher information of both channels
split_ratio=0.5; 

[~,~,~,~,~,FI1]=PSF1.CRLB(sig,bg,cam,split_ratio); %here it is important to set the detection efficiency of the channel!
[~,~,~,~,~,FI2]=PSF2.CRLB(sig,bg,cam,(1-split_ratio)); 
FI=FI1+FI2;

for m=1:size(FI,3)
    %CRBs calculated from the sum of Fisher infos -> optimal
    tmp=inv(FI(:,:,m));
    CRBx(m)=tmp(1,1);
    CRBy(m)=tmp(2,2);
    CRBz(m)=tmp(3,3);
    CRBsig(m)=tmp(4,4);
    CRBbg(m)=tmp(5,5);
end

figure(4);
plot((0:length(CRBz)-1)*PSF1.uz*1e9,sqrt(CRBz),'r.-'); xlabel('z-pos in nm'); ylabel('nm');
title(['sqrt(CRBz), cts=' num2str(sig) ', bg=' num2str(bg)]); grid on;
hold on; 
plot((0:length(CRBz)-1)*PSF1.uz*1e9,sqrt(CRBx)); xlabel('z-pos in nm'); ylabel('nm');
plot((0:length(CRBz)-1)*PSF1.uz*1e9,sqrt(CRBy),'b.-'); xlabel('z-pos in nm'); ylabel('nm');
title(['x,y,z precisions(z=red), cts='  num2str(sig) ', bg=' num2str(bg)]);
ylim([0 50]);
title('DONALD');