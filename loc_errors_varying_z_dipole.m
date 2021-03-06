%% calculating systematic localization errors caused by varying z-dipole or aberration contributions
%-)first, the PSF used for data evaluation is chosen/loaded
%-)then, molecule images are calculated based on this PSF, but with a varying
%parameter regarding dipole-emission of Zernike aberration 
%-)the molecule images are evaluated with the loaded PSF and the systematic
%errors are displayed. 

clear all; 

N=128; %pupil diameter in pixels 
Nx=13; %field size in focal plane

mode='single'; %'single' or 'biplane' or 'donald'
PSF_type='defocus';  %(important only for file/figure-naming); define type of PSF that is used for evaluation "Ast" or "Defocus"

%% --- define PSF

lens=objective('olympus 1.49');
cam=camera('orca fusion');
T=lens.transmission(N); 

%----loading PSF used for the estimation
[PSF1_name,PSF1_path,~]=uigetfile('*.mat','choose PSF file for 1st channel.');
tmp=load([PSF1_path PSF1_name]); 
PSF1=tmp.PSF;

%optionally, load a 2nd PSF:
if (strcmp(mode,'biplane')) || (strcmp(mode,'donald'))
    [PSF2_name,PSF2_path,~]=uigetfile('*.mat','choose PSF file for 2nd channel.');
    tmp=load([PSF2_path PSF2_name]); 
    PSF2=tmp.PSF;
end

%% calculating molecule images with varying parameters (e.g. aberrations or dipole orientation etc.)

noise='n'; 
%----adding noise----
sig=2000; 
bg=100; 

z_range=(0:PSF1.uz:PSF1.uz*(PSF1.Nz-1)); %z-working range for localization

%----define parameter to be varied-----------------------------------------

%marechal criterion is 72 mlambda
no_vars=11; %number of variations
var_parameter='Z11'; var_range=linspace(-72e-3*2*pi,72e-3*2*pi,no_vars);%e.g. 'Z5', 'Z1',etc. or 'z-dipole strength'
%var_parameter='Z8'; var_range=linspace(0,0.25,no_vars);%e.g. 'Z5', 'Z1',etc. or 'z-dipole strength'
%var_parameter='z-dipole strength'; var_range=linspace(0,1,no_vars); %e.g. 'Z5', 'Z1',etc. or 'z-dipole strength'

if exist('PSF2') %biplane imaging
    
    if strcmp(mode,'donald')
        %identify which PSF is the UAF-PSF: 
        if sum(sum(PSF1.data(:,:,1))) < sum(sum(PSF2.data(:,:,1)))
            mode1='UAF'; 
            mode2='';
        else
            mode2='UAF';
            mode1='';
        end
    else 
        mode1='';
        mode2='';
    end
        
    [I_test_A, FI1]=fun_varyParam_calcImages(N,var_parameter,var_range,z_range,PSF1,lens,cam,noise,sig/2,bg/2,mode1);%calculate test images with varying parameter
    [I_test_B, FI2]=fun_varyParam_calcImages(N,var_parameter,var_range,z_range,PSF2,lens,cam,noise,sig/2,bg/2,mode2);%calculate test images with varying parameter
    PSF=[PSF1,PSF2];
    
    for v=1:length(FI1)
        FI=FI1{v}+FI2{v}; 
        
        for m=1:size(FI,3) 
            tmp=inv(FI(:,:,m));
            CRBx(m,v)=tmp(1,1);
            CRBy(m,v)=tmp(2,2);
            CRBz(m,v)=tmp(3,3);
            CRBsig(m,v)=tmp(4,4);
            CRBbg(m,v)=tmp(5,5);
        end
    end
    
else %single-channel imaging
    mode='single';
    [I_test, FI]=fun_varyParam_calcImages(N,var_parameter,var_range,z_range,PSF1,lens,cam,noise,sig,bg);%calculate test images with varying parameter
    PSF=PSF1;
    
        for v=1:length(FI)
            FI_tmp=FI{v};
            for m=1:size(FI_tmp,3) 
                tmp=inv(FI_tmp(:,:,m));
                CRBx(m,v)=tmp(1,1);
                CRBy(m,v)=tmp(2,2);
                CRBz(m,v)=tmp(3,3);
                CRBsig(m,v)=tmp(4,4);
                CRBbg(m,v)=tmp(5,5);
            end
    end
end

%------------estimate----------------

z_ini=mean(z_range);
clear x_est y_est z_est N_est BG_est I error
for v=1:no_vars
    
    if (strcmp(mode,'biplane')) || (strcmp(mode,'donald'))
        I{1}=I_test_A(:,:,:,v);
        I{2}=I_test_B(:,:,:,v);
    elseif strcmp(mode,'single')
        I=I_test(:,:,:,v);
    end
        
    est=fun_MLE3(PSF,I,z_ini,cam); 
    x_est(v,:)=est(:,1)';
    y_est(v,:)=est(:,2)';
    z_est(v,:)=est(:,3)';
    N_est(v,:)=est(:,4)';
    BG_est(v,:)=est(:,5)';
    LLR(v,:)=est(:,6)'; %log-likelihood ratio (after Huang et al, suppl. mat.), should be zero for perfect fit
end
disp('done, ready for display');
    

% ----- display results------------------------------------------------
%---------------------------------------------------------------------

zlim=[-150, 150]; 
xylim=[-100, 100]; 
siglim=[-1000 1000]; 
bglim=[-50,50];
errlim=[0 50];

legend_param=num2str(var_range);
z_error=z_est-repmat(z_range,[no_vars,1]);
x_error=x_est; 
y_error=y_est; 


%------------display biases----------------------
figure(2); 
%sgtitle(['varying: ' var_parameter '; z-bias']);
subplot(3,2,1); 
plot(z_range*1e9,(z_est-z_range)*1e9,'-');
xlabel('z_{true} / nm');
ylabel('z_{est}-z_{true} / nm');
grid on; 
title('z-bias');
xlim([0 max(z_range)*1e9]);
ylim(zlim);

subplot(3,2,2); 
plot(z_range*1e9,(x_est)*1e9,'--');
hold on; 
plot(z_range*1e9,(y_est)*1e9,'-');
hold off; 
xlabel('z_{true} / nm');
ylabel('x_{est} / y_{est} / nm');
grid on; 
title('xy-bias ');
xlim([0 max(z_range)*1e9]);
ylim(xylim);

subplot(3,2,3); 
plot(z_range*1e9,N_est-sig);
xlabel('z_{true} / nm');
ylabel('N_{est}-N_{true}');
grid on; 
title(['Signal bias; sig=' num2str(sig)]);
xlim([0 max(z_range)*1e9]);
ylim(siglim);

subplot(3,2,4); 
plot(z_range*1e9,BG_est-bg);
xlabel('z_{true} / nm');
ylabel('BG_{est}-BG_{true}');
grid on; 
title(['BG bias; BG=' num2str(bg)]);
xlim([0 max(z_range)*1e9]);
ylim(bglim);

subplot(3,2,5); 
plot(z_range*1e9,LLR,'-');
xlabel('z_{true} / nm');
ylabel('error');
grid on; 
title('log-likelihood ratio');
xlim([0 max(z_range)*1e9]);
%ylim(errlim);

%---plotting CRLBs for the simulated image stacks---
%this gives an idea on the obtainable preision if there were NO model
%mismatch
clear hugo

subplot(3,2,6);
for v=1:no_vars
    plot(z_range*1e9,sqrt(CRBx(:,v)),'b', z_range*1e9,sqrt(CRBy(:,v)),'g',z_range*1e9,sqrt(CRBz(:,v)),'r'); 
    hold on; 
end
xlabel('z / nm');
ylabel('\sigma / nm')
title('CRLB-precision if no model-mismatch');
xlim([0 max(z_range)*1e9]); 
grid on; 
hold off; 

legend(legend_param);

%adding CRLB info if noise was included
if strcmp(noise,'y')
    if strcmp(mode,'single')
        [Cx,Cy,Cz,CN,Cbg,~]=PSF.CRLB(sig,bg,cam,1);
    elseif strcmp(mode,'biplane') || strcmp(mode,'donald')
        [~,~,~,~,~,FI1]=PSF1.CRLB(sig,bg,cam,0.5);
        [~,~,~,~,~,FI2]=PSF2.CRLB(sig,bg,cam,0.5);
        FI=FI1+FI2; 
        for m=1:PSF1.Nz
            tmp=inv(FI(:,:,m));
            Cx(m)=tmp(1,1);
            Cy(m)=tmp(2,2);
            Cz(m)=tmp(3,3);
            CN(m)=tmp(4,4);
            Cbg(m)=tmp(5,5);
        end
    end
    subplot(3,2,1); 
    hold on; 
    plot(z_range*1e9,sqrt([Cz]),'r--',z_range*1e9,-sqrt([Cz]),'r--');
    hold off; 
    
    subplot(3,2,2); 
    hold on; 
    plot(z_range*1e9,sqrt([Cx]),'r--',z_range*1e9,-sqrt([Cx]),'r--');
    hold off; 
    
    subplot(3,2,3); 
    hold on; 
    plot(z_range*1e9,sqrt([CN]),'r--',z_range*1e9,-sqrt([CN]),'r--');
    hold off; 

    subplot(3,2,4); 
    hold on; 
    plot(z_range*1e9,sqrt([Cbg]),'r--',z_range*1e9,-sqrt([Cbg]),'r--');
    hold off; 
end

mtit([PSF_type ', varying ' var_parameter],'fontsize',12,'xoff',0,'yoff',0.04);

hugo=legend; 
hugo.String{1}='x'; 
hugo.String{2}='y'; 
hugo.String{3}='z'; 
hugo.String(4:end)=[];
hold off; 

%mtit([PSF_type ', varying ' var_parameter '=' num2str(var_range)],'fontsize',12,'xoff',0,'yoff',0.04);

disp('all done, ready for saving data'); 

% ----calculate mean squared errors---------------------------------------

Ez=sqrt(sum((z_est-z_range).^2,2)/length(z_range)) %root mean squared z-error over the given z-range
r_est=sqrt(x_est.^2+y_est.^2);
Ex=sqrt(sum(x_est.^2,2)/length(z_range)) %root mean squared x-error over the given z-range
Ey=sqrt(sum(y_est.^2,2)/length(z_range)) %root mean squared y-error over the given z-range
Esig=sqrt(sum((N_est-sig).^2,2)/length(z_range))/sig*100; %in percent 
Ebg=sqrt(sum((BG_est-bg).^2,2)/length(z_range))/bg*100;  %in percent

figure(3); 
subplot(2,1,1);
plot(var_range,[Ez, Ex, Ey]*1e9,'o-');
xlabel('parameter');
ylabel('RMS error / nm');
legend('Ez','Ex','Ey');
subplot(2,1,2);
plot(var_range,[Esig, Ebg],'o-');
xlabel('parameter');
ylabel('RMS error / percent');
legend('Esig','Ebg');

%% saving data

%figure(2);
save([cd '\biases\' PSF_type '_' var_parameter '.mat'],'x_est','y_est','z_est','N_est','BG_est','error','CRBx','CRBy','CRBz','var_range','sig','bg','z_range');
% figure(3); 
% savefig([PSF_type '_' var_parameter '_errors.fig']);
% figure(5);
% savefig([PSF_type '_' var_parameter '_CRLBs.fig']);
