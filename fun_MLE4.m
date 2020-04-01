function est=fun_MLE4(PSF,data,z_ini,CAM)
%estimation based on 3D PSF; PSF should be finer sampled than the
%effective camera pixel size in focal plane
%inputs:
%--------------
%PSF....member of the "psf" class; If PSF is an array containing two
%PSF-class members --> biplane imaging mode
%data....structure array containing the following: 
    %data.images.... stack of single molecule imgages (e.g. size 13x13x10000)
    %data.ULC....array containing the upper left corner coordinates of the
    %molecule images with respect to the sCMOS sensor
%If PSF is a vector containing two PSFs --> biplane imaging mode
    %then, I_stack must be a structure array with two entries 
%z_ini_ind....initial estimate for molecule z position in meter
%cam...member of the camera class 
%if data contains only a stack of images (and no pixel-dependent
%noise/gain/baseline info), the properties are assumed to be uniform over
%the entire camera chip

z_ini_ind=round(z_ini/PSF(1).uz+1);
fw=(PSF(1).Nx/PSF(1).os-size(I_stack,1))/2; %frame width

%creating large grid for oversampled PSF
x=((PSF(1).os+1)/2)+fw*PSF(1).os:PSF(1).os:PSF(1).Nx-((PSF(1).os-1)/2)-fw*PSF(1).os;
y=((PSF(1).os+1)/2)+fw*PSF(1).os:PSF(1).os:PSF(1).Ny-((PSF(1).os-1)/2)-fw*PSF(1).os;
[X,Y]=ndgrid(x,y);
Z=ones(size(X));

if length(PSF)==2 %biplane imaging
    
   mode='biplane';
   ratio=1; %imbalance between the two channels can be considered here
   I_stack_a=I_stack{1}; 
   I_stack_b=I_stack{2};
   Fa=griddedInterpolant(PSF(1).data*PSF(1).os^2,'spline','nearest'); %multiplication with PSF.os^2 to preserve normalization 
   Fb=griddedInterpolant(PSF(2).data*PSF(2).os^2,'spline','nearest'); %multiplication with PSF.os^2 to preserve normalization 
   [nx,ny,no_images]=size(I_stack_a);
   
else %single-channel imaging
    
   mode='single';
   I_stack=data.images; %stack of single molecule images
   [nx,ny,no_images]=size(I_stack);
   Fa=griddedInterpolant(PSF(1).data*PSF(1).os^2,'spline','nearest'); %multiplication with PSF.os^2 to preserve normalization 

   %retrieving pixel-characteristics for every molecule image (readout noise, gain, baseline)
   var=zeros(nx,ny,no_images);
   gain=var; 
   offset=gain; 
   for m=1:no_images
       var(:,:,m)=cam.varmap(data.ULC(m,1)+(0:nx-1),data.ULC(m,2)+(0:ny-1));  %variance in ADUs
       gain(:,:,m)=cam.gainmap(data.ULC(m,1)+(0:nx-1),data.ULC(m,2)+(0:ny-1));  %gain must be in units of ADUs/electron
       offset(:,:,m)=cam.offsetmap(data.ULC(m,1)+(0:nx-1),data.ULC(m,2)+(0:ny-1)); %offset in ADUs
   end
   
end

options=optimset('fminsearch');
%options.MaxIter=100;

%------------------estimation loop------------------------
if strcmp(mode,'single')
    est=zeros(no_images,6); %initialization
elseif strcmp(mode,'biplane')
    est=zeros(no_images,8); %initialization
end

figure(1); 
for m=1:no_images
    
    %-------for single-channel imaging-------
    if strcmp(mode,'single')
    
        %sCMOS algorithm: see Huang et al. Nat.Meth. 2013
        I=(I_stack(:,:,m)-offset(:,:,m))/gain(:,:,m)+var(:,:,m)./gain(:,:,m).^2; %grab image from stack and modify according to paper of Huang et al.
        M_modif=var(:,:,m)./gain(:,:,m).^2; %model-modification for sCMOS detection 
        
        %define initial estimates and bounds
        BG0=min(I(:));
        N0=sum(I(:))-BG0*nx*ny;
        init_est=[0 0 z_ini_ind N0 BG0]; %initial estimates 
        LB=[-nx/3 -ny/3 1 0.25*N0 0]; %lower bounds
        UB=[+nx/3 +ny/3 PSF.Nz 4*N0 4*BG0]; %upper bounds

        tmp=fminsearchbnd(@fun_LLH,init_est,LB,UB,options); 
        x_est=tmp(1)*PSF.ux;
        y_est=tmp(2)*PSF.ux;
        z_est=(tmp(3)-1)*PSF.uz;
        N_est=tmp(4);
        BG_est=tmp(5);
        resnorm_MLE=sum((abs(M(:)-I(:)).^2))/(sum(I(:)))^2; 
        est(m,:)=[x_est.' y_est.' z_est.' N_est.' BG_est.' resnorm_MLE.'];

        %display every xth image
        if mod(m-1,50)==0 
            imagesc(I); 
            title([num2str(no_images-m) ' left']);
            hold on; 
            plot(round(nx/2)+y_est/PSF.ux/PSF.os,round(nx/2)+x_est/PSF.ux/PSF.os,'rx');
            hold off; 
            pause(0);
        end
        
    %-------for biplane imaging-------
    elseif strcmp(mode,'biplane')
    
        I_a=I_stack_a(:,:,m); %grab image from stack
        I_b=I_stack_b(:,:,m); %grab image from stack
        
        %define initial estimates and bounds
        BG0=min(I_a(:))+min(I_b(:)); 
        N0=sum(I_a(:)+I_b(:))-BG0*nx*ny;
        init_est=[0 0 z_ini_ind N0 BG0 0 0]; %initial estimates 
        LB=[-nx/3 -ny/3 1 0.25*N0 0 -nx/3 -ny/3]; %lower bounds
        UB=[+nx/3 +ny/3 PSF(1).Nz 4*N0 4*BG0 nx/3 ny/3]; %upper bounds
    
        tmp=fminsearchbnd(@fun_LLH_biplane,init_est,LB,UB,options); 
        x_est=tmp(1)*PSF(1).ux;
        y_est=tmp(2)*PSF(1).ux;
        z_est=(tmp(3)-1)*PSF(1).uz;
        N_est=tmp(4);
        BG_est=tmp(5);
        x2_est=tmp(6)*PSF(1).ux;
        y2_est=tmp(7)*PSF(1).ux;
        
        resnorm_MLE=sum(abs(Ma(:)-I_a(:)).^2+abs(Mb(:)-I_b(:)).^2)/(sum(I_a(:)+I_b(:)))^2; 
        est(m,:)=[x_est.' y_est.' z_est.' N_est.' BG_est.' resnorm_MLE.' x2_est.' y2_est.'];

        %display every xth image
        if mod(m-1,50)==0 
            imagesc([I_a, I_b]); 
            hold on; 
            plot(round(nx/2)+y_est/PSF(1).ux/PSF(1).os,round(nx/2)+x_est/PSF(1).ux/PSF(1).os,'rx');
            hold off;
            title([num2str(no_images-m) ' left']);
            pause(0);
        end
    end
end

%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
    
%defining log-likelihood function; works in conjunction with
%gradient-free minimum-search such as fminsearch
    function LLH=fun_LLH(v) 
        M=v(5)+v(4)*Fa(X-v(1),Y-v(2),Z*v(3)); %evaluating the model for parameter vector v
        M2=M+M_modif; 
        LLH=sum(sum(M2-I.*log(M2),1),2);
    end

    function LLH=fun_LLH_biplane(v)
        Ma=0.5*(v(5)+v(4)*Fa(X-v(1),Y-v(2),Z*v(3))); %evaluating the model for parameter vector v
        Mb=0.5/ratio*(v(5)+v(4)*Fb(X-v(6),Y-v(7),Z*v(3))); %evaluating the model for parameter vector v
        LLH=sum(sum(Ma-I_a.*log(Ma)+Mb-I_b.*log(Mb),1),2);
    end

end
