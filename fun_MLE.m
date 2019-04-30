function est=fun_MLE(I,PSF_norm,interp_incr,x_data,z_ini_info,resnorm_thresh,showimage)
%function provides position, signal and background estimate of
%single-molecule images using MLE approach
%inputs: 
%I...small image containing the molecule approximately in the center
%PSF_norm...normalized PSF model (energy normalized to 1)
%interp_incr...interpolation increment (in pixels) used to interpolate the
%PSF for x- and y-shifts 
%os...spatial oversampling of the model-PSF with respect to the camera
%pixel size (e.g. os=3); the higher the oversampling the smaller are the
%interpolation artifacts
%x_data...coordinate system used for initial Gaussfit
%z_ini...initial z-estimate (number of corresponding frame in PSF-stack)
%resnorm_thresh...quality threshold --> fits with too large resnorm are
%discarded
%output: 
%[x,y,z,signal,BG]

[nx, ny]=size(I); %size of molecule image

%defining Log-likelihood function
if ndims(PSF_norm)==3 %3D-PSF: requires interpolating in the function "fun_LLH" -> slow
    [nx0,~,nz0]=size(PSF_norm);
    fw=(nx0-nx)/2; %frame width
    fun_LLH=@(v) sum(sum(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))-I.*log(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3),'nearest')),1),2);
elseif ndims(PSF_norm)==5 %5D-PSF: no interpolation required
    [nx0,~,nz0,nxi,nyi]=size(PSF_norm);
    fw=(nx0-nx)/2; %frame width
    %fun_LLH=@(v) sum(sum(v(5)+v(4)*squeeze(PSF_norm((fw+1:nx+fw),(fw+1:nx+fw),uint8(v(3))+1,uint8(v(2))+1,uint8(v(1))+1))-I.*log(v(5)+v(4)*squeeze(PSF_norm((fw+1:nx+fw),(fw+1:nx+fw),uint8(v(3))+1,uint8(v(2))+1,uint8(v(1))+1))),1),2); %note that row and column indices are swapped
    fun_LLH=@(v) sum(sum(v(5)+v(4)*squeeze(PSF_norm(1+fw:nx0-fw,1+fw:nx0-fw,uint8(v(3))+1,uint8(v(2))+1,uint8(v(1))+1))-I.*log(v(5)+v(4)*squeeze(PSF_norm(1+fw:nx0-fw,1+fw:nx0-fw,uint8(v(3))+1,uint8(v(2))+1,uint8(v(1))+1))),1),2); %note that row and column indices are swapped
    cx=ceil((nxi+1)/2); %center pixel x-dir.
    cy=ceil((nyi+1)/2); %center pixel y-dir.
end

%defining binary mask for initial background estimation
BGmask1=zeros(nx,ny);
BGmask1(1,:)=1; BGmask1(end,:)=1; BGmask1(:,1)=1; BGmask1(:,end)=1;

BG00=sum(BGmask1(:).*I(:))/sum(BGmask1(:)); %initial background est.
       
        %----Gaussfit for initial estimates of x-y-position,BG,N0-----
        param_ini=[BG00 max(I(:))-BG00 0 0 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
        lb=[0 0.5*(max(I(:))-BG00) -nx/3 -nx/3 0.5]; %lower bounds for fit-parameters
        ub=[10*BG00 1.5*(max(I(:))-BG00) nx/3 nx/3 3]; %upper bounds
        [Param,~,~,~]=lsqcurvefit(@fun_gauss_and_offset_test,param_ini,x_data,I,lb,ub);
        %Gaussfit=fun_gauss_and_offset_test(Param,x_data);
        x0=Param(3); %in pixels (vertical direction in image I)
        y0=Param(4); %in pixels
        BG0=Param(1); %this serves as initial estimate for the MLE fit
        %N0=Param(2)*2*pi*Param(5)^2; %energy contained (see e.g. formula in Thunderstorm script)
        N0=sum(I(:))-BG0*nx*ny; %initial guess for number of photons - seems to work better 
        if isscalar(z_ini_info) %info about initial z-estimate
            z_ini=z_ini_info;
        else  %if provided z_ini_info is not scalar, then is describes a z_ini versus Gaussian-width-curve which provides a good initial z-estiamte based on the Gaussian fit width
            z_ini=interp1(z_ini_info,1:nz0,Param(5),'linear','extrap');
        end
        
            
        
        %----MLE estimation----------    
        if ndims(PSF_norm)==3
            gauss_est=[x0 y0 z_ini N0 BG0]; %initial estimates (e.g. from Gaussfit)
            LB=[-nx/3 -ny/3 0 0.25*N0 0.25*BG0]; %lower bounds
            UB=[+nx/3 +ny/3 nz0 4*N0 4*BG0]; %upper bounds
            tmp=fminsearchbnd(fun_LLH,gauss_est,LB,UB);
            x_est=tmp(1);
            y_est=tmp(2);
            z_est=tmp(3);
            N_est=tmp(4);
            BG_est=tmp(5);

        elseif ndims(PSF_norm)==5
            %-----------MLE-------------
            gauss_est=[cx+x0/interp_incr cy+y0/interp_incr z_ini-1 N0 BG0]; %initial estimates (e.g. from Gaussfit)
            %gauss_est=[x0/(interp_incr/os) y0/(interp_incr/os) z_ini N0 BG0]; %initial estimates (e.g. from Gaussfit)

            LB=[0 0 0 0.25*N0 0.25*BG0]; %lower bounds
            UB=[nxi-1 nyi-1 nz0-1 4*N0 4*BG0]; %upper bounds
            tmp=fminsearchbnd(fun_LLH,gauss_est,LB,UB);
            x_est=-(cx-(tmp(1)+1))*(interp_incr); %converting back to (camera) pixels
            y_est=-(cy-(tmp(2)+1))*(interp_incr);
            z_est=tmp(3)+1;
            N_est=tmp(4);
            BG_est=tmp(5);
          
            I_model=tmp(5)+tmp(4)*PSF_norm((fw+1:nx0-fw),(fw+1:nx0-fw),uint8(tmp(3))+1,uint8(tmp(2))+1,uint8(tmp(1))+1);
            resnorm_MLE=sum(abs(I_model(:)-I(:)).^2)/(BG_est+N_est)^2; 
        end

        est=[x_est y_est z_est N_est BG_est resnorm_MLE];
        
%         disp('initial and final x-y-z-estimates:'); disp([x0,x_est,y0,y_est,z_ini-1,z_est]);
%         pause; 
%         %disp(resnorm/N_est);
          
        %discarding erroneous localizations
        if (resnorm_MLE>resnorm_thresh) || (abs(x_est)>0.98) || (abs(y_est)>0.98) %if fitting error too large --> discard
            est=[];
            disp('localization discarded');
            disp(resnorm_MLE);
            disp(x_est);
            disp(y_est);
            
        else
            if showimage
                imagesc(I); axis equal; axis tight; title(['sig./ bg =' num2str(N_est) ' / ' num2str(BG_est)]); 
                pause(0);
            end
        end
        