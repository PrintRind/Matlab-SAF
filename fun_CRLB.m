function [Cx,Cy,Cz,CN,Cbg,FI] = fun_CRLB(PSF,N,bg,varargin)
%calculation of Cramér-Rao lower bounds for a 3D-PSF (x,y,z)
%inputs: 
%-------
%PSF....for single channel imaging: PSF is member of PSF class; 
%       for dual channel imaging: PSF is cell array whose entries are
%       members of the PSF class, e.g.: PSF{1}=PSF_tot; PSF{2}=PSF_UAF;
%       The total energies in each x-y-slice of PSF.data (or alternatively PSF{1}.data and PSF{2}.data) should reflect the
%       relative energy content compared to the input parameter N, i.e. the
%       photons entering the objective lens
%N...   signal in no. of photons; Can also be a vector, in which case each entry in this vector corresponds to the photon number contained in the corresponding PSF-z-slice
%       N is the no of photons emitted by a molecule 
%       that go into the objective lens. The number of photons that finally end on
%       the molecule image on the camera is always at least slightly smaller, which should be considered by a corresponding normalization of the PSF-stack  
%bg.... background level in no. of photons per pixel
%       for dual-channel imaging, a vector [bg1, bg2] can be provided,
%       which defines separate background levels for each channel 
%varargs: 
%readnoise ... variance due to readout noise in electrons

%handling optional arguments
if length(varargin)==1
    readnoise = varargin{1}; 
else
    readnoise = 0;
end


%TWO CHANNEL IMAGING
if iscell(PSF) && size(PSF,2)==2  

     if numel(N)>1 %if N is a vector
         N = reshape(N,[1,1,1,Nz]); 
     end
     
     %creating a 4D-PSF with a dedicated "PSF-dimension" (=dim 1)
     [Nx, Ny, Nz] = size(PSF{1}.data); 
     PSFdual = zeros([2, size(PSF{1}.data)]);
     PSFdual(1,:,:,:) = PSF{1}.data; 
     PSFdual(2,:,:,:) = PSF{2}.data;
     
     %creating interpolant of PSF, which is also downsampled to the resolution of the camera          
     F=griddedInterpolant(N.*PSFdual.*PSF{1}.os^2,'spline','spline'); %extrapolation method: "spline" is important to ensure correct CRLB values at the borders 
     %F=griddedInterpolant(PSFdual./E.*PSF{1}.os^2,'spline','spline'); %extrapolation method: "spline" is important to ensure correct CRLB values at the borders 

     x=((PSF{1}.os+1)/2):PSF{1}.os:PSF{1}.Nx-((PSF{1}.os-1)/2);
     y=((PSF{1}.os+1)/2):PSF{1}.os:PSF{1}.Ny-((PSF{1}.os-1)/2);
     z=1:Nz;
     [M,X,Y,Z]=ndgrid(1:2,x,y,z); 
     PSF_ds=F(M,X,Y,Z); %downsampled 
     
        %calculating x,y,z gradients
        delta_x=0.01; 
        delta_z=0.01; 
        Dx=(F(M,X+delta_x,Y,Z)-F(M,X-delta_x,Y,Z))/(2*delta_x*PSF{1}.ux*1e9); 
        Dy=(F(M,X,Y+delta_x,Z)-F(M,X,Y-delta_x,Z))/(2*delta_x*PSF{1}.ux*1e9); 
        Dz=(F(M,X,Y,Z+delta_z)-F(M,X,Y,Z-delta_z))/(2*delta_z*PSF{1}.uz*1e9); 
        Dbg=1; 
        Dsig=PSF_ds./N; 
        
        if numel(bg)==2 %if two values for background are provided (separate for each channel)
            tmp = ones(size(PSF_ds)); 
            tmp(1,:,:,:)=bg(1); 
            tmp(2,:,:,:)=bg(2); 
            noise=(PSF_ds+tmp)+readnoise; %shot and read noise variances (eventl. dark noise if applicable)
        else
            noise=(PSF_ds+bg)+readnoise; %shot and read noise variances (eventl. dark noise if applicable)
        end
        
        %calculating Fisher information matrix entries
        FI_zz=squeeze(sum(sum(sum(Dz.^2./noise,1),2),3));%*1e-18; %in 1/nm^2
        FI_xx=squeeze(sum(sum(sum(Dx.^2./noise,1),2),3));%*1e-18; 
        FI_yy=squeeze(sum(sum(sum(Dy.^2./noise,1),2),3));%*1e-18; 
        FI_ss=squeeze(sum(sum(sum(Dsig.^2./noise,1),2),3)); %in 1/sig^2
        FI_bb=squeeze(sum(sum(sum(ones(size(Dsig))*Dbg^2./noise,1),2),3)); %in 1/sig^2

        FI_xy=squeeze(sum(sum(sum(Dx.*Dy./noise,1),2),3));%*1e-18;
        FI_xz=squeeze(sum(sum(sum(Dx.*Dz./noise,1),2),3));%*1e-18;
        FI_yz=squeeze(sum(sum(sum(Dy.*Dz./noise,1),2),3));%*1e-18;

        FI_xs=squeeze(sum(sum(sum(Dx.*Dsig./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_ys=squeeze(sum(sum(sum(Dy.*Dsig./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_zs=squeeze(sum(sum(sum(Dz.*Dsig./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_xb=squeeze(sum(sum(sum(Dx.*Dbg./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_yb=squeeze(sum(sum(sum(Dy.*Dbg./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_zb=squeeze(sum(sum(sum(Dz.*Dbg./noise,1),2),3));%*1e-9; %in 1/nm*1/sig
        FI_bs=squeeze(sum(sum(sum(Dsig.*Dbg./noise,1),2),3)); %in 1/sig^2
        
else %SINGLE CHANNEL IMAGING
    
    [Nx, Ny, Nz] = size(PSF.data);
    
    if numel(N)>1 %if N is a vector
        N = reshape(N,[1,1,Nz]); 
    end
    
    %creating interpolant of PSF
    F=griddedInterpolant(N.*PSF.data.*PSF.os^2,'spline','spline'); %extrapolation method: "spline" is important to ensure correct CRLB values at the borders; 
                                                                   %also compensating the loss of energy caused by the downsampling by multiplication with PSF.os^2 
    x=((PSF.os+1)/2):PSF.os:PSF.Nx-((PSF.os-1)/2);
    y=((PSF.os+1)/2):PSF.os:PSF.Ny-((PSF.os-1)/2);
    z=1:Nz;
    [X,Y,Z]=ndgrid(x,y,z); 
    PSF_ds=F(X,Y,Z); %downsampled 
  
    
    %calculating x,y,z gradients
    delta_x=0.01; 
    delta_z=0.01; 
    Dx=(F(X+delta_x,Y,Z)-F(X-delta_x,Y,Z))/(2*delta_x*PSF.ux*1e9); 
    Dy=(F(X,Y+delta_x,Z)-F(X,Y-delta_x,Z))/(2*delta_x*PSF.ux*1e9); 
    Dz=(F(X,Y,Z+delta_z)-F(X,Y,Z-delta_z))/(2*delta_z*PSF.uz*1e9); 
    Dbg=1; 
    Dsig=PSF_ds./N; 

    noise=(PSF_ds+bg)+readnoise; %shot and read noise variances (eventl. dark noise if applicable)
    
    %calculating Fisher information matrix entries
    FI_zz=squeeze(sum(sum(Dz.^2./noise,1),2));%*1e-18; %in 1/nm^2
    FI_xx=squeeze(sum(sum(Dx.^2./noise,1),2));%*1e-18; 
    FI_yy=squeeze(sum(sum(Dy.^2./noise,1),2));%*1e-18; 
    FI_ss=squeeze(sum(sum(Dsig.^2./noise,1),2)); %in 1/sig^2
    FI_bb=squeeze(sum(sum(ones(size(Dsig))*Dbg^2./noise,1),2)); %in 1/sig^2

    FI_xy=squeeze(sum(sum(Dx.*Dy./noise,1),2));%*1e-18;
    FI_xz=squeeze(sum(sum(Dx.*Dz./noise,1),2));%*1e-18;
    FI_yz=squeeze(sum(sum(Dy.*Dz./noise,1),2));%*1e-18;

    FI_xs=squeeze(sum(sum(Dx.*Dsig./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_ys=squeeze(sum(sum(Dy.*Dsig./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_zs=squeeze(sum(sum(Dz.*Dsig./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_xb=squeeze(sum(sum(Dx.*Dbg./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_yb=squeeze(sum(sum(Dy.*Dbg./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_zb=squeeze(sum(sum(Dz.*Dbg./noise,1),2));%*1e-9; %in 1/nm*1/sig
    FI_bs=squeeze(sum(sum(Dsig.*Dbg./noise,1),2)); %in 1/sig^2
    
end

    FI=zeros(5,5,Nz); %Fisher matrix

    %diagonal elements
    FI(1,1,:)=FI_xx;
    FI(2,2,:)=FI_yy;
    FI(3,3,:)=FI_zz;
    FI(4,4,:)=FI_ss;
    FI(5,5,:)=FI_bb;

    %off-diagonal elements
    FI(1,2,:)=FI_xy; FI(2,1,:)=FI_xy;
    FI(1,3,:)=FI_xz; FI(3,1,:)=FI_xz;
    FI(1,4,:)=FI_xs; FI(4,1,:)=FI_xs; 
    FI(1,5,:)=FI_xb; FI(5,1,:)=FI_xb;

    FI(2,3,:)=FI_yz; FI(3,2,:)=FI_yz;
    FI(2,4,:)=FI_ys; FI(4,2,:)=FI_ys;
    FI(2,5,:)=FI_yb; FI(5,2,:)=FI_yb;

    FI(3,4,:)=FI_zs; FI(4,3,:)=FI_zs;
    FI(3,5,:)=FI_zb; FI(5,3,:)=FI_zb;

    FI(4,5,:)=FI_bs; FI(5,4,:)=FI_bs;

    %calculating CRAMER RAO bounds
    Cx = zeros(1,Nz);
    Cy = Cx; 
    Cz = Cx; 
    CN = Cx; 
    Cbg = Cx;
    
    for m=1:Nz %inverting Fisher matrices and taking diagonal values
        FI_inv=inv(FI(:,:,m));
        Cx(m)=FI_inv(1,1);
        Cy(m)=FI_inv(2,2);
        Cz(m)=FI_inv(3,3);
        CN(m)=FI_inv(4,4);
        Cbg(m)=FI_inv(5,5);
    end

    
% %display
% z = (0:Nz-1)*uz*1e9; %z-axis in nm
% 
% figure; 
% plot(z,sqrt(Cx)); 
% hold on; 
% plot(z,sqrt(Cy)); 
% plot(z,sqrt(Cz)); 
% grid on; 
% hold off; 
% xlabel('z / nm');
% ylabel('\sigma / nm');
% legend('\sigma_x','\sigma_y','\sigma_z');

end


