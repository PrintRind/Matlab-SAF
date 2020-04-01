function [I_out, FI]=fun_varyParam_calcImages(N,var_parameter,var_range,z_range,PSF,lens,cam,varargin)
%calculates single molecule images at different z-positions. 
%a variable parameter can be defined, such as aberrations or dipole
%emission characteristic
%INPUTS:
%-------
%N....size of grid (pupil diameter)
%var_parameter...param. that should be varied, e.g. Z5,Z11,.... or
%'z-dipole strength'
%var_range....vector of values for the variable parameter
%PSF....member of the psf-class
%lens...member of the objective class
%cam...member of the camera class

%-----handling optinal input arguments-----

nVarargs=length(varargin); 

if nVarargs==3
    noise=varargin{1};
    sig=varargin{2};
    bg=varargin{3};
    mode='fullpupil';
elseif nVarargs==4
    noise=varargin{1};
    sig=varargin{2};
    bg=varargin{3};
    mode=varargin{4};  
else
    noise='n';
    sig=1; 
    bg=0;
    mode='fullpupil';
end

%-------------------------------------------

no_vars=length(var_range);
T=lens.transmission(N); 

%--------
%ux=cam.pixsize/lens.M*lens.f_tube/f_tube; %effective pixel size in focal space
ux=PSF.ux*PSF.os;
Nx=PSF.Nx/PSF.os;
[~,Defocus,~] =fun_SA_RImismatch(N,lens.RI,lens.RI,lens.NA,PSF.lambda,1); %Defocus function refers to refractive index n2
uk=2*pi/PSF.lambda*2*lens.NA/N; %unit in pupil space (k-space)
[~,~,R,pupil]=create_coord(N,1,'FFT');
d_layer=0e-9; %thickness of intermediate layer (layer 2)
mu_x=1e-10; %magnitude of x-dipole (arbitrary value)
mu_y=1e-10; %magnitude of y-dipole (arbitrary value)
mu_z=1e-10; %magnitude of z-dipole (arbitrary value)
%---------

if strcmp(mode,'UAF')
    pupil=circshift(R<=((N/2)*(PSF.RI/lens.NA))*1,[0 0]); %pupil containing UAF light
end

if strcmp(var_parameter,'z-dipole strength')
    z_strength=var_range;
    a_Z=linspace(0,0,no_vars);
    m_Z=1; %Noll no. of mode
    legend_param=num2str(z_strength);
elseif strcmp(var_parameter(1),'Z')
    m_Z=str2num(var_parameter(2:end)); %Noll no. of mode
    a_Z=var_range; %magnitude of varying aberr. mode
    z_strength=linspace(1,1,no_vars);
    legend_param=num2str(a_Z);
end

Z=ZernikeCalc(m_Z,1,pupil,'NOLL');

no_vars=max(length(z_strength),length(a_Z)); %number of parameter variations
I_test=zeros(Nx,Nx,length(z_range),no_vars);

%include aberrations that are contained in the PSF for estimation
if isempty(PSF.Zernike)==0
    aberr=sum(ZernikeCalc(PSF.Zernike(1,:),PSF.Zernike(2,:)',pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded 
else
    aberr=0;
end



%--------------------------------calculate image stack------------------------------------
clear I_tmp;
for m=1:length(z_range)
    %fun_dipole_imaging(N,PSF.lambda,NA,RI,[0,0],d2,z_vec(m),f,mu_z,T); %z-dipole
    disp(length(z_range)-m);
    [Ex_Pz(:,:,m),Ey_Pz(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,[PSF.RI PSF.RI lens.RI],[0,0],d_layer,z_range(m),lens.f,mu_z,T); %z-dipole
    [Ex_Px(:,:,m),Ey_Px(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,[PSF.RI PSF.RI lens.RI],[pi/2, 0],d_layer,z_range(m),lens.f,mu_x,T); %x-dipole
    [Ex_Py(:,:,m),Ey_Py(:,:,m)]=fun_dipole_imaging(N,PSF.lambda,lens.NA,[PSF.RI PSF.RI lens.RI],[pi/2, pi/2],d_layer,z_range(m),lens.f,mu_y,T); %y-dipole
    
    for v=1:no_vars 
        
        mask=pupil.*exp(1i*aberr+1i*PSF.defocus.*Defocus+1i*a_Z(v)*Z); 
                
        I_BFP=(abs(Ex_Px(:,:,m)).^2+abs(Ex_Py(:,:,m)).^2+abs(Ey_Px(:,:,m)).^2+abs(Ey_Py(:,:,m)).^2+abs(z_strength(v)*Ex_Pz(:,:,m)).^2+abs(z_strength(v)*Ey_Pz(:,:,m)).^2);
        E_tot(m,v)=sum(sum(I_BFP));
                   
        %-----calculating total (SAF+UAF) images-----
        I_xx=abs(czt2(Ex_Px(:,:,m).*mask,uk,ux,Nx)).^2;
        I_yx=abs(czt2(Ey_Px(:,:,m).*mask,uk,ux,Nx)).^2;
        I_xy=abs(czt2(Ex_Py(:,:,m).*mask,uk,ux,Nx)).^2;
        I_yy=abs(czt2(Ey_Py(:,:,m).*mask,uk,ux,Nx)).^2;    
        I_xz=z_strength(v)*abs(czt2(Ex_Pz(:,:,m).*mask,uk,ux,Nx)).^2;
        I_yz=z_strength(v)*abs(czt2(Ey_Pz(:,:,m).*mask,uk,ux,Nx)).^2;
        I_tmp(:,:,m,v)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m,v);
%       imagesc(squeeze(I_test(:,:,m,mm))); 
%       pause(.1); 
    end
end
disp('done'); 


for v=1:no_vars
    %----calculate Fisher information of current image stack----
    tmp=PSF; %copy PSF
    tmp.data=I_tmp(:,:,:,v);
    tmp.ux=tmp.ux*tmp.os; 
    tmp.os=1;
    tmp.Nx=size(I_tmp,1);
    tmp.Ny=size(I_tmp,2);
    [~,~,~,~,~,FI{v}]=tmp.CRLB(sig,bg,cam);
end

%adding noise if selected
if strcmp(noise,'y')
    I_out=round(poissrnd(cam.QE*(sig*I_tmp+bg))/cam.amp+cam.baseline);
else
    I_out=round(cam.QE*(sig*I_tmp+bg)/cam.amp+cam.baseline);
end
