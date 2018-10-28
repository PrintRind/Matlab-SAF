function CRB=fun_PSFopt(holo_in)
     
wxy=10;
wz=20;
n_photon=1000;

E3D=fun_focalfield_3DFFT_fast(holo_in,k_sphere,Nxy_pad,Nz_pad);
PSF=embed(abs(E3D).^2,[wxy,wxy,wz]);
[CRBx,CRBy,CRBz]=fun_CRB(PSF,ux,uz,n_photon,bg);
CRB=sum(sqrt(CRBx+CRBy+CRBz)); %defining fitness parameter for iterative optimization; the lower the better
