function xc=fun_find_sphere_center(xc0,r_coord,z_theory,filtdata)
%finds x-y-center of Alexa-coated sphere
%xc0...initial guess for center xc0=[x_ini y_ini]
%r_coord=linspace(0,30,size(PSF,3))*1e-6; %radial coord. in m
%z_theory=r_sphere-sqrt(r_sphere^2-(r_coord).^2);


z_error=@(x) std(interp1(r_coord,z_theory,sqrt((filtdata(:,1)-x(1)).^2+(filtdata(:,2)-x(2)).^2)*1e-9,'linear','extrap')*1e9-filtdata(:,3));
xc=fminsearch(z_error,xc0);