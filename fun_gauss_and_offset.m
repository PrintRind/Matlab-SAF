function F = fun_gauss_and_offset(x,xdata)
 %[offset Amp x-shift y-shift sigma]
 F = x(1)+x(2).*exp(-((xdata(:,:,1)-x(3)).^2/(2*x(5)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ) );
 