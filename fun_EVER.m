function [data_corr, data_BG]=fun_EVER(data)
%EVER background correction 
%input: 
%data_stack....stack of SMLM images

V_min=min(data,[],3);
V_BG1=(V_min+27.8)/0.9513 ; 
imagesc(V_BG1); title('background');
data_corr=data-repmat(V_BG1,1,1,size(data,3));
data_BG=V_BG1;
disp('done; EVER BG-corr. used');
