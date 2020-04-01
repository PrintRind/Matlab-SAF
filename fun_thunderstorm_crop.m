function [mol,ULC]=fun_thunderstorm_crop(TS,data_stack,data_BG,w,ux)
%crops molecules from image data based on ThunderStorm localization data
%w...half width of molecule images
%ux...size of a pixel
%data_BG...background-image, e.g. obtained using EVER
%outputs: 
%-) mol... molecule image
%-) ULC...row and column position of the uppper left corner of the molecule
%image 
img_no=size(TS,1); %number of detected molecules
mol=zeros(2*w+1,2*w+1,img_no);

frame=TS(:,2); %frame number
x_pos=TS(:,3)/(ux*1e9)+0; %x-position, convert into units of pixels
y_pos=TS(:,4)/(ux*1e9)+0;

for m=1:img_no
    rows=round((y_pos(m)-w):(y_pos(m)+w));
    cols=round((x_pos(m)-w):(x_pos(m)+w));
    if (min(rows)>0) && (max(rows)<size(data_stack,1)) && (min(cols)>0) && (max(cols)<size(data_stack,2))
        img=data_stack(rows,cols,frame(m));
        bg=data_BG(rows,cols);
        mol(:,:,m)=img+mean(bg(:)); %adding mean background
        %imagesc(mol(:,:,m)); pause(0);
    end  
    ULC(m,:)=[rows(1), cols(1)]; %upper left corner (ULC) coordinates of images
end


    