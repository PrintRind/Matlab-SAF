function hero=fun_findheroes(locdata,minlength,bridgelength,binsize)
%identify "heroes" amongst the blinking molecules, i.e. molecules that are
%"on" for several subsequent frames
%returned is a structure array containing information (corresponding
%locdatas) of potential "heroes"
%minlength.....minimum number of locs contained in the blink-series
%bridgelength....no. of empty frames in between two blink-events that
%should still be recognized as one molecule
%binsize...size of a binary pixel in nm (e.g. 50), defining the distance within which a
%location is seen as the same molecule
%for testing: 
%load sphereD_substack.mat; locdata=filtdata;
%locdata=[x_est' y_est' z_est' N_est' BG_est' resnorm_MLE' loc_no' frame_no' chi2'];
%minlength=2; %minimum length of a "hero"
%bridgelength=15; 

N=length(locdata); %no. of localizations
frame_no=locdata(:,8);
x_est=locdata(:,1);
y_est=locdata(:,2);

%---create data cube x,y,t----
%binsize=100; %cube increment side length in nm 
x_min=min(x_est);
x_max=max(x_est);
y_min=min(y_est);
y_max=max(y_est);
x=x_min:binsize:x_max;
y=y_min:binsize:y_max; 
nx=length(x);
ny=length(y);
frame_min=min(frame_no);
frame_max=max(frame_no);
frame_vec=unique(frame_no);
no_frames=numel(frame_vec);

cube=boolean(zeros(nx,ny,no_frames));

%draw "1" at positions within the cube, where molecules are located
for m=1:N
    %disp(N-m);
    x_pos=ceil((x_est(m)-x_min)/binsize); %molecule position in the cube
    y_pos=ceil((y_est(m)-y_min)/binsize); 
    if x_pos==0; x_pos=1; end
    if y_pos==0; y_pos=1; end
    %frame_pos=frame_no(m)-frame_min+1; 
    frame_pos=find(frame_vec==frame_no(m));
    cube(x_pos,y_pos,frame_pos)=1; 
    linind=sub2ind([nx,ny,no_frames],x_pos,y_pos,frame_pos);
    cubepos(m,:)=[locdata(m,:),linind]; %links filtdata with position (linear index) of the respective pixel within cube
end
disp('done');

%% quality check
% 
% m=14;
% imagesc(squeeze(cube(:,m,:))'); colormap gray; 
% title('spatio-temporal slice through cube');

%% convolving with a "line" in time to connect time-separated localizations frome the same molecule

if bridgelength>1 
    cube_filt = boolean(convn(cube,ones(1,1,bridgelength),'same'));
    %remove small objects (single, isolated localization events)
    cube_filt = bwareaopen(cube_filt,bridgelength+1);
else
    cube_filt= bwareaopen(cube,2); 
end

%% quality check

% m=14;
% imagesc(squeeze(cube_filt(:,m,:))'); colormap gray; 
% title(['filtered cube; section']);

%% find "heroes", i.e. connected components in time

CC = bwconncomp(cube_filt,18); %finding binary "islands"
numPixels = cellfun(@numel,CC.PixelIdxList); %retrieve size of those islands

cube2=boolean(zeros(size(cube))); 

%deleting islands containing too few localizations
for m=1:CC.NumObjects 
    if sum(cube(CC.PixelIdxList{m}))>=minlength %if more than the defined number of locs are contained in the "island"
        cube2(CC.PixelIdxList{m})=1; 
    end
end
CC2 = bwconncomp(cube2,18); %finding binary "islands"

%% quality check
% m=14;
% imagesc(squeeze(cube2(:,m,:))'); colormap gray; 
% title(['filtered cube; section']);

%% finding loc.-numbers that belong to the heroes

%stepping through all pixels that belong to heroes and finding & pooling the
%respective locdata entries in the structure array hero.
clear hero
for m=1:CC2.NumObjects
    tmp=find(cube(CC2.PixelIdxList{m})); %pixel indices of mth hero
    pixels=CC2.PixelIdxList{m}(tmp); 
    tmp=[];
    for mm=1:numel(pixels)
        %finding localizations that are assigned to pixel mm; this might be more than one localization 
        filtdata_tmp=cubepos(cubepos(:,end)==pixels(mm),1:end-1);
        tmp=[tmp; filtdata_tmp];
    end
    hero(m).locdata=tmp;
end

disp([num2str(numel(hero)) ' heroes found']);
disp('done');


