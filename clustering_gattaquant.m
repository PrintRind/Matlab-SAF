%clustering of localization data from gattaquant nanoruler sample
%needs data produced by "MLE_fit_molecules_exp.m"

%clear all; 
%load 'filtdata.mat' into workspace

%% DB clustering 
close all; 
clear X X2; 

%take junks of the data and process them individually
batchsize=10000; %how many locs can be processed by DBscan in one go?

filtdata_sorted=sortrows(filtdata,1); %sort data according to x-coordinate

%A) first look for clusters that are formed by single molecules to get rid of noise
%when should a point be assigned to a cluster:
epsilon=25; %max. allowed distance to neighbours % if the distance is chosen to be too small, a big cluster might get separated in multiple small ones 
MinPts=30; %how many molecules should be in the neighborhood?
showNR='n'; %show nanorulers during data evaluation? (slow)
 
%parameters for pairing clusters to nanorulers
epsilon_clu=120; %max. allowed distance to next cluster
MinPts_clu=2; %min. number of neighbors within that distance

datalength=length(filtdata); 
no_batches=ceil(datalength/batchsize);

clear clu_idx NR c_msd c_mean 
NR_count=0;

for b=1:no_batches
    
    disp(no_batches-b); 
    start=(b-1)*batchsize+1;
    stop=min(start+batchsize-1,datalength);
    X=[filtdata_sorted(start:stop,1) filtdata_sorted(start:stop,2)  filtdata_sorted(start:stop,3)];
    
    [clu_idx{b} ~]=DBSCAN(double(X),epsilon,MinPts); %returns the cluster number to which the coordinates belong to; idx=0 --> noise
    disp([num2str(max(clu_idx{b})) ' molecule clusters found.']);

    %extract mol.coords that belong to the same cluster; 
    %calculate mean positions and size of clusters and their distances to neighbor clusters 
    %clear NR c_mean c_msd NR_vec NR_length NR_angle
    for mm=1:max(clu_idx{b})
        X2{mm}=X(clu_idx{b}==mm,:); %extract mol.coords that belong to the same cluster 
        c_mean(mm,:)=mean(X2{mm},1); %mean position of cluster
        c_msd(mm)=mean(sum((double(X2{mm}(:,:))-c_mean(mm,:)).^2,2)); %mean squared distance from center-of-mass
        %X(clu_idx{b}==mm,:)
        %pause; 
        %scatter3(X2{m}(:,1),X2{m}(:,2),X2{m}(:,3));
        %axis equal; 
        %axis tight; 
    end
        
    NR_idx=DBSCAN(c_mean,epsilon_clu,MinPts_clu); %creating numbers for clusters
    NR_count=NR_count+max(NR_idx);
    disp([num2str(max(NR_idx)) ' nanorulers found.']);

    for m=1:max(NR_idx) %step through the nanorulers

        idx=find(NR_idx==m); %find cluster indices that form a nanoruler

        if length(idx)==2
            NR(b).vec(m,:)=c_mean(idx(2),:)-c_mean(idx(1),:);
            NR(b).length(m)=(norm(NR(b).vec(m,:)));
            NR(b).msd(m)=max(c_msd(idx(2)),c_msd(idx(1)));
            NR(b).angle(m)=acosd(abs(NR(b).vec(m,3))/NR(b).length(m));

            %showing nanorulers
            if strcmp(showNR,'y')
                for mm=1:length(idx)
                    scatter3(X2{idx(mm)}(:,1), X2{idx(mm)}(:,2), X2{idx(mm)}(:,3));
                    title(['NR length = ' num2str(NR(b).length(m),3) ' nm']);
                    hold on; 
                    axis equal; axis tight; 
                end
                pause; 
                hold off; 
            end
        else
            NR(b).vec(m,:)=[0 0 0];
            NR(b).length(m)=0;
            NR(b).msd(m)=0;
            NR(b).angle(m)=0;
        end

        %pause;
    end

    if exist('NR')
        %remove zeroes (represents NR with more or less than two clusters)
        zero_idx=find(NR(b).length==0);
        NR(b).vec(zero_idx,:)=[];
        NR(b).length(zero_idx)=[];
        NR(b).msd(zero_idx)=[];
        NR(b).angle(zero_idx)=[];

        %display
        figure(2); 
        plot(1:length(NR(b).length),NR(b).length,'b.');
        hold on; 
        plot(1:length(NR(b).length),sqrt(NR(b).msd),'r.');
        xlabel('no. of nanoruler');
        ylabel('nm');
        title(['NR length (blue) and width of larger cluster (RMSD); mean length=' num2str(mean(NR(b).length),3)]); 

        figure(3); 
        plot(NR(b).angle,NR(b).length,'bo');
        xlabel('NR angle / deg.');
        ylabel('NR length / nm');
        title('nanoruler length vs. angle to z-axis');
        xlim([0 90]);
        ylim([0 150]);
        hold on; 

        plot(0:90,90,'r.');
        plot(0:90,70,'r.');
    end
    
end

figure(2); hold off; 
figure(3); hold off; 

disp('done');
disp([num2str(NR_count) ' nanorulers found']);




