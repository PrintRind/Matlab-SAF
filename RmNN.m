function TS_sorted = RmNN(NrIMG,TS_,Treshold,Sigma,Uncertainty)
    % NrIMG is ther number of images in stack to be analyzed
    % TS_ = Thunderstorm file
    % Treshold = treshold in nm
    % Sigma = with of gaussian (~500nm)
    % Uncertainty of localisation (eg. < 25nm)
    % anyone who finds typos can keep them ;)

    ImgNr=NrIMG;
    clear Select
    Select = find(TS_(:,2)<=ImgNr);
    TS=TS_(Select,:);
    clear Select

    % use cosine similarity to fliter close neighbores
    clear steps
    clc

    TS_sorted=[];
    Delete_=[];

    for j=1:max(TS(:,2))
        Delete_=[];
        % steps gives starting point
        steps=1;
        doIt=0;

        % select single image
        clear Select TS_new
        Select = find(TS(:,2)==j);
        TS_new=TS(Select,:);
       if not(isempty(TS_new))
        clear Select

        % runs through all localisations
        while doIt==0
           Switch = 0;
           LenTS = length(TS_new(:,2));

           % do while steps is smaler than LenTS
           if LenTS >= steps
               Delete=[];

               % initialize starting vector
               V=[TS_new(steps,3),TS_new(steps,4)];
               % list of number
               elements=[1:LenTS];
               % Exclude current element
               elements = elements(elements~=steps);

               % calculate distance to every other vector 
               for i = elements
                  V2=[TS_new(i,3),TS_new(i,4)];
                  distance = norm(V2-V);

                  % mark position if smaller than treshold
                  if distance<Treshold | TS_new(i,5)> Sigma | TS_new(i,10)> Uncertainty
                        Delete_=[Delete_ i];
                        Switch = 1;
                  end      
               end

               Delete = [Delete Delete_];
               steps = steps+1;
           else
               doIt=1;
           end

        end
        % Delete close neighbores
        TS_new(unique(Delete_','rows'),:)=[];
        TS_sorted=[TS_sorted; TS_new];
       end
    end

    [before, ~]=size(TS);
    [after, ~]=size(TS_sorted);
    
    disp(sprintf('%1.f%% of the data sorted out',100-100*after/before))