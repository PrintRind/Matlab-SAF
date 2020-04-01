function [XX_Corr,YY_Corr]=fun_DriftCorrection(XX,YY,ZZ,TS,BinSize,nCorrStep,z)
    % xy-drift correctio
    % XX = list of x-positions
    % XX = list of y-positions
    % TS = list of corresponding image numbers
    % BinSize
    % # of correction steps
    % z switches z correction on

    nIMG = max(TS);

    steps = floor(linspace(1,nIMG,nCorrStep+1));

    LbX=min(XX(:));
    LbY=min(YY(:));
    UbX=max(XX(:));
    UbY=max(YY(:));

    BinSizeX = BinSize;
    BinSizeY = BinSize;
    BinNrX = floor((UbX-LbX)/BinSizeX);
    BinNrY = floor((UbY-LbY)/BinSizeY);
    
    % Select by imagenumber range
    for i=1:nCorrStep
        Select(i,:) = (TS(:)>steps(i))& (TS(:)<steps(i+1));
    end

    
    % Calculate histograms
    for i=1:nCorrStep
        H_(i,:,:)=(hist2d([XX(Select(i,:)),YY(Select(i,:))],linspace(LbX,UbX,BinNrX),linspace(LbY,UbY,BinNrY)));
    end
    H=imgaussfilt(H_,1);

    % Reference
    bla=fftshift(ifftn(fftn(fftshift(H(1,:,:))).*conj(fftn(fftshift(H(1,:,:)))),'symmetric'));
    [Ra,Rb,Rc] = ind2sub(size(bla),find(bla == max(bla(:))));

    % Calculate correction
    for i=1:nCorrStep
        bla=fftshift(ifftn(fftn(fftshift(H(1,:,:))).*conj(fftn(fftshift(H(i,:,:)))),'symmetric'));
        [a,b,c] = ind2sub(size(bla),find(bla == max(bla(:))));
        CorrFac_(i,:)=[a,b,c];
    end

    % Correction factor
    CorrFac = (CorrFac_-[Ra,Rb,Rc])*BinSizeX;
    
    % spline interp
    % X-dir
       XspX = steps;
       YspX = [0 CorrFac(:,2)'];
       XXspX = floor(linspace(1,nIMG,nIMG));
       YYspX = spline(XspX,YspX,XXspX);
    % Y-dir
       XspY = steps;
       YspY = [0 CorrFac(:,3)'];
       XXspY = floor(linspace(1,nIMG,nIMG));
       YYspY = spline(XspY,YspY,XXspY);
       
       
    % Correct data
    XX_Corr=XX;
    YY_Corr=YY;

    XX_Corr = XX_Corr+YYspX(TS(:))';
    YY_Corr = YY_Corr+YYspY(TS(:))';

    figure(33)
    scatter(XX_Corr ,YY_Corr,5,1:length(XX_Corr),'filled')
    title('Corrected data')

    figure(34)
    scatter(XX ,YY,5,1:length(XX),'filled')
    title('uncorrected data')

    figure(35)
    plot(XspX,YspX,'bo',XXspX,YYspX,'b-.',XspY,YspY,'mo',XXspY,YYspY,'m--')
    legend('X-drift','Y-drift')
    title('relative drift XY')
    
    
    
        % z-correction
    if z
        nCorrStepZ=3;

        % slice histogramms in z in z
        LdataZ = length(ZZ);
        stepsZ = floor(linspace(1,LdataZ,nCorrStep+1));

        numBunch_ = diff(stepsZ);
        numBunch = numBunch_(1)+10;

        BinNrZ = 25;
        SumSize = floor((UbY-LbY)/BinNrZ);

        % generate x projected slices

        for j =1:nCorrStepZ
            Hz=[];
            for i =1:BinNrZ 
                clear SelectTMP
                SelectTMP = (abs(YY_Corr(Select(j,:),:)-SumSize/2-SumSize*(i-1))<SumSize/2);
                Hz_(:,:)=(hist2d([XX_Corr(SelectTMP),ZZ(SelectTMP)],linspace(LbX,UbX,BinNrX),linspace(0,250,250/BinNrZ)));
                Hz = [Hz ,Hz_'];
            end
            Hz=imgaussfilt(Hz,1);
            HZ(j,:,:)=Hz;
        end

        %
        clear CorrFacZ_ CorrFacZ

        % Reference
        bla=fftshift(ifftn(fftn(fftshift(HZ(1,:,:))).*conj(fftn(fftshift(HZ(1,:,:)))),'symmetric'));
        [Ra,Rb,Rc] = ind2sub(size(bla),find(bla == max(bla(:))));

        % Calculate correction
        for i=1:nCorrStepZ
            bla=fftshift(ifftn(fftn(fftshift(HZ(1,:,:))).*conj(fftn(fftshift(HZ(i,:,:)))),'symmetric'));
            [a,b,c] = ind2sub(size(bla),find(bla == max(bla(:))));
            CorrFacZ_(i,:)=[a,b,c];
        end

        % Correction factor
        CorrFacZ = (CorrFacZ_-[Ra,Rb,Rc])*BinNrZ


        for i=1:nCorrStepZ
            ZZ_Corr(Select(i))=ZZ(stepsZ(Select(i)))+CorrFacZ(i,3);
        end
    end
end
