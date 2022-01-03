fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

binSize = 200;
dsRate = 5;
trlWindow = {[-800 1500]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
lfpWindow = [16 32];
numPerms = 100;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

frThresh = 0.1/(binSize/1000);
% frThresh = 0;
%% Data Vectors
popVects = cell(length(fileDirs),1);
popVectsSortVect = cell(length(fileDirs),1);
popVectsThreshVect = cell(length(fileDirs),1);
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    [trlSpikes, timeVect] = mlb.PP_TrialMatrix_Spiking(trlWindow{1}, alignment{1});
    nsmblSpks = nan(size(trlSpikes,2),size(trlSpikes,1), mlb.seqLength);
    peakNdx = nan(size(trlSpikes,2),mlb.seqLength);
    peakRt = nan(size(trlSpikes,2),mlb.seqLength);
    for op = 1:mlb.seqLength
        tempNsmbl = mean(trlSpikes(:,:,mlb.fiscTrials(op,:)),3);
        for uni = 1:size(tempNsmbl,2)
            peakNdx(uni,op) = find(tempNsmbl(:,uni)==max(tempNsmbl(:,uni)),1,'first');
            peakRt(uni,op) = max(tempNsmbl(:,uni));
            tempNsmbl(:,uni) = tempNsmbl(:,uni)./max(tempNsmbl(:,uni));
        end
        nsmblSpks(:,:,op) = tempNsmbl';
    end
    popVects{ani} = nsmblSpks;
    popVectsSortVect{ani} = peakNdx;
    popVectsThreshVect{ani} = peakRt;
end
mlb.SetLikes_FISC;
grpPV = cell2mat(popVects);
grpPVsort = cell2mat(popVectsSortVect);
grpPVthresh = cell2mat(popVectsThreshVect);
%% Plot PopVects
figure;
ndxCorr = nan(mlb.seqLength);
pvCorr = nan(mlb.seqLength);
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        sortedPV = sortrows([grpPVsort(:,op1), grpPVthresh(:,op2), grpPV(:,:,op2)]);
        subplot(mlb.seqLength,mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
        imagesc(mlb.likeTimeVect, 1:size(sortedPV,1), sortedPV(sortedPV(:,2)>=frThresh,3:end), [0 1]);
        ndxCorr(op1,op2) = corr(grpPVsort(:,op1), grpPVsort(:,op2), 'rows', 'pairwise');
        tempX = grpPV(grpPVthresh(:,op2)>=frThresh,:,op1);
        tempY = grpPV(grpPVthresh(:,op2)>=frThresh,:,op2);
        pvCorr(op1,op2) = corr(tempX(:), tempY(:), 'Rows', 'pairwise');
    end
end
figure; 
subplot(1,2,1);
imagesc(ndxCorr, [0 1]);
subplot(1,2,2);
imagesc(pvCorr, [0 1]);