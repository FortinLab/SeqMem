tets = [];
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% setupSeqLength = 4; 

% % CA1 Data
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% % tets = [1,22,17,18,17]; % Lateral/Proximal
% tets = [7,3,1,5,5]; % Medial/Distal
% setupSeqLength = 5;

% All well-trained files PFC
fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
setupSeqLength = 4; 

% 
% % Only animals that show "normal" beta PFC
% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
% setupSeqLength = 4; 
binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
numPerms = 10;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
% Behavior Variables
fiscPokeOutLat = cell(length(fileDirs),1);
fiscRwdDelivLat = cell(length(fileDirs),1);
smi = nan(length(fileDirs),1);
dPrm = nan(length(fileDirs),1);
ri = nan(length(fileDirs),1);
smiByOP = nan(length(fileDirs),setupSeqLength,2);
dPrmByOP = nan(length(fileDirs),setupSeqLength,2);
riByOP = nan(length(fileDirs),setupSeqLength,2);
% Neural Variables
aniPosts = cell(size(fileDirs));
aniOdorDecodes = cell(size(fileDirs));
aniTimeDecodes = cell(size(fileDirs));
aniPosAccLog = cell(size(fileDirs));
aniTimeAccLog = cell(size(fileDirs));
aniTrlPosIDs = cell(size(fileDirs));
aniLFP = cell(size(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
% %     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(tets)
        mlb.lfpRefTet = find(strcmp(mlb.lfpMatrixColIDs, sprintf('T%i', tets(ani))));
    end
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.numPerms = numPerms;
    mlb.ssProportion = ssProportion;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    %% Extract Behavioral Variables
    if strcmp(alignment{1}, 'PokeIn')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
    elseif strcmp(alignment{1}, 'PokeOut')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeInIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
    end
    smi(ani) = mlb.smi;
    dPrm(ani) = mlb.dPrime;
    ri(ani) = mlb.ri;
    for op = 1:mlb.seqLength
        smiByOP(ani,:,1) = mlb.smiByPos;
        smiByOP(ani,:,2) = mlb.smiByOdr;
        dPrmByOP(ani,:,1) = mlb.dPrimeByPos;
        dPrmByOP(ani,:,2) = mlb.dPrimeByOdr;
        riByOP(ani,:,1) = mlb.riByPos;
        riByOP(ani,:,2) = mlb.riByOdr;
    end
    %% Extract LFP
    [betaPhase, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
    [thetaPhase, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
    %% Decode ISC via Sub-Sampling
%     mlb.bayesType = 1;  % Comment In to decode using Gaussian rather than Poisson
    mlb.SetLikes_SubSample;
%     mlb.bayesType = 3;  % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_Observes;
    tempPosts = cell2mat(reshape(mlb.post, [1,1,numel(mlb.post)]));
    postTrlIDs = cell2mat(reshape(mlb.postTrlIDs, [1,1,numel(mlb.postTrlIDs)]));
    tempOdorDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,3));
    tempPosAccuracy = arrayfun(@(a,b)a==b,repmat(tempOdorDecode,[1,1,mlb.seqLength]), repmat(permute(1:mlb.seqLength,[3,1,2]), size(tempOdorDecode,1),size(tempOdorDecode,2),1));
    tempTimeDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    tempTimeAccuracy = arrayfun(@(a,b)a==b,repmat(tempTimeDecode,[1,1,length(timePoints)]), repmat(permute(timePoints,[3,2,1]), size(tempTimeDecode,1),size(tempTimeDecode,2),1));
    tempLFP = permute([betaPower(:,:,postTrlIDs), thetaPower(:,:,postTrlIDs)], [1,3,2]);
    %% Outputs
    aniPosts{ani} = cell2mat(reshape(mlb.post, [1,1,mlb.numPerms]));
    aniOdorDecodes{ani} = tempOdorDecode;
    aniTimeDecodes{ani} = tempTimeDecode;
    aniPosAccLog{ani} = tempPosAccuracy;
    aniTimeAccLog{ani} = tempTimeAccuracy;
    aniTrlPosIDs{ani} = [mlb.trialInfo(postTrlIDs).Position];
    aniLFP{ani} = tempLFP;
end
pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

%%
if strcmp(alignment{1}, 'PokeIn')
    trlTimeLog = mlb.obsvTimeVect>0;
else
    trlTimeLog = mlb.obsvTimeVect>-0.5;
end
grpTrlPosts = cell2mat(reshape(aniPosts, [1,1,length(fileDirs)]));
grpTrlLFP = cell2mat(aniLFP);
grpTrlPos = cell2mat(aniTrlPosIDs);
grpTrlPosAccLog = cell2mat(aniPosAccLog);
grpTrlTimeAccLog = cell2mat(aniTimeAccLog);
grpTrlDecodes = cell2mat(aniOdorDecodes);
betaTrlMean = mean(grpTrlLFP(trlTimeLog,:,1));
thetaTrlMean = mean(grpTrlLFP(trlTimeLog ,:,2));

lfpThresh = nan(mlb.seqLength,2,2);
threshPosts = cell(mlb.seqLength,2,2);
threshAccPos = cell(mlb.seqLength,2,2);
threshAccTime = cell(mlb.seqLength,2,2);

figure;
for p = 1:mlb.seqLength
    subplot(mlb.seqLength,2, sub2ind([2,mlb.seqLength], 1,p));
%     betaThresh = [mean(betaTrlMean(grpTrlPos==p))-(std(betaTrlMean(grpTrlPos==p))*1), mean(betaTrlMean(grpTrlPos==p))+(std(betaTrlMean(grpTrlPos==p))*1)]; % INDIVIDUAL threshold (Mean+/-STD)
%     betaThresh = [mean(betaTrlMean(grpTrlPos))-(std(betaTrlMean(grpTrlPos))*1), mean(betaTrlMean(grpTrlPos))+(std(betaTrlMean(grpTrlPos))*1)]; % SESSION threshold (Mean +/-STD)
    sortedBeta = sort(betaTrlMean(grpTrlPos==p)); betaThresh = [sortedBeta(ceil(length(sortedBeta)*0.25)); sortedBeta(floor(length(sortedBeta)*0.75))]; % INDIVIDUAL threshold (<25/>75)
%     sortedBeta = sort(betaTrlMean(grpTrlPos)); betaThresh = [sortedBeta(ceil(length(sortedBeta)*0.25)); sortedBeta(floor(length(sortedBeta)*0.75))]; % SESSION threshold (<25/>75)
    lfpThresh(p,:,1) = betaThresh;
    histogram(betaTrlMean(grpTrlPos==p));    
    title(sprintf('Beta %i', p));
    threshPosts{p,1,1} = grpTrlPosts(:,:,grpTrlPos==p & betaTrlMean>betaThresh(2));
    threshPosts{p,2,1} = grpTrlPosts(:,:,grpTrlPos==p & betaTrlMean<betaThresh(1));
        
    subplot(mlb.seqLength,2, sub2ind([2,mlb.seqLength], 2,p));
%     thetaThresh = [mean(thetaTrlMean(grpTrlPos==p))-(std(thetaTrlMean(grpTrlPos==p))*1), mean(thetaTrlMean(grpTrlPos==p))+(std(thetaTrlMean(grpTrlPos==p))*1)]; % INDIVIDUAL threshold (Mean+/-STD)
%     thetaThresh = [mean(thetaTrlMean(grpTrlPos))-(std(thetaTrlMean(grpTrlPos))*1), mean(thetaTrlMean(grpTrlPos))+(std(thetaTrlMean(grpTrlPos))*1)]; % SESSION threshold (Mean +/-STD)
    sortedTheta = sort(thetaTrlMean(grpTrlPos==p)); thetaThresh = [sortedTheta(ceil(length(sortedTheta)*0.25)); sortedTheta(floor(length(sortedTheta)*0.75))]; % INDIVIDUAL threshold (<25/>75)
%     sortedTheta = sort(thetaTrlMean(grpTrlPos)); thetaThresh = [sortedTheta(ceil(length(sortedTheta)*0.25)); sortedTheta(floor(length(sortedTheta)*0.75))]; % SESSION threshold (<25/>75)
    lfpThresh(p,:,2) = thetaThresh;
    histogram(thetaTrlMean(grpTrlPos==p));
    title(sprintf('Theta %i', p));
    threshPosts{p,1,2} = grpTrlPosts(:,:,grpTrlPos==p & thetaTrlMean>thetaThresh(2));    
    threshPosts{p,2,2} = grpTrlPosts(:,:,grpTrlPos==p & thetaTrlMean<thetaThresh(1));
    
    
    tempLowBetaAccPos = nan(length(mlb.obsvTimeVect),mlb.seqLength);
    tempHighBetaAccPos = nan(length(mlb.obsvTimeVect),mlb.seqLength);
    tempLowThetaAccPos = nan(length(mlb.obsvTimeVect),mlb.seqLength);
    tempHighThetaAccPos = nan(length(mlb.obsvTimeVect),mlb.seqLength);
    for o = 1:mlb.seqLength
        tempLowBetaAccPos(:,o) = mean(grpTrlPosAccLog(:,grpTrlPos==p & betaTrlMean<betaThresh(1),o),2);
        tempHighBetaAccPos(:,o) = mean(grpTrlPosAccLog(:,grpTrlPos==p & betaTrlMean>betaThresh(2),o),2);
        tempLowThetaAccPos(:,o) = mean(grpTrlPosAccLog(:,grpTrlPos==p & thetaTrlMean<thetaThresh(1),o),2);
        tempHighThetaAccPos(:,o) = mean(grpTrlPosAccLog(:,grpTrlPos==p & thetaTrlMean>thetaThresh(2),o),2);
    end
    threshAccPos{p,1,1} = tempHighBetaAccPos;
    threshAccPos{p,2,1} = tempLowBetaAccPos;
    threshAccPos{p,1,2} = tempHighThetaAccPos;
    threshAccPos{p,2,2} = tempLowThetaAccPos;
        
    
    tempLowBetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempHighBetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempLowThetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempHighThetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    for t = 1:length(timePoints)
        tempLowBetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,grpTrlPos==p & betaTrlMean<betaThresh(1),t),2, 'omitnan');
        tempHighBetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,grpTrlPos==p & betaTrlMean>betaThresh(2),t),2);
        tempLowThetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,grpTrlPos==p & thetaTrlMean<thetaThresh(1),t),2);
        tempHighThetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,grpTrlPos==p & thetaTrlMean>thetaThresh(2),t),2);
    end
    threshAccTime{p,1,1} = tempHighBetaAccTime;
    threshAccTime{p,2,1} = tempLowBetaAccTime;
    threshAccTime{p,1,2} = tempHighThetaAccTime;
    threshAccTime{p,2,2} = tempLowThetaAccTime;
end
linkaxes;

%%
for band = 1:size(threshPosts,3)
    tempThresh = lfpThresh(:,:,band);
    if band==1
        tempLFPtrl = betaTrlMean;
    else
        tempLFPtrl = thetaTrlMean;
    end
    tempBandDecodes = cell(2,mlb.seqLength);
    tempTrlIDs = cell(2,mlb.seqLength);
    for pos = 1:mlb.seqLength
        tempBandDecodes{1,pos} = grpTrlDecodes(:,grpTrlPos==pos & tempLFPtrl<tempThresh(pos,1));
        tempTrlIDs{1,pos} = grpTrlPos(grpTrlPos==pos & tempLFPtrl<tempThresh(pos,1));
        tempBandDecodes{2,pos} = grpTrlDecodes(:,grpTrlPos==pos & tempLFPtrl>tempThresh(pos,2));
        tempTrlIDs{2,pos} = grpTrlPos(grpTrlPos==pos & tempLFPtrl>tempThresh(pos,2));
    end
    lowBandD= mlb.CalcDprmFromDecode(cell2mat(tempBandDecodes(1,:)), cell2mat(tempTrlIDs(1,:)));
    highBandD = mlb.CalcDprmFromDecode(cell2mat(tempBandDecodes(2,:)), cell2mat(tempTrlIDs(2,:)));
    dDiff = highBandD-lowBandD;
    figure; 
    h = gca;
    hold(h, 'on');
    for pos = 1:mlb.seqLength
        plot(mlb.obsvTimeVect,dDiff(:,pos), 'color', mlb.PositionColors(pos,:));
    end
    h.XRuler.FirstCrossoverValue = 0;
    h.XRuler.SecondCrossoverValue = 0;
    if band==1
        title('Beta');
    else
        title('Theta');
    end
end
        
%%
piNdx = find(abs(mlb.likeTimeVect)==min(abs(mlb.likeTimeVect)))+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
for band = 1:size(threshPosts,3)
    if band ==1
        bnm = 'Beta';
    elseif band == 2
        bnm = 'Theta';
    end
    tempPosts = threshPosts(:,:,band);
    curHighPowTrlPost = cell2mat(cellfun(@(a)mean(a,3),tempPosts(:,1),'uniformoutput',0)); cLim = [0 0.05];
%     curHighPowTrlPost = curHighPowTrlPost./max(curHighPowTrlPost(:));  cLim = [0 0.25];
%     curHighPowTrlPost = (curHighPowTrlPost-(mean(curHighPowTrlPost(:))))./std(curHighPowTrlPost(:));  cLim = [-6 6];    
    curLowPowTrlPost = cell2mat(cellfun(@(a)mean(a,3),tempPosts(:,2),'uniformoutput',0));
%     curLowPowTrlPost = curLowPowTrlPost./max(curLowPowTrlPost(:));
%     curLowPowTrlPost = (curLowPowTrlPost-(mean(curLowPowTrlPost(:))))./std(curLowPowTrlPost(:));
    imsp = nan(1,3);
    pasp = nan(1,3);
    tasp = nan(1,3);
    tdsp = nan(1,3);
    
    taClim = [0 0.15];
    
    figure;
    tasp(1) = subplot(3,9, [1,10]);
    highPowTimeAcc = cell2mat(threshAccTime(:,1,band));
    imagesc(highPowTimeAcc, taClim);
    set(gca, 'yticklabel', [], 'xticklabel', []);
    title('Time Decode');
    hold on;
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    for ndx = 1:length(piNdx)
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end    
    tdsp(1) = subplot(3,9,19);
    highPowTimeAccDiag = mean(cell2mat(cellfun(@(a){diag(a)}, permute(threshAccTime(:,1,band), [2,3,1]))),3);
    plot(highPowTimeAccDiag, 'k');
    hold on;    
    set(gca, 'ylim', [0 0.25], 'xticklabel', []);
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);    
    title('Mean Diag');
    imsp(1) = subplot(3,9,[2:3,11:12]);
    imagesc(curHighPowTrlPost, cLim);
    hold on;
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end
    title(sprintf('High Power %s', bnm));
    set(gca, 'yticklabel', [], 'xticklabel', []);
    pasp(1) = subplot(3,9,20:21);
    tempHighThresh = cell2mat(threshAccPos(:,1,band));
    hold on;
    for op = 1:mlb.seqLength
        plot(tempHighThresh(:,op), 'color', mlb.PositionColors(op,:));
    end
    axis tight;
    set(gca, 'ylim', [0 1], 'xticklabel', []);
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    title(sprintf('High %s Accuracy', bnm));
    
        
    
    tasp(2) = subplot(3,9, [4,13]);
    lowPowTimeAcc = cell2mat(threshAccTime(:,2,band));
    imagesc(lowPowTimeAcc, taClim);
    title('Time Decode');
    set(gca, 'yticklabel', [], 'xticklabel', []);
    hold on;
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    for ndx = 1:length(piNdx)
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end
    tdsp(2) = subplot(3,9,22);
    lowPowTimeAccDiag = mean(cell2mat(cellfun(@(a){diag(a)}, permute(threshAccTime(:,2,band), [2,3,1]))),3);
    plot(lowPowTimeAccDiag, 'k');
    hold on;    
    set(gca, 'ylim', [0 0.25], 'xticklabel', []);
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    title('Mean Diag');
    imsp(2) = subplot(3,9,[5:6,14:15]);
    imagesc(curLowPowTrlPost, cLim);
    hold on;
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end
    title(sprintf('Low Power %s', bnm));
    set(gca, 'yticklabel', [], 'xticklabel', []);
    pasp(2) = subplot(3,9,23:24);
    tempLowThresh = cell2mat(threshAccPos(:,2,band));
    hold on;
    for op = 1:mlb.seqLength
        plot(tempLowThresh(:,op), 'color', mlb.PositionColors(op,:));
    end
    axis tight;
    set(gca, 'ylim', [0 1], 'xticklabel', []);
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    title(sprintf('Low %s Accuracy', bnm));
    
    
    tasp(3) = subplot(3,9, [7,16]);
    timeAccDiff = highPowTimeAcc - lowPowTimeAcc;
    imagesc(timeAccDiff, [-0.05 0.05]);
    title('Time Decode');
    set(gca, 'yticklabel', [], 'xticklabel', []);
    hold on;
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    for ndx = 1:length(piNdx)
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end
    tdsp(3) = subplot(3,9,25);
    plot(highPowTimeAccDiag-lowPowTimeAccDiag, 'k');
    hold on;    
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(get(gca, 'xlim'), [0 0], '-k', 'linewidth', 1);
    set(gca, 'xticklabel', []);
    title('TimeDiag Diff');
    imsp(3) = subplot(3,9,[8:9,17:18]);
    imagesc(curHighPowTrlPost-curLowPowTrlPost, [max(cLim)*-1 max(cLim)]);
    hold on;
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 1);
        plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
        end
    end
    title(sprintf('%s Difference', bnm));
    set(gca, 'yticklabel', [], 'xticklabel', []);
    pasp(3) = subplot(3,9,26:27);
    tempHighMatch = nan(size(tempHighThresh));
    tempLowMatch = nan(size(tempLowThresh));
    boundNdx = [find(mlb.likeTimeVect==min(mlb.likeTimeVect)); length(mlb.likeTimeVect)+1];
    for p = 1:mlb.seqLength
        tempHighMatch(boundNdx(p):boundNdx(p+1)-1,p) = tempHighThresh(boundNdx(p):boundNdx(p+1)-1,p);
        tempLowMatch(boundNdx(p):boundNdx(p+1)-1,p) = tempLowThresh(boundNdx(p):boundNdx(p+1)-1,p);
    end
    tempDiff = tempHighMatch-tempLowMatch;
    hold on;
    for op = 1:mlb.seqLength
        plot(tempDiff(:,op), 'color', mlb.PositionColors(op,:), 'linewidth', 2);
        tmpLn = plot(tempHighThresh(:,op)-tempLowThresh(:,op), 'linewidth', 1);
        tmpLn.Color = [mlb.PositionColors(op,:), 0.5];
    end
    axis tight;
    set(gca, 'ylim', [-0.5 0.5], 'xticklabel', []);
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    plot(get(gca, 'xlim'), [0 0], '-k', 'linewidth', 1);
    title(sprintf('%s Accuracy Diff', bnm));
    
    linkaxes(pasp, 'x');
    linkaxes(tasp, 'xy');
    linkaxes(pasp(1:2), 'y');   
    linkaxes(tdsp(1:2), 'y');
    linkaxes([imsp, pasp], 'x');
    linkaxes([tasp, tdsp], 'x');
    linkaxes(imsp, 'xy');

    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s Split; %s aligned; Trial Window = (%.0fms:%.0fms); %i Perms', bnm, alignment{1}, trlWindow{1}(1), trlWindow{1}(2),numPerms),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
    