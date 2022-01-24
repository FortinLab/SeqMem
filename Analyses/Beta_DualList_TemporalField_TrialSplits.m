clear all; %#ok<CLALL>
%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'},...
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];

fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual List Sessions\GE11_Session146'},...
    {'D:\WorkBigDataFiles\PFC\Dual List Sessions\GE13_Session103'},...
    {'D:\WorkBigDataFiles\PFC\Dual List Sessions\GE17_Session110'}];

binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
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

%% Data Outputs
% Behavior Variables
fiscPokeOutLat = cell(length(fileDirs),1);
fiscRwdDelivLat = cell(length(fileDirs),1);
% smi = nan(length(fileDirs),2);
% dPrm = nan(length(fileDirs),2);
% ri = nan(length(fileDirs),2);
% smiByOP = nan(length(fileDirs),8,2);
% dPrmByOP = nan(length(fileDirs),8,2);
% riByOP = nan(length(fileDirs),8,2);
% Neural Variables
aniPosts = cell(size(fileDirs));
aniOdorDecodes = cell(size(fileDirs));
aniTimeDecodes = cell(size(fileDirs));
aniPosAccLog = cell(size(fileDirs));
aniTimeAccLog = cell(size(fileDirs));
aniTrlPosIDs = cell(size(fileDirs));
aniTrlOdrIDs = cell(size(fileDirs));
aniLFP = cell(size(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%     smi(ani,:) = mlb.smi;
%     dPrm(ani,:) = mlb.dPrime;
%     ri(ani,:) = mlb.ri;
%     for op = 1:mlb.seqLength
%         smiByOP(ani,:,1) = mlb.smiByPos;
%         smiByOP(ani,:,2) = mlb.smiByOdr;
%         dPrmByOP(ani,:,1) = mlb.dPrimeByPos;
%         dPrmByOP(ani,:,2) = mlb.dPrimeByOdr;
%         riByOP(ani,:,1) = mlb.riByPos;
%         riByOP(ani,:,2) = mlb.riByOdr;
%     end
    %% Extract LFP
    [~, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
    [~, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
        %% Decode ISC via Sub-Sampling
%     mlb.bayesType = 1;  % Comment In to decode using Gaussian rather than Poisson
    mlb.SetLikes_SubSample;
%     mlb.bayesType = 3;  % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_Observes;
    tempPosts = cell2mat(reshape(mlb.post, [1,1,numel(mlb.post)]));
    postTrlIDs = cell2mat(reshape(mlb.postTrlIDs, [1,1,numel(mlb.postTrlIDs)]));
    tempOdorDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,4));
    odrVect = unique(mlb.decodeIDvects{1}(:,4));
    tempPosAccuracy = arrayfun(@(a,b)a==b,repmat(tempOdorDecode,[1,1,length(odrVect)]), repmat(permute(odrVect,[3,2,1]), size(tempOdorDecode,1),size(tempOdorDecode,2),1));
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
    aniTrlOdrIDs{ani} = [mlb.trialInfo(postTrlIDs).Odor];
    aniLFP{ani} = tempLFP;
    %% Free up memory
    clear tempPosts tempOdorDecode tempTimeDecode tempPosAccuracy tempTimeAccuracy tempLFP
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
grpTrlOdr = cell2mat(aniTrlOdrIDs);
grpTrlPos = cell2mat(aniTrlPosIDs);
grpTrlPosAccLog = cell2mat(aniPosAccLog);
grpTrlTimeAccLog = cell2mat(aniTimeAccLog);
betaTrlMean = mean(grpTrlLFP(trlTimeLog,:,1));
thetaTrlMean = mean(grpTrlLFP(trlTimeLog ,:,2));

lfpThresh = nan(mlb.seqLength,2,2);
threshPosts = cell(mlb.seqLength,2,2);
threshAccPos = cell(mlb.seqLength,2,2);
threshAccTime = cell(mlb.seqLength,2,2);

figure;
for o = 1:length(odrVect)
    curOdr = odrVect(o);
    [~,curPos] = find(curOdr==mlb.odrSeqs);
    curOdrLog = grpTrlOdr==curOdr;
    curPosLog = grpTrlPos==curPos;
    subplot(length(odrVect),2, sub2ind([2,length(odrVect)], 1,o));
%     betaThresh = [mean(betaTrlMean(curPosLog))-(std(betaTrlMean(curPosLog))*1), mean(betaTrlMean(curPosLog))+(std(betaTrlMean(curPosLog))*1)];              % POSITION threshold (Mean+/-STD)
%     betaThresh = [mean(betaTrlMean(curOdrLog))-(std(betaTrlMean(curOdrLog))*1), mean(betaTrlMean(curOdrLog))+(std(betaTrlMean(curOdrLog))*1)];              % ODOR threshold (Mean+/-STD)
%     betaThresh = [mean(betaTrlMean)-(std(betaTrlMean(grpTrlOdr))*1), mean(betaTrlMean)+(std(betaTrlMean)*1)];                                               % SESSION threshold (Mean +/-STD)
%     sortedBeta = sort(betaTrlMean(curPosLog)); betaThresh = [sortedBeta(ceil(length(sortedBeta)*0.25)); sortedBeta(floor(length(sortedBeta)*0.75))];        % POSITION threshold (<25/>75)
    sortedBeta = sort(betaTrlMean(curOdrLog)); betaThresh = [sortedBeta(ceil(length(sortedBeta)*0.25)); sortedBeta(floor(length(sortedBeta)*0.75))];        % ODOR threshold (<25/>75)
%     sortedBeta = sort(betaTrlMean); betaThresh = [sortedBeta(ceil(length(sortedBeta)*0.25)); sortedBeta(floor(length(sortedBeta)*0.75))];                   % SESSION threshold (<25/>75)
    lfpThresh(o,:,1) = betaThresh;
    histogram(betaTrlMean(curOdrLog), 'Normalization', 'cdf'); 
    hold on;
    title(sprintf('Beta %i', curOdr));
    threshPosts{o,1,1} = grpTrlPosts(:,:,curOdrLog & betaTrlMean>betaThresh(2));
    threshPosts{o,2,1} = grpTrlPosts(:,:,curOdrLog & betaTrlMean<betaThresh(1));
        
    subplot(length(odrVect),2, sub2ind([2,length(odrVect)], 2,o));
%     thetaThresh = [mean(thetaTrlMean(curPosLog))-(std(thetaTrlMean(curPosLog))*1), mean(thetaTrlMean(curPosLog))+(std(thetaTrlMean(curPosLog))*1)];         % POSITION threshold (Mean+/-STD)
%     thetaThresh = [mean(thetaTrlMean(curOdrLog))-(std(thetaTrlMean(curOdrLog))*1), mean(thetaTrlMean(curOdrLog))+(std(thetaTrlMean(curOdrLog))*1)];         % ODOR threshold (Mean+/-STD)
%     thetaThresh = [mean(thetaTrlMean)-(std(thetaTrlMean)*1), mean(thetaTrlMean)+(std(thetaTrlMean)*1)];                                                     % SESSION threshold (Mean +/-STD)
%     sortedTheta = sort(thetaTrlMean(curPosLog)); thetaThresh = [sortedTheta(ceil(length(sortedTheta)*0.25)); sortedTheta(floor(length(sortedTheta)*0.75))]; % POSITION threshold (<25/>75)
    sortedTheta = sort(thetaTrlMean(curOdrLog)); thetaThresh = [sortedTheta(ceil(length(sortedTheta)*0.25)); sortedTheta(floor(length(sortedTheta)*0.75))]; % ODOR threshold (<25/>75)
%     sortedTheta = sort(thetaTrlMean); thetaThresh = [sortedTheta(ceil(length(sortedTheta)*0.25)); sortedTheta(floor(length(sortedTheta)*0.75))];            % SESSION threshold (<25/>75)
    lfpThresh(o,:,2) = thetaThresh;
    histogram(thetaTrlMean(curOdrLog), 'Normalization', 'cdf');
    hold on;
    title(sprintf('Theta %i', curOdr));
    threshPosts{o,1,2} = grpTrlPosts(:,:,curOdrLog & thetaTrlMean>thetaThresh(2));    
    threshPosts{o,2,2} = grpTrlPosts(:,:,curOdrLog & thetaTrlMean<thetaThresh(1));
    
    
    tempLowBetaAccPos = nan(length(mlb.obsvTimeVect),length(odrVect));
    tempHighBetaAccPos = nan(length(mlb.obsvTimeVect),length(odrVect));
    tempLowThetaAccPos = nan(length(mlb.obsvTimeVect),length(odrVect));
    tempHighThetaAccPos = nan(length(mlb.obsvTimeVect),length(odrVect));
    for odr = 1:length(odrVect)
        tempLowBetaAccPos(:,odr) = mean(grpTrlPosAccLog(:,curOdrLog & betaTrlMean<betaThresh(1),odr),2);
        tempHighBetaAccPos(:,odr) = mean(grpTrlPosAccLog(:,curOdrLog & betaTrlMean>betaThresh(2),odr),2);
        tempLowThetaAccPos(:,odr) = mean(grpTrlPosAccLog(:,curOdrLog & thetaTrlMean<thetaThresh(1),odr),2);
        tempHighThetaAccPos(:,odr) = mean(grpTrlPosAccLog(:,curOdrLog & thetaTrlMean>thetaThresh(2),odr),2);
    end
    threshAccPos{o,1,1} = tempHighBetaAccPos;
    threshAccPos{o,2,1} = tempLowBetaAccPos;
    threshAccPos{o,1,2} = tempHighThetaAccPos;
    threshAccPos{o,2,2} = tempLowThetaAccPos;
        
    
    tempLowBetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempHighBetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempLowThetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempHighThetaAccTime = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    for t = 1:length(timePoints)
        tempLowBetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,curOdrLog & betaTrlMean<betaThresh(1),t),2, 'omitnan');
        tempHighBetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,curOdrLog & betaTrlMean>betaThresh(2),t),2);
        tempLowThetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,curOdrLog & thetaTrlMean<thetaThresh(1),t),2);
        tempHighThetaAccTime(:,t) = mean(grpTrlTimeAccLog(:,curOdrLog & thetaTrlMean>thetaThresh(2),t),2);
    end
    threshAccTime{o,1,1} = tempHighBetaAccTime;
    threshAccTime{o,2,1} = tempLowBetaAccTime;
    threshAccTime{o,1,2} = tempHighThetaAccTime;
    threshAccTime{o,2,2} = tempLowThetaAccTime;
end
linkaxes;

for o = 1:length(odrVect)
    subplot(length(odrVect),2, sub2ind([2,length(odrVect)], 1,o));
    plot(repmat(lfpThresh(o,1,1),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
    plot(repmat(lfpThresh(o,2,1),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
    subplot(length(odrVect),2, sub2ind([2,length(odrVect)], 2,o));
    plot(repmat(lfpThresh(o,1,2),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
    plot(repmat(lfpThresh(o,2,2),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);    
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
    for op = 1:length(odrVect)
        [curSeq,curPos] = find(odrVect(op)==mlb.odrSeqs);
        if curSeq==1
            plot(tempHighThresh(:,op), 'color', mlb.PositionColors(curPos,:));
        else
            plot(tempHighThresh(:,op), 'color', mlb.PositionColors(curPos,:), 'linestyle', '--');
        end
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
    for op = 1:length(odrVect)
        [curSeq,curPos] = find(odrVect(op)==mlb.odrSeqs);
        if curSeq==1
            plot(tempLowThresh(:,op), 'color', mlb.PositionColors(curPos,:));
        else
            plot(tempLowThresh(:,op), 'color', mlb.PositionColors(curPos,:), 'linestyle', '--');
        end
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
    for p = 1:length(odrVect)
        tempHighMatch(boundNdx(p):boundNdx(p+1)-1,p) = tempHighThresh(boundNdx(p):boundNdx(p+1)-1,p);
        tempLowMatch(boundNdx(p):boundNdx(p+1)-1,p) = tempLowThresh(boundNdx(p):boundNdx(p+1)-1,p);
    end
    tempDiff = tempHighMatch-tempLowMatch;
    hold on;
    for op = 1:length(odrVect)
        [curSeq,curPos] = find(odrVect(op)==mlb.odrSeqs);
        if curSeq==1
            plot(tempDiff(:,op), 'color', mlb.PositionColors(curPos,:), 'linewidth', 2);
            tmpLn = plot(tempHighThresh(:,op)-tempLowThresh(:,op), 'linewidth', 1);
            tmpLn.Color = [mlb.PositionColors(curPos,:), 0.5];
        else            
            plot(tempDiff(:,op), 'color', mlb.PositionColors(curPos,:), 'linestyle', '--', 'linewidth', 2);
            tmpLn = plot(tempHighThresh(:,op)-tempLowThresh(:,op), 'linewidth', 1);
            tmpLn.Color = [mlb.PositionColors(curPos,:), 0.5];
        end
        
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
    
    figure;    
    sps = nan(mlb.seqLength,2);
    for pos = 1:mlb.seqLength
        seq1OdrLog = odrVect==mlb.odrSeqs(1,pos);
        seq2OdrLog = odrVect==mlb.odrSeqs(2,pos);
        
        seq1 = cell2mat(cellfun(@(a)a(:,seq1OdrLog), threshAccPos(seq1OdrLog,:,band), 'uniformoutput', 0)) - cell2mat(cellfun(@(a)a(:,seq2OdrLog), threshAccPos(seq1OdrLog,:,band), 'uniformoutput', 0));
        seq2 = cell2mat(cellfun(@(a)a(:,seq1OdrLog), threshAccPos(seq2OdrLog,:,band), 'uniformoutput', 0)) - cell2mat(cellfun(@(a)a(:,seq2OdrLog), threshAccPos(seq2OdrLog,:,band), 'uniformoutput', 0));
        
        sps(pos,1) = subplot(2,mlb.seqLength, sub2ind([mlb.seqLength, 2], pos,1));
        h = gca;
        h.XRuler.FirstCrossoverValue = 0;
        h.XRuler.SecondCrossoverValue = 0;
        hold(h, 'on');
        plot(mlb.obsvTimeVect, seq1(:,1), 'color', mlb.PositionColors(pos,:), 'linestyle', '-');        
        plot(mlb.obsvTimeVect, seq1(:,2), 'color', [mlb.PositionColors(pos,:), 0.5], 'linestyle', '-');       
        plot(mlb.obsvTimeVect, seq2(:,1), 'color', mlb.PositionColors(pos,:), 'linestyle', '--');        
        plot(mlb.obsvTimeVect, seq2(:,2), 'color', [mlb.PositionColors(pos,:), 0.5], 'linestyle', '--');      
        
        sps(pos,2) = subplot(2,mlb.seqLength, sub2ind([mlb.seqLength, 2], pos,2));
        h = gca;
        h.XRuler.FirstCrossoverValue = 0;
        h.XRuler.SecondCrossoverValue = 0;
        hold(h, 'on');
        plot(mlb.obsvTimeVect, diff(fliplr(seq1),1,2), 'color', 'k', 'linestyle', '-');        
        plot(mlb.obsvTimeVect, diff(fliplr(seq2),1,2), 'color', 'k', 'linestyle', '--');             
    end
    linkaxes(sps);
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s Split; %s aligned; Trial Window = (%.0fms:%.0fms); %i Perms', bnm, alignment{1}, trlWindow{1}(1), trlWindow{1}(2),numPerms),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
    