clear all;
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
setupSeqLength = 4; 

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];
% tets = [1,22,17,18,17];
% setupSeqLength = 5;

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4; 
% 
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4;  

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
% fiscLOO analysis
fiscLOO_Posts = cell(size(fileDirs));
fiscLOO_Decodes = cell(3,1,length(fileDirs));
fiscLOO_DecodeTime = [];
fiscLOO_ObserveTime = [];
% fiscISC analysis
fiscISC_Posts = cell(size(fileDirs));
fiscISC_Decodes = cell(3,1,length(fileDirs));
fiscISC_DecodeTime = [];
fiscISC_ObserveTime = [];
% iscBS analysis
iscBS_Posts = cell(size(fileDirs));
iscBS_Decodes = cell(3,1,length(fileDirs));
iscBS_DecodeTime = [];
iscBS_ObserveTime = [];
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
    uniInfo = mlb.unitInfo;
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
    mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
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
    fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
    fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
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
    %% Decode FISC via Leave-1-Out
%     trialIDs = [{'Time'}, {'Window'}, {'Position'}, {'Odor'}];
    mlb.SetLikes_FISC;
%     mlb.bayesType = 3; % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_LikelyLOO;
    fiscLOO_Posts(ani) = mlb.post;
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints));
    for t = 1:length(timePoints)
        timeAccuracy(:,t) = mean(timeDecode==timePoints(t),2, 'omitnan');
    end
    fiscLOO_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions));
    for p = 1:length(positions)
        posAccuracy(:,p) = mean(posDecode==positions(p),2, 'omitnan');
    end
    fiscLOO_Decodes{2,1,ani} = posAccuracy;
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors));
    for p = 1:length(odors)
        odrAccuracy(:,p) = mean(odrDecode==odors(p),2, 'omitnan');
    end
    fiscLOO_Decodes{3,1,ani} = odrAccuracy;    
    if ani==1
        fiscLOO_DecodeTime = mlb.likeTimeVect;
    end
%     
    %% Decode ISC via FISC 
    mlb.Process_Observes;
    trlOdrVect = [mlb.trialInfo(reshape(mlb.postTrlIDs{1}, [1,numel(mlb.postTrlIDs{1})])).Odor];
    trlPosVect = [mlb.trialInfo(reshape(mlb.postTrlIDs{1}, [1,numel(mlb.postTrlIDs{1})])).Position];
    trlPoss = unique(trlPosVect);
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints), length(trlPoss));
    timePosts = nan(size(mlb.post{1},1), size(mlb.post{1},2), length(trlPoss));
    for pos = 1:length(trlPoss)
        timePosts(:,:,pos) = mean(mlb.post{1}(:,:,trlPosVect==trlPoss(pos)),3, 'omitnan');
        for t = 1:length(timePoints)
            timeAccuracy(:,t,pos) = mean(timeDecode(:,trlPosVect==trlPoss(pos))==timePoints(t),2, 'omitnan');
        end
    end
    fiscISC_Posts{ani} = timePosts;
    fiscISC_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions), length(positions));
    for pos = 1:length(positions)
        for p = 1:length(positions)
            posAccuracy(:,p,pos) = mean(posDecode(:,trlPosVect==positions(pos))==positions(p),2, 'omitnan');
        end
    end
    fiscISC_Decodes{2,1,ani} = posAccuracy; 
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors), length(odors));
    for odr = 1:length(odors)
        for p = 1:length(odors)
            odrAccuracy(:,p,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(p),2, 'omitnan');
        end
    end
    fiscISC_Decodes{3,1,ani} = odrAccuracy;    
    
    if ani==1
        fiscISC_DecodeTime = mlb.likeTimeVect;
        fiscISC_ObserveTime = mlb.obsvTimeVect;
    end
    
    %% Decode ISC via Sub-Sampling
%     mlb.bayesType = 1;  % Comment In to decode using Gaussian rather than Poisson
    mlb.SetLikes_SubSample;
%     mlb.bayesType = 3;  % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_Observes;
    tempPosts = cell2mat(reshape(mlb.post, [1,1,numel(mlb.post)]));
    postTrlIDs = cell2mat(reshape(mlb.postTrlIDs, [1,1,numel(mlb.postTrlIDs)]));
    
    trlOdrVect = [mlb.trialInfo(postTrlIDs).Odor];
    trlOdrs = unique(trlOdrVect);
    trlPosVect = [mlb.trialInfo(postTrlIDs).Position];
    trlPoss = unique(trlPosVect);
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints), length(trlPoss));
    timePosts = nan(size(tempPosts,1), size(tempPosts,2), length(trlPoss));
    for pos = 1:length(trlPoss)
        timePosts(:,:,pos) = mean(tempPosts(:,:,trlPosVect==trlPoss(pos)),3,'omitnan'); 
        for t = 1:length(timePoints)
            timeAccuracy(:,t,pos) = mean(timeDecode(:,trlPosVect==trlPoss(pos))==timePoints(t),2, 'omitnan');
        end
    end
    iscBS_Posts{ani} = timePosts;
    iscBS_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions), length(positions));
    for pos = 1:length(positions)
        for p = 1:length(positions)
            posAccuracy(:,p,pos) = mean(posDecode(:,trlPosVect==positions(pos))==positions(p),2, 'omitnan');
        end
    end
    iscBS_Decodes{2,1,ani} = posAccuracy; 
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors), length(odors));
    for odr = 1:length(odors)
        for p = 1:length(odors)
            odrAccuracy(:,p,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(p),2, 'omitnan');
        end
    end
    iscBS_Decodes{3,1,ani} = odrAccuracy;     
    
    if ani==1
        iscBS_DecodeTime = mlb.likeTimeVect;
        iscBS_ObserveTime = mlb.obsvTimeVect;
    end
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%% Plot FISC Leave-One-Out
figure;
sp1 = subplot(5,5,[2:5,7:10,12:15, 17:20]);
imagesc(mean(cell2mat(reshape(fiscLOO_Posts, [1,1,numel(fiscLOO_Posts)])),3, 'omitnan'), postCLim);
set(gca, 'ydir', 'normal', 'ytick', [], 'xtick', []);
hold on
piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
for ndx = 1:length(piNdx)
    plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);

    if ndx<length(piNdx)
        plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
    end
end
title('Posteriors');

sp2 = subplot(5,5,1:5:16);
imagesc(mlb.obsvTimeVect, 1:length(mlb.likeTimeVect), mean(cell2mat(fiscLOO_Decodes(1,1,:)),3), decodeCLim);
hold on;
plot([0 0], get(gca,'ylim'), '--k', 'linewidth', 2);
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
for ndx = 1:length(piNdx)
    plot(get(gca, 'xlim'), repmat(piNdx(ndx), [1,2]), '--k', 'linewidth', 2);
    plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
    if ndx<length(piNdx)
        plot(get(gca, 'xlim'), repmat(posNdx(ndx),[1,2]), '-k', 'linewidth', 2);
    end
end    
set(gca, 'ydir', 'normal', 'ytick', []);
title('Temporal Accuracy');
ylabel('True Time');

sp3 = subplot(5,5,22:25);
tempPosAcc = cell2mat(fiscLOO_Decodes(2,1,:));
for p = 1:size(fiscLOO_Decodes{2,1,1},2)
    meanPosAcc = mean(tempPosAcc(:,p,:),3, 'omitnan');
    varPosAcc = mlb.SEMcalc(tempPosAcc(:,p,:),0,3);
    plot(meanPosAcc, 'color', mlb.PositionColors(p,:));
    hold on;
    patch('XData', [1:size(tempPosAcc,1), size(tempPosAcc,1):-1:1],...
        'YData', [(meanPosAcc+varPosAcc)', flipud(meanPosAcc-varPosAcc)'], 'edgecolor', mlb.PositionColors(p,:),...
        'facecolor', mlb.PositionColors(p,:), 'facealpha', 0.5);
    piNdx = find(mlb.likeTimeVect==0);
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]),[0 1], '--k','linewidth', 1);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 1);
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
end
set(gca, 'xtick', [], 'color', 'none');
xlabel('Decoded Time');
title('Position Accuracy');
box off
axis tight
linkaxes([sp1 sp3],'x');
linkaxes([sp1 sp2],'y');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("FISC Leave-One-Out Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%%
opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    tempTimeNdx = find(mlb.likeTimeVect==mlb.obsvTimeVect(t));
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = tempPosAcc(tempTimeNdx(op),op,ani);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,1), [numel(smiByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,1), [numel(dPrmByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,1), [numel(riByOP(:,:,1)),1]), 'Rows', 'pairwise');
end
figure; 
subplot(3,1,1)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Position');


opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    tempTimeNdx = find(mlb.likeTimeVect==mlb.obsvTimeVect(t));
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = tempPosAcc(tempTimeNdx(op),op,ani);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,2), [numel(smiByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,2), [numel(dPrmByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,2), [numel(riByOP(:,:,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,2)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odor');

opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength-1);
    tempTimeNdx = find(mlb.likeTimeVect==mlb.obsvTimeVect(t));
    for ani = 1:length(fileDirs)
        for op = 2:mlb.seqLength
            tempDecode(ani,op-1) = tempPosAcc(tempTimeNdx(op),op,ani);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,2:end,2), [numel(smiByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,2:end,2), [numel(dPrmByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,2:end,2), [numel(riByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,3)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odors B-D');

linkaxes;
for sp = 1:3
    subplot(3,1,sp)
    plot(get(gca, 'xlim'), [0 0], '-k');
    plot([0 0], get(gca, 'ylim'), '--k');
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 1);
end
legend([{'SMI'}, {'d'''}, {'RI'}]);


annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("FISC Leave-One-Out Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot FISC ISC Analysis
figure;
tempPosAcc = fiscISC_Decodes(2,1,:);
sps = nan(1,4);
for p = 1:size(fiscISC_Posts{1},3)
    tempPost = mean(cell2mat(reshape(cellfun(@(a)a(:,:,p),fiscISC_Posts,'uniformoutput', 0), [1,1,length(fileDirs)])),3, 'omitnan');
    subplot(5,5,(2:5)+(5*(size(fiscISC_Posts{1},3)-p)));
    imagesc(tempPost, postCLim);
    hold on
    plot(get(gca, 'xlim'), repmat(find(mlb.obsvTimeVect==0)+0.5, [1,2]), '--k', 'linewidth', 2);
    set(gca, 'ydir', 'normal', 'xtick', [], 'ytick', []);
    piNdx = find(mlb.likeTimeVect==0)+0.5;
    poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
    rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    if p==size(fiscISC_Posts{1},3)
        title('Posteriors');
    end
    
    subplot(5,5,5*(size(fiscISC_Posts{1},3)-p)+1)
    imagesc(mlb.obsvTimeVect,mlb.obsvTimeVect,mean(cell2mat(cellfun(@(a)a(:,:,p), fiscISC_Decodes(1,1,:), 'uniformoutput',0)),3), decodeCLim);
    hold on;
    plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '--k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
    ylabel(sprintf('%s/%i', mlb.Rosetta{p}, p));
    if p==size(fiscISC_Posts{1},3)
        title('Temporal Accuracy');
    end

    sps(p) = subplot(5,5,(5*4)+p+1);
    tempAcc = cell2mat(reshape(cellfun(@(a)a(:,:,p),tempPosAcc, 'uniformoutput', 0), [1,1,length(fileDirs)]));
    for o = 1:size(tempAcc,2)
        meanPosAcc = mean(tempAcc(:,o,:),3, 'omitnan');
        varPosAcc = mlb.SEMcalc(tempAcc(:,o,:),0,3);
        plot(mlb.obsvTimeVect, meanPosAcc, 'color', mlb.PositionColors(o,:));
        hold on;
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)]',...
            'YData', [(meanPosAcc+varPosAcc)', flipud(meanPosAcc-varPosAcc)'], 'edgecolor', mlb.PositionColors(o,:),...
            'facecolor', mlb.PositionColors(o,:), 'facealpha', 0.5);
        set(gca, 'xtick', [], 'color', 'none');
        box off
        axis tight
        title(sprintf('%s/%i', mlb.Rosetta{p}, p));
        for ndx = 1:length(piNdx)
            plot([0 0], [0 1], '--k');
            plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
            plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 1);
        end
    end
end
linkaxes(sps, 'xy');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("FISC Decode ISC Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');


%%
opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = fiscISC_Decodes{2,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,1), [numel(smiByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,1), [numel(dPrmByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,1), [numel(riByOP(:,:,1)),1]), 'Rows', 'pairwise');
end
figure; 
subplot(3,1,1)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Position');


opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = fiscISC_Decodes{3,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,2), [numel(smiByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,2), [numel(dPrmByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,2), [numel(riByOP(:,:,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,2)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odor');

opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength-1);
    for ani = 1:length(fileDirs)
        for op = 2:mlb.seqLength
            tempDecode(ani,op-1) = fiscISC_Decodes{3,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,2:end,2), [numel(smiByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,2:end,2), [numel(dPrmByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,2:end,2), [numel(riByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,3)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odors B-D');

linkaxes;
for sp = 1:3
    subplot(3,1,sp)
    plot(get(gca, 'xlim'), [0 0], '-k');
    plot([0 0], get(gca, 'ylim'), '--k');
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 1);
end
legend([{'SMI'}, {'d'''}, {'RI'}]);


annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("FISC Decode ISC Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot ISC Bootstrapped analysis
figure;
tempPosAcc = iscBS_Decodes(2,1,:);
sps = nan(1,4);
for p = 1:size(iscBS_Posts{1},3)
    tempPost = mean(cell2mat(reshape(cellfun(@(a)a(:,:,p),iscBS_Posts,'uniformoutput', 0), [1,1,length(fileDirs)])),3, 'omitnan');
    subplot(5,5,(2:5)+(5*(size(fiscISC_Posts{1},3)-p)));
    imagesc(tempPost, postCLim);
    hold on
    plot(get(gca, 'xlim'), repmat(find(mlb.obsvTimeVect==0)+0.5, [1,2]), '--k', 'linewidth', 2);
    set(gca, 'ydir', 'normal', 'xtick', [], 'ytick', []);
    ylabel(sprintf('%s/%i', mlb.Rosetta{p}, p));
    piNdx = find(mlb.likeTimeVect==0)+0.5;
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    if p==size(fiscISC_Posts{1},3)
        title('Posteriors');
    end

    subplot(5,5,5*(size(iscBS_Posts{1},3)-p)+1)
    imagesc(mlb.obsvTimeVect,mlb.obsvTimeVect,mean(cell2mat(cellfun(@(a)a(:,:,p), iscBS_Decodes(1,1,:), 'uniformoutput',0)),3), decodeCLim);
    hold on;
    plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '--k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
    ylabel(sprintf('%s/%i', mlb.Rosetta{p}, p));
    if p==size(iscBS_Posts{1},3)
        title('Temporal Accuracy');
    end

    sps(p) = subplot(5,5,(5*4)+p+1);
    tempAcc = cell2mat(reshape(cellfun(@(a)a(:,:,p),tempPosAcc, 'uniformoutput', 0), [1,1,length(fileDirs)]));
    for o = 1:size(tempAcc,2)
        meanPosAcc = mean(tempAcc(:,o,:),3, 'omitnan');
        varPosAcc = mlb.SEMcalc(tempAcc(:,o,:),0,3);
        plot(mlb.obsvTimeVect, meanPosAcc, 'color', mlb.PositionColors(o,:));
        hold on;
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)]',...
            'YData', [(meanPosAcc+varPosAcc)', flipud(meanPosAcc-varPosAcc)'], 'edgecolor', mlb.PositionColors(o,:),...
            'facecolor', mlb.PositionColors(o,:), 'facealpha', 0.5);
        set(gca, 'xtick', [], 'color', 'none');
        box off
        axis tight
        title(sprintf('%s/%i', mlb.Rosetta{p}, p));
        for ndx = 1:length(piNdx)
            plot([0 0], [0 1], '--k');
            plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
            plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 1);
        end
    end
end
linkaxes(sps, 'xy');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Subsampled ISC Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');


%%
opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = iscBS_Decodes{2,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,1), [numel(smiByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,1), [numel(dPrmByOP(:,:,1)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,1), [numel(riByOP(:,:,1)),1]), 'Rows', 'pairwise');
end
figure; 
subplot(3,1,1)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Position');


opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength);
    for ani = 1:length(fileDirs)
        for op = 1:mlb.seqLength
            tempDecode(ani,op) = iscBS_Decodes{3,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,:,2), [numel(smiByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,:,2), [numel(dPrmByOP(:,:,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,:,2), [numel(riByOP(:,:,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,2)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odor');

opSMIcorrVect = nan(size(mlb.obsvTimeVect));
opDPcorrVect = nan(size(mlb.obsvTimeVect));
opRIcorrVect = nan(size(mlb.obsvTimeVect));
for t = 1:length(mlb.obsvTimeVect)
    tempDecode = nan(length(fileDirs), mlb.seqLength-1);
    for ani = 1:length(fileDirs)
        for op = 2:mlb.seqLength
            tempDecode(ani,op-1) = iscBS_Decodes{3,1,ani}(t,op,op);
        end
    end    
    opSMIcorrVect(t) = corr(tempDecode(:), reshape(smiByOP(:,2:end,2), [numel(smiByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opDPcorrVect(t) = corr(tempDecode(:), reshape(dPrmByOP(:,2:end,2), [numel(dPrmByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
    opRIcorrVect(t) = corr(tempDecode(:), reshape(riByOP(:,2:end,2), [numel(riByOP(:,2:end,2)),1]), 'Rows', 'pairwise');
end
subplot(3,1,3)
plot(mlb.obsvTimeVect, opSMIcorrVect, '-k', 'linewidth', 2);
hold on;
plot(mlb.obsvTimeVect, opDPcorrVect, '--k', 'linewidth', 2);
plot(mlb.obsvTimeVect, opRIcorrVect, ':k', 'linewidth', 2);
title('Behavior Correlates vs Odors B-D');

linkaxes;
for sp = 1:3
    subplot(3,1,sp)
    plot(get(gca, 'xlim'), [0 0], '-k');
    plot([0 0], get(gca, 'ylim'), '--k');
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 1);
end
legend([{'SMI'}, {'d'''}, {'RI'}]);

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Subsampled ISC Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');


