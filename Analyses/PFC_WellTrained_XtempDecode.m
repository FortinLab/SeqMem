% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% 
fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% 
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
ssProportion = 0.5;
numPerms = 10;
ssType = 0; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
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
smiByOP = nan(length(fileDirs),4,2);
dPrmByOP = nan(length(fileDirs),4,2);
riByOP = nan(length(fileDirs),4,2);
% Neural Variables
grpPost = cell(size(fileDirs));
grpDecode = cell(size(fileDirs));
grpAccLog = cell(size(fileDirs));
grpLFP = cell(size(fileDirs));
grpTrlPos = cell(size(fileDirs));
grpDecodeMean = cell(4,4,length(fileDirs));
grpTPdecode = cell(1,4,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
    uniInfo = mlb.unitInfo;
    betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
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
    %% Process Neural Data
    [betaPhase, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
    [thetaPhase, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
%     mlb.SetLikes_FISC;
    mlb.SetLikes_SubSample;
    mlb.Process_IterativeObserves;
    grpPost{ani} = mlb.post{1};
    [decodes, maxDecode] = mlb.DecodeBayesPost(mlb.post{1});   
    grpDecode{ani} = decodes;
    trlLFP = nan(size(decodes,1),4,size(decodes,3));
    tempAccLog = nan(size(decodes,1),size(decodes,2), size(decodes,3), mlb.seqLength);
    for trl = 1:size(decodes,3)
        trlLFP(:,1,trl) = betaPower(:,:,mlb.obsvTrlIDs{1}(trl));
        trlLFP(:,2,trl) = thetaPower(:,:,mlb.obsvTrlIDs{1}(trl));
        for op = 1:mlb.seqLength
            tempAccLog(:,:,trl,op) = decodes(:,:,trl)==op;
        end
    end
    grpLFP{ani} = trlLFP;
    grpTrlPos{ani} = [mlb.trialInfo(mlb.obsvTrlIDs{1}).Position];
    grpAccLog{ani} = tempAccLog;
    figure;
    for tPos = 1:mlb.seqLength
        tPosLog = [mlb.trialInfo(mlb.obsvTrlIDs{1}).Position]==tPos;
        for dPos = 1:mlb.seqLength
            grpDecodeMean{tPos,dPos,ani} = mean(decodes(:,:,tPosLog)==dPos,3);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], dPos, tPos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, grpDecodeMean{tPos,dPos,ani}, [0 0.5]);
            set(gca, 'ydir', 'normal');
            drawnow;
            title(sprintf('%i Decode %i', tPos, dPos));
        end
    end
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', mlb.nsmblMatFile,...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    for op = 1:mlb.seqLength
        tempPosDecode = grpDecodeMean(:,op,ani);
        tempHit = tempPosDecode{op};
        tempFA = mean(cell2mat(reshape(tempPosDecode((1:mlb.seqLength)~=op), [1,1,mlb.seqLength-1])),3);
        grpTPdecode{1,op,ani} = tempHit-tempFA;
    end
%     figure;
%     sps = nan(1,mlb.seqLength);
%     spss = nan(1,mlb.seqLength);
%     spsss = nan(1,mlb.seqLength);
%     for p = 1:mlb.seqLength
%         tempPos = grpDecodeMean(:,p,ani);
%         tempMatch = tempPos{p};
%         tempFA = mean(cell2mat(reshape(tempPos((1:mlb.seqLength)~=p), [1,1,mlb.seqLength-1])),3);
%         subplot(4,mlb.seqLength,sub2ind([mlb.seqLength,3],p,1));
%         imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempMatch-tempFA, [0 0.5]);
%         grpTPdecode{1,p,ani} = tempMatch-tempFA;
%         set(gca, 'ydir', 'normal');
%         tpInt = mlb.IntegrateAntiDiagonal(tempMatch-tempFA);
% %         tpInt = mlb.QuantPersist(tempMatch-tempFA,0.15);
% %         faInt = mlb.IntegrateAntiDiagonal(tempFA);
%         sps(p) = subplot(4,mlb.seqLength, sub2ind([mlb.seqLength,4],p,2));
%         yyaxis(gca, 'left');
%         plot(mlb.obsvTimeVect, tpInt, '-k');
%         hold on;
%         yyaxis(gca, 'right');
%         plot(mlb.obsvTimeVect, mean(trlLFP(:,2,trlLFP(1,4,:)==p),3), '-r');
%         plot(mlb.obsvTimeVect, mean(trlLFP(:,3,trlLFP(1,4,:)==p),3), '--r');
%         
%         tempTrlTI = trlLFP(mlb.obsvTimeVect>=0,1,trlLFP(1,4,:)==p);
% %         tempTrlTI = tempTrlTI-repmat(faInt', [1,1,size(tempTrlTI,3)]);
%         tempTrlBeta = trlLFP(mlb.obsvTimeVect>=0,2,trlLFP(1,4,:)==p);
%         tempTrlTheta = trlLFP(mlb.obsvTimeVect>=0,3,trlLFP(1,4,:)==p);
%         spss(p) = subplot(4,mlb.seqLength, sub2ind([mlb.seqLength,4],p,3));
%         corrScatPlot(tempTrlBeta(:), tempTrlTI(:), 'Beta Power', 'TI');
%         spsss(p) = subplot(4,mlb.seqLength, sub2ind([mlb.seqLength,4],p,4));
%         corrScatPlot(tempTrlTheta(:), tempTrlTI(:), 'Theta Power', 'TI');
% %     end
%     linkaxes(sps, 'xy');
%     linkaxes(spss, 'xy');
%     linkaxes(spsss, 'xy');
end    

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%%
figure;
for p = 1:mlb.seqLength
    subplot(1,4,p);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(grpTPdecode(1,p,:)), 3, 'omitnan'), [0 0.5]);
    set(gca, 'ydir', 'normal');
    hold on;
    plot(get(gca, 'xlim'), get(gca, 'ylim'), '-k', 'linewidth', 1);
end
annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("Cross-Temporal Hit-FA: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

figure; 
imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(reshape(grpTPdecode(1,2:4,:), [1,1,numel(grpTPdecode(1,2:4,:))])), 3, 'omitnan'), [0 0.35]);
set(gca, 'ydir', 'normal')
hold on;
plot([0 0], get(gca, 'ylim'), '--k');
plot(get(gca, 'xlim'),[0 0], '--k');
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
plot(get(gca, 'ylim'),repmat(nearestPOtime, [1,2]), '--k', 'linewidth', 1);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
plot(get(gca, 'ylim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
plot(get(gca, 'xlim'), get(gca, 'ylim'), '-k', 'linewidth', 1);

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("Cross-Temporal Hit-FA Mean (nonA): binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%%
grpBetaPwrSplit = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
grpThetaPwrSplit = cell(mlb.seqLength,mlb.seqLength, length(fileDirs));
for ani = 1:length(fileDirs)
    curAniPost = grpPost{ani};
    curAniDecode = grpDecode{ani};
    curAniLFP = grpLFP{ani};
    curTrlPos = grpTrlPos{ani};
    
    %% Evaluate animal trial-wise LFP power distributions
    trlPrdLog = mlb.obsvTimeVect>=0;
    trlLFP = mean(curAniLFP(trlPrdLog,:,:),1);
    figure; 
    subplot(1,2,1);
    histogram(trlLFP(1,1,:), -1:0.1:3);
    hold on;
    plot(repmat(median(trlLFP(1,1,:)), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(mean(trlLFP(1,1,:)), [1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
    title('Beta Power');
    subplot(1,2,2);
    histogram(trlLFP(1,2,:), -1:0.1:3);
    hold on;
    plot(repmat(median(trlLFP(1,2,:)), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(repmat(mean(trlLFP(1,2,:)), [1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
    title('Theta Power');
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    %%
%     highBetaTrlLog = reshape(trlLFP(1,1,:)>median(trlLFP(1,1,:)), [1, size(trlLFP,3)]); lowBetaTrlLog = reshape(trlLFP(1,1,:)<median(trlLFP(1,1,:)), [1, size(trlLFP,3)]);
%     highBetaTrlLog = reshape(trlLFP(1,1,:)>(mean(trlLFP(1,1,:))+std(trlLFP(1,1,:))), [1, size(trlLFP,3)]);    lowBetaTrlLog = reshape(trlLFP(1,1,:)<(mean(trlLFP(1,1,:))-std(trlLFP(1,1,:))), [1, size(trlLFP,3)]);
    highBetaTrlLog = reshape(trlLFP(1,1,:)>median(trlLFP(1,1,trlLFP(1,1,:)>median(trlLFP(1,1,:)))), [1, size(trlLFP,3)]);    lowBetaTrlLog = reshape(trlLFP(1,1,:)<median(trlLFP(1,1,:)), [1, size(trlLFP,3)]);
%     highThetaTrlLog = reshape(trlLFP(1,2,:)>median(trlLFP(1,2,:)), [1,size(trlLFP,3)]); lowThetaTrlLog = reshape(trlLFP(1,2,:)<median(trlLFP(1,2,:)), [1,size(trlLFP,3)]);
%     highThetaTrlLog = reshape(trlLFP(1,2,:)>(mean(trlLFP(1,2,:))), [1,size(trlLFP,3)]);    lowThetaTrlLog = reshape(trlLFP(1,2,:)<(mean(trlLFP(1,2,:))-std(trlLFP(1,2,:))), [1,size(trlLFP,3)]);
    highThetaTrlLog = reshape(trlLFP(1,2,:)>median(trlLFP(1,2,trlLFP(1,2,:)>median(trlLFP(1,2,:)))), [1,size(trlLFP,3)]);    lowThetaTrlLog = reshape(trlLFP(1,2,:)<median(trlLFP(1,2,:)), [1,size(trlLFP,3)]);
    betaPwrSplitDecode = cell(mlb.seqLength,mlb.seqLength,2);
    betaDiff = cell(mlb.seqLength);
    thetaPwrSplitDecode = cell(mlb.seqLength,mlb.seqLength,2);
    thetaDiff = cell(mlb.seqLength);
    betaD = figure; 
    betaH = figure;
    betaL = figure;
    thetaD = figure; 
    thetaH = figure; 
    thetaL = figure; 
    for pos = 1:4            
        for odr = 1:4
            betaPwrSplitDecode{odr,pos,1} = mean(curAniPost(:,:,odr, curTrlPos==pos & highBetaTrlLog),4);
            betaPwrSplitDecode{odr,pos,2} = mean(curAniPost(:,:,odr, curTrlPos==pos & lowBetaTrlLog),4);
            betaDiff{odr,pos} =  betaPwrSplitDecode{odr,pos,1}-betaPwrSplitDecode{odr,pos,2};
            figure(betaD);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, betaPwrSplitDecode{odr,pos,1}-betaPwrSplitDecode{odr,pos,2}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
            figure(betaH);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, betaPwrSplitDecode{odr,pos,1}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
            figure(betaL);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, betaPwrSplitDecode{odr,pos,2}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
            
            thetaPwrSplitDecode{odr,pos,1} = mean(curAniPost(:,:,odr, curTrlPos==pos & highThetaTrlLog),4);
            thetaPwrSplitDecode{odr,pos,2} = mean(curAniPost(:,:,odr, curTrlPos==pos & lowThetaTrlLog),4);
            thetaDiff{odr,pos} =  thetaPwrSplitDecode{odr,pos,1}-thetaPwrSplitDecode{odr,pos,2};
            figure(thetaD);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, thetaPwrSplitDecode{odr,pos,1}-thetaPwrSplitDecode{odr,pos,2}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
            figure(thetaH);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, thetaPwrSplitDecode{odr,pos,1}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
            figure(thetaL);
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
            imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, thetaPwrSplitDecode{odr,pos,2}, [-0.5 0.5]);
            title(sprintf('%i Decode %i', pos, odr));
            set(gca, 'ydir', 'normal');
        end
    end    
    annotation(betaD,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Beta Diff: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(betaH,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('High Beta Trials: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(betaL,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Low Beta Trials: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(thetaD,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Theta Diff: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(thetaH,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('High Theta Trials: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(thetaL,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Low Theta Trials: %s', fileDirs{ani}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    grpBetaPwrSplit(:,:,ani) = betaDiff;
    grpThetaPwrSplit(:,:,ani) = thetaDiff;
end

%%
% close all;
betaD = figure;
thetaD = figure;
for pos = 1:4
    for odr = 1:4
        figure(betaD);
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(grpBetaPwrSplit(odr,pos,:)),3, 'omitnan'), [-0.2 0.2]);
        set(gca, 'ydir', 'normal');
        title(sprintf('%i Decode %i', pos, odr));
        
        figure(thetaD);
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength,mlb.seqLength], odr, pos));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(grpThetaPwrSplit(odr,pos,:)),3, 'omitnan'), [-0.2 0.2]);
        set(gca, 'ydir', 'normal');
        title(sprintf('%i Decode %i', pos, odr));
    end
end
annotation(betaD,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Beta Diff Average',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
annotation(thetaD,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Theta Diff average',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%%
aniMatchDiffBeta = cell(1,1,length(fileDirs));
aniMatchDiffTheta = cell(1,1,length(fileDirs));
betaAni = figure;
thetaAni = figure;
for ani = 1:length(fileDirs)
    aniPOtime = mean(fiscPokeOutLat{ani});
    aniRWDtime = mean(fiscRwdDelivLat{ani});
    tempBeta = grpBetaPwrSplit(:,:,ani);
    tempTheta = grpThetaPwrSplit(:,:,ani);
    aniMatchDiffBeta{ani} = mean(cell2mat(reshape(tempBeta(logical(eye(mlb.seqLength))), [1,1,mlb.seqLength])),3, 'omitnan');
%     tempTP = cell(1,1,mlb.seqLength);
%     for op = 1:mlb.seqLength
%         tempPosDecode = tempBeta(:,op);
%         tempHit = tempPosDecode{op};
%         tempFA = mean(cell2mat(reshape(tempPosDecode((1:mlb.seqLength)~=op), [1,1,mlb.seqLength-1])),3, 'omitnan');
%         tempTP{1,op} = tempHit-tempFA;
%     end
%     aniMatchDiffBeta{ani} = mean(cell2mat(tempTP),3);
    figure(betaAni);
    subplot(1,length(fileDirs), ani);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, aniMatchDiffBeta{ani}, [-0.2 0.2]);
    hold on;
    plot([0 0], get(gca, 'ylim'), '--k');
    plot(get(gca, 'xlim'),[0 0], '--k');
    plot(repmat(aniPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
    plot(get(gca, 'ylim'),repmat(aniPOtime, [1,2]), '--k', 'linewidth', 1);
    plot(repmat(aniRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(get(gca, 'ylim'), repmat(aniRWDtime, [1,2]), ':k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
    
    aniMatchDiffTheta{ani} = mean(cell2mat(reshape(tempTheta(logical(eye(mlb.seqLength))), [1,1,mlb.seqLength])),3, 'omitnan');
%     tempTP = cell(1,1,mlb.seqLength);
%     for op = 1:mlb.seqLength
%         tempPosDecode = tempTheta(:,op);
%         tempHit = tempPosDecode{op};
%         tempFA = mean(cell2mat(reshape(tempPosDecode((1:mlb.seqLength)~=op), [1,1,mlb.seqLength-1])),3, 'omitnan');
%         tempTP{1,op} = tempHit-tempFA;
%     end
%     aniMatchDiffTheta{ani} = mean(cell2mat(tempTP),3);
    figure(thetaAni);
    subplot(1,length(fileDirs), ani);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, aniMatchDiffTheta{ani}, [-0.2 0.2]);
    hold on;
    plot([0 0], get(gca, 'ylim'), '--k');
    plot(get(gca, 'xlim'),[0 0], '--k');
    plot(repmat(aniPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
    plot(get(gca, 'ylim'),repmat(aniPOtime, [1,2]), '--k', 'linewidth', 1);
    plot(repmat(aniRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(get(gca, 'ylim'), repmat(aniRWDtime, [1,2]), ':k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
end
annotation(betaAni,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Animal Beta Diff Average',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
annotation(thetaAni,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Animal Theta Diff average',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%%
figure;
subplot(1,2,1);
% imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,mean(cell2mat(aniMatchDiffBeta),3).*std(cell2mat(aniMatchDiffBeta),0,3), [-0.01 0.01]);
imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,mean(cell2mat(aniMatchDiffBeta),3, 'omitnan'), [-0.2 0.2]);
hold on;
plot([0 0], get(gca, 'ylim'), '--k');
plot(get(gca, 'xlim'),[0 0], '--k');
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
plot(get(gca, 'ylim'),repmat(nearestPOtime, [1,2]), '--k', 'linewidth', 1);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
plot(get(gca, 'ylim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
plot(get(gca, 'xlim'),get(gca, 'ylim'), '-k', 'linewidth', 1);
set(gca, 'ydir', 'normal');
title('Beta');
subplot(1,2,2);
% imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,mean(cell2mat(aniMatchDiffTheta),3).*std(cell2mat(aniMatchDiffTheta),0,3), [-0.01 0.01]);
imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,mean(cell2mat(aniMatchDiffTheta),3, 'omitnan'), [-0.2 0.2]);
hold on;
plot([0 0], get(gca, 'ylim'), '--k');
plot(get(gca, 'xlim'),[0 0], '--k');
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 1);
plot(get(gca, 'ylim'),repmat(nearestPOtime, [1,2]), '--k', 'linewidth', 1);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
plot(get(gca, 'ylim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
plot(get(gca, 'xlim'),get(gca, 'ylim'), '-k', 'linewidth', 1);
set(gca, 'ydir', 'normal');
title('Theta')

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf('%s aligned', alignment{1}),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%%
