fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

binSize = 100;
dsRate = 5;
trlWindow = [0 1200];
numRandPerm = 50;
rng(sum(clock));
%%
zKEEPodorDECODEodor = cell(1,1,length(fileDirs));
zKEEPtimeDECODEodor = cell(1,1,length(fileDirs));
zKEEPnadaDECODEodor = cell(1,1,length(fileDirs));

zKEEPodorDECODEtime = cell(1,1,length(fileDirs));
zKEEPtimeDECODEtime = cell(1,1,length(fileDirs));
zKEEPnadaDECODEtime = cell(1,1,length(fileDirs));

zKEEPodorDECODEndx = cell(1,1,length(fileDirs));
zKEEPtimeDECODEndx = cell(1,1,length(fileDirs));
zKEEPnadaDECODEndx = cell(1,1,length(fileDirs));

for fl = 1:length(fileDirs)
    pfcMLB = MLB_SM(fileDirs{fl});
    pfcMLB.binSize = binSize;
    pfcMLB.dsRate = dsRate;
    [trlSpikes, trlTimeVect] = pfcMLB.PP_TrialMatrix_Spiking(trlWindow, 'PokeIn'); 
    pfcMLB.PP_IdentifyFISCseqs;
    
    %% Compile FISC likelihoods
    fiscLikes = nan(size(trlSpikes,1)*4, size(trlSpikes,2), size(pfcMLB.fiscTrials,2));
    for seq = 1:size(pfcMLB.fiscTrials,2)
        fiscLikes(:,:,seq) = [trlSpikes(:,:,pfcMLB.fiscTrials(1,seq));...
            trlSpikes(:,:,pfcMLB.fiscTrials(2,seq));...
            trlSpikes(:,:,pfcMLB.fiscTrials(3,seq));...
            trlSpikes(:,:,pfcMLB.fiscTrials(4,seq))];
    end
    timeVectLog = repmat(trlTimeVect, [4,1]);
    odorVectLog = [ones(size(trlTimeVect));...
        ones(size(trlTimeVect))*2;...
        ones(size(trlTimeVect))*3;...
        ones(size(trlTimeVect))*4];
    odorTimeNdxVect = (1:length(odorVectLog))';
    
    %% FISC Trials Real Decoding
    % Calculate FISC posteriors via Leave 1 Out
    fiscPosts = nan(size(timeVectLog,1), size(timeVectLog,1), size(pfcMLB.fiscTrials,2));
    for seq = 1:size(pfcMLB.fiscTrials,2)
        tempLikes = fiscLikes;
        tempLikes(:,:,seq) = [];
        fiscPosts(:,:,seq) = pfcMLB.CalcStaticBayesPost(nanmean(tempLikes,3), fiscLikes(:,:,seq));
    end
            
    % Decode FISC posteriors
    fiscOdorDecode = pfcMLB.DecodeBayesPost(fiscPosts, odorVectLog);
    fiscOdorDecoded = nan(size(timeVectLog,1),4);
    for op = 1:4
        fiscOdorDecoded(:,op) = mean(fiscOdorDecode==op,2);
    end
    fiscTimeDecoded = mean(pfcMLB.DecodeBayesPost(fiscPosts, timeVectLog) - repmat(timeVectLog, [1, size(pfcMLB.fiscTrials,2)]), 2, 'omitnan');
    fiscNdxDecoded = mean(pfcMLB.DecodeBayesPost(fiscPosts, odorTimeNdxVect) - repmat(odorTimeNdxVect, [1, size(pfcMLB.fiscTrials,2)]), 2, 'omitnan');
        
    %% FISC Trials Chance Decoding
    fiscKEEPodorDECODEodor = nan(size(timeVectLog,1), 4, numRandPerm);
    fiscKEEPtimeDECODEodor = nan(size(timeVectLog,1), 4, numRandPerm);
    fiscKEEPnadaDECODEodor = nan(size(timeVectLog,1), 4, numRandPerm);
    
    fiscKEEPodorDECODEtime = nan(size(timeVectLog,1), numRandPerm);
    fiscKEEPtimeDECODEtime = nan(size(timeVectLog,1), numRandPerm);
    fiscKEEPnadaDECODEtime = nan(size(timeVectLog,1), numRandPerm);
    
    fiscKEEPodorDECODEndx = nan(size(timeVectLog,1), numRandPerm);
    fiscKEEPtimeDECODEndx = nan(size(timeVectLog,1), numRandPerm);
    fiscKEEPnadaDECODEndx = nan(size(timeVectLog,1), numRandPerm);
    
    for randLikes = 1:numRandPerm
        % Calculate Time Shuffled Vector (preserve odor & shuffle time)
        %   Prediction: Should estimate temporal 
        tempRandLikesOdor = fiscLikes;
        tempRand_KEEPodor = nan(size(odorTimeNdxVect));
        for op = 1:4
            tempOdrNdxs = sortrows([(randperm(sum(odorVectLog==op)))',odorTimeNdxVect(odorVectLog==op)]);
            tempRand_KEEPodor(odorVectLog==op) = tempOdrNdxs(:,2);
        end
            
        % Calculate Odor Shuffled Vector (preserve time & shuffle odor ID)
        %   Prediction: Temporal coding falls to chance while Odor coding is preserved
        tempRandLikesTime = fiscLikes;
        tempRand_KEEPtime = nan(size(odorTimeNdxVect));
        for t = 1:length(trlTimeVect)
            tempTimeVectNdxs = find(timeVectLog==trlTimeVect(t));
            tempOdrIdVect = randperm(4);
            for op = 1:4
                tempRand_KEEPtime(tempTimeVectNdxs(op)) = odorTimeNdxVect(timeVectLog==trlTimeVect(t) & odorVectLog==tempOdrIdVect(op));
            end
        end
        
        % Calculate Fully Shuffled Vector (preserve nothing)
        %   Prediction: Should disrupt both odor and temporal coding 
        tempRandLikesNada = fiscLikes;
        tempRand_KEEPnada = randperm(length(odorTimeNdxVect))';
                
        % Apply the shuffled vectors to the relevant likelihoods
        for seq = 1:size(fiscLikes,3)
            tempKEEPodor = sortrows([tempRand_KEEPodor, tempRandLikesOdor(:,:,seq)]);
            tempRandLikesOdor(:,:,seq) = tempKEEPodor(:,2:end);
            
            tempKEEPtime = sortrows([tempRand_KEEPtime, tempRandLikesTime(:,:,seq)]);
            tempRandLikesTime(:,:,seq) = tempKEEPtime(:,2:end);
            
            tempKEEPnada = sortrows([tempRand_KEEPnada, tempRandLikesNada(:,:,seq)]);
            tempRandLikesNada(:,:,seq) = tempKEEPnada(:,2:end);
        end
        
        tempFISCpostKEEPodor = nan(size(fiscPosts));
        tempFISCpostKEEPtime = nan(size(fiscPosts));
        tempFISCpostKEEPnada = nan(size(fiscPosts));
        for seq = 1:size(fiscLikes,3)
            tempLikesKEEPodor = tempRandLikesOdor;
            tempLikesKEEPodor(:,:,seq) = [];
            tempFISCpostKEEPodor(:,:,seq) = pfcMLB.CalcStaticBayesPost(nanmean(tempLikesKEEPodor,3), fiscLikes(:,:,seq));
            
            tempLikesKEEPtime = tempRandLikesTime;
            tempLikesKEEPtime(:,:,seq) = [];
            tempFISCpostKEEPtime(:,:,seq) = pfcMLB.CalcStaticBayesPost(nanmean(tempLikesKEEPtime,3), fiscLikes(:,:,seq));
            
            tempLikesKEEPnada = tempRandLikesNada;
            tempLikesKEEPnada(:,:,seq) = [];
            tempFISCpostKEEPnada(:,:,seq) = pfcMLB.CalcStaticBayesPost(nanmean(tempLikesKEEPnada,3), fiscLikes(:,:,seq));            
        end
        
        % Decode FISC chance posteriors
        % Odor
        tempKEEPodorDECODEodor = pfcMLB.DecodeBayesPost(tempFISCpostKEEPodor, odorVectLog);
        tempKEEPtimeDEOCDEodor = pfcMLB.DecodeBayesPost(tempFISCpostKEEPtime, odorVectLog);
        tempKEEPnadaDECODEodor = pfcMLB.DecodeBayesPost(tempFISCpostKEEPnada, odorVectLog);
        for op = 1:4
            fiscKEEPodorDECODEodor(:,op,randLikes) = mean(tempKEEPodorDECODEodor==op,2);
            fiscKEEPtimeDECODEodor(:,op,randLikes) = mean(tempKEEPtimeDEOCDEodor==op,2);
            fiscKEEPnadaDECODEodor(:,op,randLikes) = mean(tempKEEPnadaDECODEodor==op,2);
        end
        % Time
        tempKEEPodorDECODEtime = pfcMLB.DecodeBayesPost(tempFISCpostKEEPodor, timeVectLog) - repmat(timeVectLog, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPodorDECODEtime(:,randLikes) = nanmean(tempKEEPodorDECODEtime,2);
        tempKEEPtimeDECODEtime = pfcMLB.DecodeBayesPost(tempFISCpostKEEPtime, timeVectLog) - repmat(timeVectLog, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPtimeDECODEtime(:,randLikes) = nanmean(tempKEEPtimeDECODEtime,2);
        tempKEEPnadaDECODEtime = pfcMLB.DecodeBayesPost(tempFISCpostKEEPnada, timeVectLog) - repmat(timeVectLog, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPnadaDECODEtime(:,randLikes) = nanmean(tempKEEPnadaDECODEtime,2);
        % Index
        tempKEEPodorDECODEndx = pfcMLB.DecodeBayesPost(tempFISCpostKEEPodor, odorTimeNdxVect) - repmat(odorTimeNdxVect, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPodorDECODEndx(:,randLikes) = nanmean(tempKEEPodorDECODEndx,2);
        tempKEEPtimeDECODEndx = pfcMLB.DecodeBayesPost(tempFISCpostKEEPtime, odorTimeNdxVect) - repmat(odorTimeNdxVect, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPtimeDECODEndx(:,randLikes) = nanmean(tempKEEPtimeDECODEndx,2);
        tempKEEPnadaDECODEndx = pfcMLB.DecodeBayesPost(tempFISCpostKEEPnada, odorTimeNdxVect) - repmat(odorTimeNdxVect, [1, size(pfcMLB.fiscTrials,2)]);
        fiscKEEPnadaDECODEndx(:,randLikes) = nanmean(tempKEEPnadaDECODEndx,2);
    end
    %% Z-Score
    zTempKoDo = nan(size(fiscOdorDecoded));
    zTempKtDo = nan(size(fiscOdorDecoded));
    zTempKnDo = nan(size(fiscOdorDecoded));
    for t = 1:size(fiscOdorDecoded,1)
        for op = 1:size(fiscOdorDecoded,2)
            meanKoDo = mean(fiscKEEPodorDECODEodor(t,op,:),3, 'omitnan');
            stdKoDo = std(fiscKEEPodorDECODEodor(t,op,:),0,3,'omitnan');
            zTempKoDo(t,op) = (fiscOdorDecoded(t,op)-meanKoDo)/stdKoDo;
            
            meanKtDo = mean(fiscKEEPtimeDECODEodor(t,op,:),3,'omitnan');
            stdKtDo = std(fiscKEEPtimeDECODEodor(t,op,:),0,3,'omitnan');
            zTempKtDo(t,op) = (fiscOdorDecoded(t,op)-meanKtDo)/stdKtDo;
            
            meanKnDo = mean(fiscKEEPnadaDECODEodor(t,op,:),3,'omitnan');
            stdKnDo = std(fiscKEEPnadaDECODEodor(t,op,:),0,3,'omitnan');
            zTempKnDo(t,op) = (fiscOdorDecoded(t,op)-meanKnDo)/stdKnDo;
        end
    end
    zKEEPodorDECODEodor{fl} = zTempKoDo;
    zKEEPtimeDECODEodor{fl} = zTempKtDo;
    zKEEPnadaDECODEodor{fl} = zTempKnDo;
    
    zTempKoDt = nan(size(fiscOdorDecoded));
    zTempKtDt = nan(size(fiscOdorDecoded));
    zTempKnDt = nan(size(fiscOdorDecoded));
    for t = 1:size(fiscTimeDecoded,1)
        meanKoDt = mean(fiscKEEPodorDECODEtime(t,:), 'omitnan');
        stdKoDt = std(fiscKEEPodorDECODEtime(t,:), 'omitnan');
        zTempKoDt(t) = (fiscTimeDecoded(t)-meanKoDt)/stdKoDt;
        
        meanKtDt = mean(fiscKEEPtimeDECODEtime(t,:), 'omitnan');
        stdKtDt = std(fiscKEEPtimeDECODEtime(t,:), 'omitnan');
        zTempKtDt(t) = (fiscTimeDecoded(t)-meanKtDt)/stdKtDt;
        
        meanKnDt = mean(fiscKEEPnadaDECODEtime(t,:), 'omitnan');
        stdKnDt = std(fiscKEEPnadaDECODEtime(t,:), 'omitnan');
        zTempKnDt(t) = (fiscTimeDecoded(t)-meanKnDt)/stdKnDt;
    end    
    zKEEPodorDECODEtime{fl} = zTempKoDt;
    zKEEPtimeDECODEtime{fl} = zTempKtDt;
    zKEEPnadaDECODEtime{fl} = zTempKnDt;
    
    zTempKoDn = nan(size(fiscOdorDecoded));
    zTempKtDn = nan(size(fiscOdorDecoded));
    zTempKnDn = nan(size(fiscOdorDecoded));
    for t = 1:size(fiscTimeDecoded,1)
        meanKoDn = mean(fiscKEEPodorDECODEndx(t,:), 'omitnan');
        stdKoDn = std(fiscKEEPodorDECODEndx(t,:), 'omitnan');
        zTempKoDn(t) = (fiscNdxDecoded(t)-meanKoDn)/stdKoDn;
        
        meanKtDn = mean(fiscKEEPtimeDECODEndx(t,:), 'omitnan');
        stdKtDn = std(fiscKEEPtimeDECODEndx(t,:), 'omitnan');
        zTempKtDn(t) = (fiscNdxDecoded(t)-meanKtDn)/stdKtDn;
        
        meanKnDn = mean(fiscKEEPnadaDECODEndx(t,:), 'omitnan');
        stdKnDn = std(fiscKEEPnadaDECODEndx(t,:), 'omitnan');
        zTempKnDn(t) = (fiscNdxDecoded(t)-meanKnDn)/stdKnDn;        
    end
    zKEEPodorDECODEndx{fl} = zTempKoDn;
    zKEEPtimeDECODEndx{fl} = zTempKtDn;
    zKEEPnadaDECODEndx{fl} = zTempKnDn;
    
    
    %% Calculate ISC posteriors
    odorLog = [pfcMLB.trialInfo.Odor];
    posLog = [pfcMLB.trialInfo.Position];
    perfLog = [pfcMLB.trialInfo.Performance]==1;
    iscLog = (odorLog==posLog) & perfLog;
    iscTrls = trlSpikes(:,:,iscLog);
    iscPosts = nan(size(iscTrls,1),size(timeVectLog,1),size(iscTrls,3));
    for trl = 1:size(iscTrls,3)
        iscPosts(:,:,trl) = pfcMLB.CalcStaticBayesPost(nanmean(fiscLikes,3), iscTrls(:,:,trl));
    end
    iscOdorLog = odorLog(iscLog);
    iscPosLog = posLog(iscLog);
    %% Decode FISC posteriors
    iscOdorDecode = pfcMLB.DecodeBayesPost(iscPosts, odorVectLog);
    iscTimeDecode = pfcMLB.DecodeBayesPost(iscPosts, timeVectLog) - repmat(trlTimeVect, [1, size(iscPosts,3)]);
end
%%
zKEEPodorDECODEodor = cell2mat(zKEEPodorDECODEodor);
zKEEPtimeDECODEodor = cell2mat(zKEEPtimeDECODEodor);
zKEEPnadaDECODEodor = cell2mat(zKEEPnadaDECODEodor);

zKEEPodorDECODEtime = cell2mat(zKEEPodorDECODEtime);
zKEEPtimeDECODEtime = cell2mat(zKEEPtimeDECODEtime);
zKEEPnadaDECODEtime = cell2mat(zKEEPnadaDECODEtime);

zKEEPodorDECODEndx = cell2mat(zKEEPodorDECODEndx);
zKEEPtimeDECODEndx = cell2mat(zKEEPtimeDECODEndx);
zKEEPnadaDECODEndx = cell2mat(zKEEPnadaDECODEndx);

%% Plot Group Plots
figure; 
spO(1) = subplot(3,3,1);
for o = 1:4
    meanDecode = mean(zKEEPodorDECODEodor(:,o,:),3, 'omitnan');
    stdDecode = std(zKEEPodorDECODEodor(:,o,:),0,3,'omitnan')./sqrt(size(zKEEPodorDECODEodor,3)-1);
    plot(meanDecode, 'color', pfcMLB.PositionColors(o,:), 'linewidth', 2);
    hold on;
    patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
        pfcMLB(1).PositionColors(o,:), 'FaceAlpha', 0.5, 'EdgeColor', pfcMLB(1).PositionColors(o,:));
end
title('Keep Odor - Permute Time');
ylabel([{'Odor/Position Decoding'};{'Z-Norm to Chance'}]);

spO(2) = subplot(3,3,2);
for o = 1:4
    meanDecode = mean(zKEEPtimeDECODEodor(:,o,:),3, 'omitnan');
    stdDecode = std(zKEEPtimeDECODEodor(:,o,:),0,3, 'omitnan')./sqrt(size(zKEEPtimeDECODEodor,3)-1);
    plot(meanDecode, 'color', pfcMLB.PositionColors(o,:), 'linewidth', 2);
    hold on;
    patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
        pfcMLB(1).PositionColors(o,:), 'FaceAlpha', 0.5, 'EdgeColor', pfcMLB(1).PositionColors(o,:));
end
title('Keep Time - Permute Odor');

spO(3) = subplot(3,3,3);
for o = 1:4
    meanDecode = mean(zKEEPnadaDECODEodor(:,o,:),3, 'omitnan');
    stdDecode = std(zKEEPnadaDECODEodor(:,o,:),0,3, 'omitnan')./sqrt(size(zKEEPnadaDECODEodor,3)-1);
    plot(meanDecode, 'color', pfcMLB.PositionColors(o,:), 'linewidth', 2);
    hold on;
    patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
        pfcMLB(1).PositionColors(o,:), 'FaceAlpha', 0.5, 'EdgeColor', pfcMLB(1).PositionColors(o,:));
end
title('Permute Odor & Time');

linkaxes(spO, 'xy');
for sp = 1:length(spO)
    zeroNdxs = find(timeVectLog==0);
    for op = 1:4
        line(spO(sp),repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
    end
end

spT(1) = subplot(3,3,4);
plot(mean(zKEEPodorDECODEtime,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;
ylabel([{'Temporal Discrepancy'};{'Z-Norm to Chance'}]);

spT(2) = subplot(3,3,5);
plot(mean(zKEEPtimeDECODEtime,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;

spT(3) = subplot(3,3,6);
plot(mean(zKEEPnadaDECODEtime,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;

linkaxes(spT, 'xy');
for sp = 1:length(spT)
    zeroNdxs = find(timeVectLog==0);
    for op = 1:4
        line(spT(sp),repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
    end
end

spN(1) = subplot(3,3,7);
plot(mean(zKEEPodorDECODEndx,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;
zeroNdxs = find(timeVectLog==0);
for op = 1:4
    line(repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
end
ylabel([{'Absolute Discrepancy'};{'Z-Norm to Chance'}]);

spN(2) = subplot(3,3,8);
plot(mean(zKEEPtimeDECODEndx,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;
zeroNdxs = find(timeVectLog==0);
for op = 1:4
    line(repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
end

spN(3) = subplot(3,3,9);
plot(mean(zKEEPnadaDECODEndx,3, 'omitnan'), 'color', 'k', 'linewidth', 2);
hold on;
zeroNdxs = find(timeVectLog==0);
for op = 1:4
    line(repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
end

linkaxes(spN, 'xy');
for sp = 1:length(spN)
    zeroNdxs = find(timeVectLog==0);
    for op = 1:4
        line(spN(sp),repmat(zeroNdxs(op), [1,2]), get(gca, 'ylim'), 'linestyle', '--', 'color', 'k', 'linewidth', 2);
    end
end
