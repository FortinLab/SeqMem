% PFC_WellTrained_Group_PositionChance
%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% setupSeqLength = 4; 

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];
% tets = [1,22,17,18,17];
% setupSeqLength = 5;

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
setupSeqLength = 4; 
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
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

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
% Posteriors
realPost_L1O = cell(setupSeqLength,1,length(fileDirs));
chancePost_L1O = cell(setupSeqLength,numChancePerms,length(fileDirs));
realPost_ISC = cell(setupSeqLength,1,length(fileDirs));
chancePost_ISC = cell(setupSeqLength, numChancePerms,length(fileDirs));

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
    mlb.bayesType = bayesType;
    
    %% Extract Behavioral Variables
    fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeInIndex])'/1000;
    fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeInIndex])'/1000;
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
    %% FISC via Leave-1-Out
    mlb.SetLikes_FISC;
%     mlb.bayesType = 3; % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_LikelyL1O;
    realPost_L1O(:,:,ani) = mlb.post;
      
    %% ISC via FISC 
    mlb.Process_Observes;
    trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
    trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
    trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
    tempPosts = cell(mlb.seqLength,1);
    for pos = 1:mlb.seqLength
        curISClog = trlOdrVect==pos & trlPosVect==pos & trlPerfVect==1;
        tempPosts{pos} = mlb.post(:,:,curISClog);
    end
    realPost_ISC(:,:,ani) = tempPosts;
    %% Calculate Chance
    for perm = 1:numChancePerms 
        chancePerm = permute(reshape(randperm(numel(mlb.fiscTrials)), size(mlb.fiscTrials)), [1,3,2]);
        mlb.likeTrlSpikes = mlb.likeTrlSpikes(chancePerm);
        mlb.Process_LikelyL1O;
        chancePost_L1O(:,perm,ani) = mlb.post;
        mlb.Process_Observes;
        trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
        trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
        trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
        tempPosts = cell(mlb.seqLength,1);
        for pos = 1:mlb.seqLength
            curISClog = trlOdrVect==pos & trlPosVect==pos & trlPerfVect==1;
            tempPosts{pos} = mlb.post(:,:,curISClog);
        end
        chancePost_ISC(:,perm,ani) = tempPosts;
    end
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

save('PFC_WT_Trial_Group_Posts.mat', 'realPost_L1O', 'chancePost_L1O', 'realPost_ISC', 'chancePost_ISC', 'mlb', '-v7.3');
% 
% %% Collapse Trials Across Animals
% realPostL1O = cell(mlb.seqLength,1,3);
% realPostMargTimeL1O = cell(mlb.seqLength,1,3);
% realPostMargPosL1O = cell(mlb.seqLength,1,3);
% 
% realPostISC = cell(mlb.seqLength,1,2);
% realPostMargTimeISC = cell(mlb.seqLength,1,2);
% realPostMargPosISC = cell(mlb.seqLength,1,2);
% 
% chancePostL1O = cell(mlb.seqLength,1,permChance);
% chancePostMargTimeL1O = cell(mlb.seqLength,1,permChance);
% chancePostMargPosL1O = cell(mlb.seqLength,1,permChance);
% 
% chancePostISC = cell(mlb.seqLength,1,permChance);
% chancePostMargTimeISC = cell(mlb.seqLength,1,permChance);
% chancePostMargPosISC = cell(mlb.seqLength,1,permChance);
% for pos = 1:mlb.seqLength
%     tempPostL1O = cell2mat(realPost_L1O(pos,:,:));
%     realPostL1O{pos,1,1} = mean(tempPostL1O,3);
%     realPostL1O{pos,1,2} = std(tempPostL1O,0,3);
%     realPostL1O{pos,1,3} = sum(~isnan(tempPostL1O),3);
%     
%     tempPostISC = cell2mat(realPost_ISC(pos,:,:));
%     realPostISC{pos,1,1} = mean(tempPostISC,3);
%     realPostISC{pos,1,2} = std(tempPostISC,0,3);
%     realPostISC{pos,1,3} = sum(~isnan(tempPostL1O),3);
%     
%     for perm = 1:numChancePerms 
%         chancePostL1O{pos,1,perm} = mean(cell2mat(chancePost_L1O(pos,perm,:)),3);
%         chancePostISC{pos,1,perm} = mean(cell2mat(chancePost_ISC(pos,perm,:)),3);
%     end
% end
% 
% 
% realDecodeL1O = nan(length(mlb.likeTimeVect));
% realDecodeTimeL1O = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect));
% realDecodePosL1O = nan(length(mlb.likeTimeVect), mlb.seqLength);
% 
% realDecodeISC = nan(length(mlb.likeTimeVect));
% realDecodeTimeISC = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect));
% realDecodePosISC = nan(length(mlb.likeTimeVect), mlb.seqLength);
% 
% chancePostL1O = nan(length(mlb.likeTimeVect),length(mlb.likeTimeVect), 2);
% chancePostMargTimeL10 = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), 2);
% chancePostMargPosL1O = nan(length(mlb.likeTimeVect), mlb.seqLength, 2);
% chanceDecodeL1O = nan(length(mlb.likeTimeVect),length(mlb.likeTimeVect),2);
% chanceDecodeTimeL1O = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), 2);
% chanceDecodePosL1O = nan(length(mlb.likeTimeVect), mlb.seqLength, 2);
% 
% chancePostISC = nan(length(mlb.likeTimeVect),length(mlb.likeTimeVect), 2);
% chancePostMargTimeISC = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), 2);
% chancePostMargPosISC = nan(length(mlb.likeTimeVect), mlb.seqLength, 2);
% chanceDecodeISC = nan(length(mlb.likeTimeVect),length(mlb.likeTimeVect),2);
% chanceDecodeTimeISC = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), 2);
% chanceDecodePosISC = nan(length(mlb.likeTimeVect), mlb.seqLength, 2);
% 
% 
% 
% 
% %% Process L1O Chance Data
% % FISC L1O Posterior Values
% fiscChancePost = cell(mlb.seqLength, numChancePerms );
% for c = 1:numChancePerms 
%     for o = 1:mlb.seqLength
%         fiscChancePost{o,c} = mean(cell2mat(chance_fiscL1O_Post(o,c,:)),3);
%     end
% end
% clear chance_fiscL1O_Post
% 
% % FISC L1O Time-Position Coding
% fiscChanceTimePosAcc = nan(length(mlb.likeTimeVect), length(mlb.likeTimeVect), numChancePerms );
% for c = 1:numChancePerms 
%     tempChance = cell2mat(chance_fiscL1O_TimePos(:,:,c));
%     for t = 1:length(mlb.likeTimeVect)
%         for ndx = 1:length(mlb.likeTimeVect)
%             fiscChanceTimePosAcc(t,ndx,c) = mean(tempChance(t,:)==ndx);
%         end
%     end
% end
% clear chance_fiscL1O_TimePos
% 
% % FISC L1O Time Coding
% fiscChanceTimeAcc = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), numChancePerms );
% for c = 1:numChancePerms 
%     tempChance = cell2mat(chance_fiscL1O_MargTime(:,:,c));
%     for t = 1:length(mlb.likeTimeVect)
%         for ndx = 1:length(mlb.obsvTimeVect)
%             fiscChanceTimeAcc(t,ndx,c) = mean(tempChance(t,:)==mlb.obsvTimeVect(ndx));
%         end
%     end
% end
% clear chance_fiscL1O_MargTime
% 
% % FISC L1O Odor Coding
% fiscChancePosAcc = nan(length(mlb.likeTimeVect), mlb.seqLength, numChancePerms );
% for c = 1:numChancePerms 
%     tempChance = cell2mat(chance_fiscL1O_MargPos(:,:,c));
%     for t = 1:length(mlb.likeTimeVect)
%         for ndx = 1:mlb.seqLength
%             fiscChancePosAcc(t,ndx,c) = mean(tempChance(t,:)==ndx);
%         end
%     end
% end
% clear chance_fiscL1O_MargPos 
% 
% %% Plot FISC-L1O Observed Data Relative to Chance 
% figure; 
% piNdx = find(mlb.likeTimeVect==0)+0.5;
% poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
% rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
% posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
% 
% fiscPosts = mean(cell2mat(fiscL1O_Posts),3);
% 
% fiscTP = cell2mat(fiscL1O_TimePosDecode);
% fiscTPacc = nan(size(mlb.likeTimeVect,1));
% for t = 1:size(fiscTP,1)
%     curTime = fiscTP(t,:);
%     for ndx = 1:size(mlb.likeTimeVect,1)
%         fiscTPacc(t,ndx) = mean(curTime==ndx);
%     end
% end
% sp(1) = subplot(3,3,[2:3, 5:6]);
% z_fiscTPacc = (fiscTPacc-mean(fiscChanceTimePosAcc,3))./std(fiscChanceTimePosAcc,0,3);
% z_fiscTPacc(z_fiscTPacc<2 & z_fiscTPacc>-2) = 0;
% imagesc(z_fiscTPacc, [-2 2]);
% hold on;
% for ndx = 1:length(piNdx)
%     plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
% 
%     if ndx<length(piNdx)
%         plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
%     end
% end
% title('Time-In-Position Decoding Accuracy');
% set(sp(1), 'xticklabel', [], 'yticklabel', []);
% 
% fiscTime = cell2mat(fiscL1O_MargDecodeTime);
% timeVect = unique(mlb.decodeIDvects(:,1));
% fiscTimeAcc = nan(size(fiscTime,1), length(timeVect));
% for t = 1:size(fiscTime,1)
%     curTime = fiscTime(t,:);
%     for ndx = 1:length(timeVect)
%         fiscTimeAcc(t,ndx) = mean(curTime==timeVect(ndx));
%     end
% end
% sp(2) = subplot(3,3,[1,4]);
% z_fiscTimeAcc = (fiscTimeAcc-mean(fiscChanceTimeAcc,3))./std(fiscChanceTimeAcc,0,3);
% z_fiscTimeAcc(z_fiscTimeAcc<2 & z_fiscTimeAcc>-2) = 0;
% imagesc(mlb.obsvTimeVect,1:length(mlb.likeTimeVect), z_fiscTimeAcc, [-2 2]);
% hold on;
% plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
% for ndx = 1:length(piNdx)
%     plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
% 
%     if ndx<length(piNdx)
%         plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
%     end
% end
% title('Temporal Decoding Accuracy');
% set(sp(2), 'xtick', [0, nearestPOtime, nearestRWDtime], 'xticklabel', [{'PI'}, {'PO'}, {'RWD'}], 'xticklabelrotation', 45, 'yticklabel', [], 'TickDir', 'out');
% cb = colorbar;
% cb.Location = 'westoutside';
% cb.Label.String = [{'Z-Score'};{'Norm to Chance'}];
% cb.Label.Position(1) = 0;
% cb.Label.FontWeight = 'Bold';
% cb.Ticks = [-2 2];
% 
% fiscPos = cell2mat(fiscL1O_MargDecodePos);
% fiscPosAcc = nan(size(fiscPos,1), mlb.seqLength);
% for t = 1:size(fiscPos,1)
%     curTime = fiscPos(t,:);
%     for op = 1:mlb.seqLength
%         fiscPosAcc(t,op) = mean(curTime==op);
%     end
% end
% sp(3) = subplot(3,3,8:9);
% cpID = nan(1,mlb.seqLength);
% for op = 1:mlb.seqLength
%     cpID(op) = plot((fiscPosAcc(:,op)-mean(fiscChancePosAcc(:,op,:),3))./std(fiscChancePosAcc(:,op,:),0,3), 'color', mlb.PositionColors(op,:), 'linewidth', 2);
%     hold on;
%     plot(repmat(piNdx(op), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
%     plot(repmat(poNdx(op), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
%     plot(repmat(rwdNdx(op), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
%     if op<length(piNdx)
%         plot(repmat(posNdx(op), [1,2]), get(gca, 'ylim'), '-k', 'linewidth',2);
%     end
% end
% axis tight;
% plot(sp(3), get(gca, 'xlim'), [2 2], ':k', 'linewidth', 2);
% plot(sp(3), get(gca, 'xlim'), [-2 -2], ':k', 'linewidth', 2);
% title('Odor Decoding Accuracy');
% set(sp(3), 'xticklabel', [], 'ytick', [-10 -5 -2 0 2 5 10]);
% sp(3).TickDir = 'out';
% sp(3).XRuler.FirstCrossoverValue = 0;
% sp(3).XRuler.SecondCrossoverValue = 0;
% legend(sp(3), cpID, arrayfun(@(a){sprintf('Pos%i', a)}, 1:mlb.seqLength), 'location', 'north', 'NumColumns', mlb.seqLength);
% box 'off'
% ylabel([{'Z-Score'};{'Norm to Chance'}]);
% linkaxes(sp([1,3]), 'x');
% 
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', sprintf("FISC Leave-One-Out Decode FISC: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
%     binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
% %% Process ISC Chance Data
% % FISC ISC Posterior Values
% iscChancePost = cell(mlb.seqLength, numChancePerms );
% for c = 1:numChancePerms 
%     for o = 1:mlb.seqLength
%         iscChancePost{o,c} = mean(cell2mat(chance_fiscISC_Post(o,c,:)),3);
%     end
% end
% clear chance_fiscISC_Post
% 
% % FISC ISC Time-Position Coding
% iscChanceTimePosAcc = nan(length(mlb.likeTimeVect), length(mlb.likeTimeVect), numChancePerms );
% for c = 1:numChancePerms 
%     curPerm = chance_fiscISC_TimePos(:,:,c);
%     tempAcc = cell(mlb.seqLength,1);
%     for op = 1:mlb.seqLength
%         tempChance = cell2mat(curPerm(op,:));
%         tempTPacc = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
%         for t = 1:length(mlb.obsvTimeVect)
%             for ndx = 1:length(mlb.likeTimeVect)
%                 tempTPacc(t,ndx) = mean(tempChance(t,:)==ndx);
%             end
%         end
%         tempAcc{op} = tempTPacc;
%     end
%     iscChanceTimePosAcc(:,:,c) = cell2mat(tempAcc);
% end
% clear chance_fiscISC_TimePos
% 
% % FISC ISC Time Coding
% iscChanceTimeAcc = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), numChancePerms );
% for c = 1:numChancePerms 
%     curPerm = chance_fiscISC_MargTime(:,:,c);
%     tempAcc = cell(mlb.seqLength,1);
%     for op = 1:mlb.seqLength
%         tempChance = cell2mat(curPerm(op,:));
%         tempTacc = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
%         for t = 1:length(mlb.obsvTimeVect)
%             for ndx = 1:length(mlb.obsvTimeVect)
%                 tempTacc(t,ndx) = mean(tempChance(t,:)==mlb.obsvTimeVect(ndx));
%             end
%         end
%         tempAcc{op} = tempTacc;
%     end
%     iscChanceTimeAcc(:,:,c) = cell2mat(tempAcc);
% end
% clear chance_fiscISC_MargTime
% 
% % FISC ISC Odor Coding
% iscChancePosAcc = nan(length(mlb.likeTimeVect), mlb.seqLength, numChancePerms );
% for c = 1:numChancePerms 
%     curPerm = chance_fiscISC_MargPos(:,:,c);
%     tempAcc = cell(mlb.seqLength,1);
%     for op = 1:mlb.seqLength
%         tempChance = cell2mat(curPerm(op,:));
%         tempPacc = nan(length(mlb.obsvTimeVect), mlb.seqLength);
%         for t = 1:length(mlb.obsvTimeVect)
%             for ndx = 1:mlb.seqLength
%                 tempPacc(t,ndx) = mean(tempChance(t,:)==ndx);
%             end
%         end
%         tempAcc{op} = tempPacc;
%     end
%     iscChancePosAcc(:,:,c) = cell2mat(tempAcc);
% end
% clear chance_fiscISC_MargPos 
% 
% %% Plot FISC-ISC Observed Data Relative to Chance 
% figure; 
% piNdx = find(mlb.likeTimeVect==0)+0.5;
% poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
% rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
% posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
% 
% iscPosts = cell(mlb.seqLength,1);
% for pos = 1:mlb.seqLength
%     iscPosts{pos} = mean(cell2mat(fiscISC_Posts(pos,:,:)),3);
% end
% 
% iscTPacc = cell(mlb.seqLength,1);
% for op = 1:mlb.seqLength
%     curISCtp = cell2mat(fiscISC_TimePosDecode(op,:));
%     tempAcc = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
%     for t = 1:length(mlb.obsvTimeVect)
%         for ndx = 1:length(mlb.likeTimeVect)
%             tempAcc(t,ndx) = mean(curISCtp(t,:)==ndx);
%         end
%     end
%     iscTPacc{op} = tempAcc;
% end
% iscTPacc = cell2mat(iscTPacc);
% sp(1) = subplot(3,3,[2:3 5:6]);
% z_iscTPacc = (iscTPacc-mean(iscChanceTimePosAcc,3))./std(iscChanceTimePosAcc,0,3);
% z_iscTPacc(z_iscTPacc<2 & z_iscTPacc>-2) = 0;
% imagesc(z_iscTPacc, [-2 2]);
% hold on; 
% for ndx = 1:length(piNdx)
%     plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
% 
%     if ndx<length(piNdx)
%         plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
%     end
% end
% title('Time-In-Position Decoding Accuracy');
% set(sp(1), 'xticklabel', [], 'yticklabel', []);
% 
% iscTimeAcc = cell(mlb.seqLength, 1);
% for op = 1:mlb.seqLength
%     curISCt = cell2mat(fiscISC_MargDecodeTime(op,:));
%     tempAcc = nan(length(mlb.obsvTimeVect));
%     for t = 1:length(mlb.obsvTimeVect)
%         for ndx = 1:length(mlb.obsvTimeVect)
%             tempAcc(t,ndx) = mean(curISCt(t,:)==mlb.obsvTimeVect(ndx));
%         end
%     end
%     iscTimeAcc{op} = tempAcc;
% end
% iscTimeAcc = cell2mat(iscTimeAcc);
% sp(2) = subplot(3,3,[1,4]);
% z_iscTimeAcc = (iscTimeAcc-mean(iscChanceTimeAcc,3))./std(iscChanceTimeAcc,0,3);
% z_iscTimeAcc(z_iscTimeAcc<2 & z_iscTimeAcc>-2) = 0;
% imagesc(mlb.obsvTimeVect, 1:length(mlb.likeTimeVect), z_iscTimeAcc, [-2 2]);
% hold on;
% plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
% for ndx = 1:length(piNdx)
%     plot(get(gca, 'xlim'),repmat(piNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(poNdx(ndx), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(rwdNdx(ndx), [1,2]), ':k','linewidth', 2);
% 
%     if ndx<length(piNdx)
%         plot(get(gca, 'xlim'),repmat(posNdx(ndx), [1,2]), '-k','linewidth', 2);
%     end
% end
% title('Temporal Decoding Accuracy');
% set(sp(2), 'xtick', [0, nearestPOtime, nearestRWDtime], 'xticklabel', [{'PI'}, {'PO'}, {'RWD'}], 'xticklabelrotation', 45, 'yticklabel', [], 'TickDir', 'out');
% cb = colorbar;
% cb.Location = 'westoutside';
% cb.Label.String = [{'Z-Score'};{'Norm to Chance'}];
% cb.Label.Position(1) = 0;
% cb.Label.FontWeight = 'Bold';
% cb.Ticks = [-2 2];
% 
% iscPosAcc = cell(mlb.seqLength,1);
% for op = 1:mlb.seqLength
%     curISCp = cell2mat(fiscISC_MargDecodePos(op,:));
%     tempAcc = nan(length(mlb.obsvTimeVect),mlb.seqLength);
%     for t = 1:length(mlb.obsvTimeVect)
%         for pos = 1:mlb.seqLength
%             tempAcc(t,pos) = mean(curISCp(t,:)==pos);
%         end
%     end
%     iscPosAcc{op} = tempAcc;
% end
% iscPosAcc = cell2mat(iscPosAcc);
% sp(3) = subplot(3,3,8:9);
% cpID = nan(1,mlb.seqLength);
% for op = 1:mlb.seqLength
%     cpID(op) = plot((iscPosAcc(:,op)-mean(iscChancePosAcc(:,op,:),3))./std(iscChancePosAcc(:,op,:),0,3), 'color', mlb.PositionColors(op,:), 'linewidth', 2);
%     hold on;
%     plot(repmat(piNdx(op), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
%     plot(repmat(poNdx(op), [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
%     plot(repmat(rwdNdx(op), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
%     if op<length(piNdx)
%         plot(repmat(posNdx(op), [1,2]), get(gca, 'ylim'), '-k', 'linewidth',2);
%     end
% end
% axis tight;
% plot(sp(3), get(gca, 'xlim'), [2 2], ':k', 'linewidth', 2);
% plot(sp(3), get(gca, 'xlim'), [-2 -2], ':k', 'linewidth', 2);
% title('Odor Decoding Accuracy');
% set(sp(3), 'xticklabel', [], 'ytick', [-10 -5 -2 0 2 5 10]);
% sp(3).TickDir = 'out';
% sp(3).XRuler.FirstCrossoverValue = 0;
% sp(3).XRuler.SecondCrossoverValue = 0;
% legend(sp(3), cpID, arrayfun(@(a){sprintf('Pos%i', a)}, 1:mlb.seqLength), 'location', 'north', 'NumColumns', mlb.seqLength);
% box 'off'
% ylabel([{'Z-Score'};{'Norm to Chance'}]);
% linkaxes(sp([1,3]), 'x');
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', sprintf("FISC Decode ISC: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
%     binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
