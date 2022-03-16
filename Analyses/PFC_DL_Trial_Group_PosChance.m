% PFC_DL_Trial_Group_PosChance
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
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%%
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb = MLB_SM(fileDirs{ani});
    % Create Analysis Variables
    if ani == 1 
        % Behavior Variables
        fiscPokeOutLat = cell(length(fileDirs),1);
        fiscRwdDelivLat = cell(length(fileDirs),1);
        smi = nan(length(fileDirs),size(mlb.odrSeqs,1));
        dPrm = nan(length(fileDirs),size(mlb.odrSeqs,1));
        ri = nan(length(fileDirs),size(mlb.odrSeqs,1));
        smiByOP = nan(length(fileDirs),numel(mlb.odrSeqs),2);
        dPrmByOP = nan(length(fileDirs),numel(mlb.odrSeqs),2);
        riByOP = nan(length(fileDirs),numel(mlb.odrSeqs),2);
        % Posteriors
        realPost = cell(numel(mlb.odrSeqs), 1, length(fileDirs));
        chancePost = cell(numel(mlb.odrSeqs), numChancePerms, length(fileDirs));
    end
        
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
    smi(ani,:) = mlb.smi;
    dPrm(ani,:) = mlb.dPrime;
    ri(ani,:) = mlb.ri;
    for op = 1:mlb.seqLength
        smiByOP(ani,:,1) = reshape(mlb.smiByPos', [numel(mlb.smiByPos),1]);
        smiByOP(ani,:,2) = reshape(mlb.smiByOdr', [numel(mlb.smiByPos),1]);
        dPrmByOP(ani,:,1) = reshape(mlb.dPrimeByPos', [numel(mlb.smiByPos),1]);
        dPrmByOP(ani,:,2) = reshape(mlb.dPrimeByOdr', [numel(mlb.smiByPos),1]);
        riByOP(ani,:,1) = reshape(mlb.riByPos', [numel(mlb.smiByPos),1]);
        riByOP(ani,:,2) = reshape(mlb.riByOdr', [numel(mlb.smiByPos),1]);
    end
    %% Process Observations
    % FISC via Leave-1-Out
    mlb.SetLikes_ISC;
    mlb.Process_LikelyL1O;
    realPost(:,1,ani) = mlb.post;
    %% Process Chance
    for perm = 1:numChancePerms 
        chancePerm = sortrows([randperm(numel(mlb.postTrlIDs(~isnan(mlb.postTrlIDs))))',find(~isnan(mlb.postTrlIDs))]);
        mlb.likeTrlSpikes(mlb.postTrlIDs(~isnan(mlb.postTrlIDs))) = mlb.likeTrlSpikes(chancePerm(:,2));
        mlb.Process_LikelyL1O;
        chancePost(:,perm,ani) = mlb.post;
    end
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%% Collapse Trials Across Animals
groupPostTIP_Real = cell(length(mlb.odrVect),1);            groupPostTIP_Chance = cell(length(mlb.odrVect),2);
groupPostOdr_Real = cell(length(mlb.odrVect),1);            groupPostOdr_Chance = cell(length(mlb.odrVect),2);
groupPostPos_Real = cell(length(mlb.odrVect),1);            groupPostPos_Chance = cell(length(mlb.odrVect),2);
groupPostSeqDiff_Real = cell(length(mlb.odrVect),1);        groupPostSeqDiff_Chance = cell(length(mlb.odrVect),2);
for odr = 1:length(mlb.odrVect)
    %% Real Data
    groupPostTIP_Real{odr} = cell2mat(realPost(odr,:,:));
    groupPostOdr_Real{odr} = mlb.TabulateBayesPost(groupPostTIP_Real{odr}, mlb.decodeIDvects(:,4));
    groupPostPos_Real{odr} = mlb.TabulateBayesPost(groupPostTIP_Real{odr}, mlb.decodeIDvects(:,3));
    tempSeqDiff = nan(length(mlb.obsvTimeVect),mlb.seqLength, size(groupPostOdr_Real{odr},3));
    for pos = 1:mlb.seqLength
        tempSeqDiff(:,pos,:) = diff(groupPostOdr_Real{odr}(:,[pos, pos+mlb.seqLength],:),1,2);
    end
    groupPostSeqDiff_Real{odr} = tempSeqDiff;
    %% Chance
    tempPostTIP_Chance = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
    tempPostPos_Chance = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    tempPostOdr_Chance = nan(length(mlb.obsvTimeVect), length(mlb.odrVect), numChancePerms);
    tempPostSeqDiff_Chance = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    for perm = 1:numChancePerms
        curChanceTIP = cell2mat(chancePost(odr,perm,:));
        tempPostTIP_Chance(:,:,perm) = mean(curChanceTIP,3, 'omitnan');
        tempPostPos_Chance(:,:,perm) = mean(mlb.TabulateBayesPost(curChanceTIP, mlb.decodeIDvects(:,3)),3, 'omitnan');
        tempOdorTab = mlb.TabulateBayesPost(curChanceTIP, mlb.decodeIDvects(:,4));
        tempPostOdr_Chance(:,:,perm) = mean(tempOdorTab,3, 'omitnan');
        tempSeqDiff = nan(length(mlb.obsvTimeVect),mlb.seqLength, size(groupPostOdr_Real{odr},3));
        for pos = 1:mlb.seqLength
            tempSeqDiff(:,pos,:) = diff(tempOdorTab(:,[pos, pos+mlb.seqLength],:),1,2);
        end
        tempPostSeqDiff_Chance(:,:,perm) = mean(tempSeqDiff,3,'omitnan');
    end
    groupPostTIP_Chance{odr,1} = mean(tempPostTIP_Chance,3,'omitnan');
    groupPostTIP_Chance{odr,2} = std(tempPostTIP_Chance,0,3,'omitnan');
    
    groupPostOdr_Chance{odr,1} = mean(tempPostOdr_Chance,3,'omitnan');
    groupPostOdr_Chance{odr,2} = std(tempPostOdr_Chance,0,3,'omitnan');
    
    groupPostPos_Chance{odr,1} = mean(tempPostPos_Chance,3,'omitnan');
    groupPostPos_Chance{odr,2} = std(tempPostPos_Chance,0,3,'omitnan');
    
    groupPostSeqDiff_Chance{odr,1} = mean(tempPostSeqDiff_Chance,3,'omitnan');
    groupPostSeqDiff_Chance{odr,2} = std(tempPostSeqDiff_Chance,0,3,'omitnan');
end

%% Plot things
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);

piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;

figure;
colormap(cMap);
z = colormap;
for odr = 1:length(mlb.odrVect)     
    figure
    subplot(length(mlb.odrVect)+2,length(mlb.odrVect),sub2ind([length(mlb.odrVect),length(mlb.odrVect)+2],1:length(mlb.odrVect),repmat(odr, [1,length(mlb.odrVect)])));
    imagesc(mean(groupPostTIP_Real{odr},3, 'omitnan'),[0 0.05]);
    hold on;
    plot(get(gca, 'xlim'),repmat(piNdx(1), [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(poNdx(1), [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(rwdNdx(1), [1,2]), ':k','linewidth', 2);
    for ndx = 1:length(piNdx)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
    if odr==length(mlb.odrVect)
        cb = colorbar;
        cb.Location = 'westoutside';
        cb.Label.String = 'P(T,Odor|Spks)';
        cb.Label.Position(1) = 0;
        cb.Label.FontWeight = 'Bold';
        cb.Ticks = [0 0.05];
    end
    set(gca, 'xticklabel', [], 'yticklabel', []);
    title(sprintf('Time in Odor Info Odor %i', mlb.odrVect(odr)));
    tempTIPReal = mean(groupPostTIP_Real{odr},3, 'omitnan')-(tinv(0.975,sum(~isnan(groupPostTIP_Real{odr}(1,1,:)))-1).*mlb.SEMcalc(groupPostTIP_Real{odr},0,3));
    tempTIPThresh = groupPostTIP_Chance{odr,1}+(tinv(0.975,numChancePerms-1).*(groupPostTIP_Chance{odr,1}./sqrt(numChancePerms-1)));
    abvThresh = tempTIPReal>tempTIPThresh;
    bounds = bwboundaries(abvThresh);
    for b = 1:length(bounds)
        plot(bounds{b}(:,2), bounds{b}(:,1), 'k', 'linewidth', 2);
    end
    
    subplot(length(mlb.odrVect)+2,length(mlb.odrVect),sub2ind([length(mlb.odrVect),length(mlb.odrVect)+2],odr,length(mlb.odrVect)+1));
    tempChanceMean = groupPostOdr_Chance{odr,1}(:,odr);
    tempChanceSEM = groupPostOdr_Chance{odr,2}(:,odr)./sqrt(numChancePerms-1);
    tempChanceCI = tinv(0.975, numChancePerms-1).*tempChanceSEM;
    tempReal = groupPostOdr_Real{odr};
    tempRealMean = mean(tempReal,3, 'omitnan');
    tempRealSEM = mlb.SEMcalc(tempReal,0,3);
    tempRealCI = tinv(0.975, sum(~isnan(tempReal(1,1,:)),3)-1).*tempRealSEM;
    for o = 1:length(mlb.odrVect)
        [seq,p] = find(mlb.odrSeqs==mlb.odrVect(o));
        if seq==1
            plot(mlb.obsvTimeVect, tempRealMean(:,o), 'color', mlb.PositionColors(p,:), 'linewidth', 1.5);
        else
            plot(mlb.obsvTimeVect, tempRealMean(:,o), 'linestyle', '--', 'color', mlb.PositionColors(p,:), 'linewidth', 1.5);
        end
        hold on;
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempRealMean(:,o)+tempRealSEM(:,o))', flipud(tempRealMean(:,o)-tempRealSEM(:,o))'],...
            'linestyle', 'none', 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0.25);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempRealMean(:,o)+tempRealCI(:,o))', flipud(tempRealMean(:,o)-tempRealCI(:,o))'],...
            'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
        
        plot(mlb.obsvTimeVect, tempChanceMean, 'color', 'k', 'linewidth', 1.5);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempChanceMean+tempChanceSEM)', flipud(tempChanceMean-tempChanceSEM)'],...
            'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempChanceMean+tempChanceCI)', flipud(tempChanceMean-tempChanceCI)'],...
            'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
        set(gca, 'ylim', [0 1]);
        
        plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        if odr==1
            sortedYticks = sortrows([{0};{tempChanceMean(1)}; {0.5};{1}]);
            sortedYtickLabels = sortedYticks;
            sortedYtickLabels{cellfun(@(a)a==tempChanceMean(1),sortedYticks)} = 'Chance';
            set(gca, 'ytick', cell2mat(sortedYticks), 'yticklabel', sortedYtickLabels);
            ylabel('P(Odr|Spks,Pos)');
        else
            set(gca, 'yticklabel', []);
        end
    end
    title(sprintf('Odr %i', odr));
    
    subplot(length(mlb.odrVect)+2,length(mlb.odrVect),sub2ind([length(mlb.odrVect),length(mlb.odrVect)+2],odr,length(mlb.odrVect)+2));
    [seq,pos] = find(mlb.odrSeqs==mlb.odrVect(odr));
    tempChanceMean = groupPostPos_Chance{odr,1}(:,pos);
    tempChanceSEM = groupPostPos_Chance{odr,2}(:,pos)./sqrt(numChancePerms-1);
    tempChanceCI = tinv(0.975, numChancePerms-1).*tempChanceSEM;
    tempReal = groupPostPos_Real{odr};
    tempRealMean = mean(tempReal,3, 'omitnan');
    tempRealSEM = mlb.SEMcalc(tempReal,0,3);
    tempRealCI = tinv(0.975, sum(~isnan(tempReal(1,1,:)),3)-1).*tempRealSEM;
    for p = 1:mlb.seqLength
        plot(mlb.obsvTimeVect, tempRealMean(:,p), 'color', mlb.PositionColors(p,:), 'linewidth', 1.5);
        hold on;
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempRealMean(:,p)+tempRealSEM(:,p))', flipud(tempRealMean(:,p)-tempRealSEM(:,p))'],...
            'linestyle', 'none', 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0.25);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempRealMean(:,p)+tempRealCI(:,p))', flipud(tempRealMean(:,p)-tempRealCI(:,p))'],...
            'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
        
        plot(mlb.obsvTimeVect, tempChanceMean, 'color', 'k', 'linewidth', 1.5);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempChanceMean+tempChanceSEM)', flipud(tempChanceMean-tempChanceSEM)'],...
            'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(tempChanceMean+tempChanceCI)', flipud(tempChanceMean-tempChanceCI)'],...
            'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
        set(gca, 'ylim', [0 1]);
        
        plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        if odr==1
            sortedYticks = sortrows([{0};{tempChanceMean(1)}; {0.5};{1}]);
            sortedYtickLabels = sortedYticks;
            sortedYtickLabels{cellfun(@(a)a==tempChanceMean(1),sortedYticks)} = 'Chance';
            set(gca, 'ytick', cell2mat(sortedYticks), 'yticklabel', sortedYtickLabels);
            ylabel('P(Pos|Spks,Pos)');
        else
            set(gca, 'yticklabel', []);
        end
    end
    title(sprintf('Pos %i', odr));
end
annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("Dual-List ISC-L1O; binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot Odor Difference
figure;
allTrls_Seq1 = cell(1,1,mlb.seqLength);
allTrls_Seq2 = cell(1,1,mlb.seqLength);
allTrls_Chance = cell(1,3);
for pos = 1:mlb.seqLength
    h = subplot(5,4,sub2ind([mlb.seqLength,mlb.seqLength+1],pos,mlb.seqLength+1));
    hold on;
    tempReal_Seq1 = groupPostSeqDiff_Real{pos}(:,pos,:);
    seq1Mean = mean(tempReal_Seq1,3,'omitnan');
    seq1SEM = mlb.SEMcalc(tempReal_Seq1,0,3);
    seq1CI = tinv(0.975, sum(~isnan(tempReal_Seq1(1,1,:)),3)-1).*seq1SEM;
    plot(mlb.obsvTimeVect, seq1Mean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(seq1Mean+seq1SEM)', flipud(seq1Mean-seq1SEM)'],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(seq1Mean+seq1CI)', flipud(seq1Mean-seq1CI)'],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    allTrls_Seq1{pos} = tempReal_Seq1;
    
    tempReal_Seq2 = groupPostSeqDiff_Real{pos+mlb.seqLength}(:,pos,:);
    seq2Mean = mean(tempReal_Seq2,3,'omitnan');
    seq2SEM = mlb.SEMcalc(tempReal_Seq2,0,3);
    seq2CI = tinv(0.975, sum(~isnan(tempReal_Seq2(1,1,:)),3)-1).*seq2SEM;
    plot(mlb.obsvTimeVect, seq2Mean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5, 'linestyle', '--');
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(seq2Mean+seq2SEM)', flipud(seq2Mean-seq2SEM)'],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(seq2Mean+seq2CI)', flipud(seq2Mean-seq2CI)'],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    allTrls_Seq2{pos} = tempReal_Seq2;
    
    tempChanceMean = groupPostSeqDiff_Chance{pos,1}(:,pos);
    tempChanceSEM = groupPostSeqDiff_Chance{pos,2}(:,pos)./sqrt(numChancePerms-1);
    tempChanceCI = tinv(0.975, numChancePerms-1).*tempChanceSEM;
    plot(mlb.obsvTimeVect, tempChanceMean, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(tempChanceMean+tempChanceSEM)', flipud(tempChanceMean-tempChanceSEM)'],...
        'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(tempChanceMean+tempChanceCI)', flipud(tempChanceMean-tempChanceCI)'],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    allTrls_Chance{pos} = tempReal_Seq2;
    
    plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
    
    
    h.XRuler.FirstCrossoverValue = 0;
    h.XRuler.SecondCrossoverValue = 0;
end
h = subplot(5,4,1:16);
hold on;
seq1 = cell2mat(allTrls_Seq1);
seq1Mean = mean(seq1,3,'omitnan');
seq1SEM = mlb.SEMcalc(tempReal_Seq1,0,3);
seq1CI = tinv(0.975, sum(~isnan(tempReal_Seq1(1,1,:)),3)-1).*seq1SEM;
plot(mlb.obsvTimeVect, seq1Mean, 'color', 'k', 'linewidth', 1.5);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(seq1Mean+seq1SEM)', flipud(seq1Mean-seq1SEM)'],...
    'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(seq1Mean+seq1CI)', flipud(seq1Mean-seq1CI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0);

seq2 = cell2mat(allTrls_Seq2);
seq2Mean = mean(seq2,3,'omitnan');
seq2SEM = mlb.SEMcalc(tempReal_Seq2,0,3);
seq2CI = tinv(0.975, sum(~isnan(tempReal_Seq2(1,1,:)),3)-1).*seq2SEM;
plot(mlb.obsvTimeVect, seq2Mean, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(seq2Mean+seq2SEM)', flipud(seq2Mean-seq2SEM)'],...
    'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(seq2Mean+seq2CI)', flipud(seq2Mean-seq2CI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0);
allTrls_Seq2{pos} = tempReal_Seq2;

plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);

h.XRuler.FirstCrossoverValue = 0;
h.XRuler.SecondCrossoverValue = 0;
linkaxes
    
%%
clear chancePost
save('PFC_DualList_PosChance.mat', '-v7.3');