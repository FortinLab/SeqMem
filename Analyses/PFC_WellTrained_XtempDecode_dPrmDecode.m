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
% tets = [1,22,17,18,17]; % Lateral/Distal
% % tets = [7,3,1,5,5]; % Medial/Proximal
% setupSeqLength = 5;

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
setupSeqLength = 4; 

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
ssProportion = 0.5;
numPerms = 100;
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
grpDprmDecodePosts = cell(1,length(fileDirs));
grpDprmDecodePostsChance = cell(1,length(fileDirs));
grpDecodes = cell(1,length(fileDirs));
grpDecodeTrlNums = cell(1,length(fileDirs));
grpDecodeTrlPosIDs = cell(1,length(fileDirs));
grpDecodePermIDs = cell(1,length(fileDirs));
grpLFP = cell(1,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    nsmblType = 'All Cells';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
% %     betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; nsmblType = 'Mod Only';                      % only MODULATED cells
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; nsmblType = 'NonMod Only';                      % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.numPerms = 1;
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
    [~, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
    [~, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
    grpLFP{ani} = [betaPower, thetaPower];
%     [~, gammaPower] = mlb.PP_TrialMatrix_LFP([40 100], trlWindow{1}, alignment{1});
%     grpLFP{ani} = [betaPower, thetaPower, gammaPower];
    
    dPrmPerm = cell(1,numPerms);
    dPrmPermChance = cell(1,numPerms);
    tempDecodes = cell(1,1,numPerms);
    tempTrlNums = cell(1,numPerms);
    tempTrlIDs = cell(1,numPerms);
    tempPermIDs = cell(1,numPerms);
    for perm = 1:numPerms
        fprintf('Iteration #%i', perm);
        mlb.SetLikes_SubSample;
        mlb.Process_IterativeObserves;
        [decodes, ~] = mlb.DecodeBayesPost(mlb.post{1});
        trlPosVect = [mlb.trialInfo(mlb.obsvTrlIDs{1}).Position];
        [dPrmPerm{perm},dPrmPermChance{perm}] = mlb.CalcDprmMtxFromDecode(decodes, trlPosVect);
        tempDecodes{perm} = decodes;
        tempTrlNums{perm} = squeeze(mlb.postTrlIDs{1})';
        tempTrlIDs{perm} = [mlb.trialInfo(mlb.postTrlIDs{1}).Position];
        tempPermIDs{perm} = squeeze(ones(size(mlb.postTrlIDs{1})))'*perm;
        fprintf(' complete\n');
    end
    grpDprmDecodePosts{ani} = dPrmPerm;
    grpDprmDecodePostsChance{ani} = dPrmPermChance;
    grpDecodes{ani}	= cell2mat(tempDecodes);
    grpDecodeTrlNums{ani} = cell2mat(tempTrlNums);
    grpDecodeTrlPosIDs{ani} = cell2mat(tempTrlIDs);
    grpDecodePermIDs{ani} = cell2mat(tempPermIDs);
    
    %%
    figure; 
    for t = 1:mlb.seqLength
        subplot(ceil(mlb.seqLength/2),floor(mlb.seqLength/2),t); 
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(permute(cellfun(@(a){a(:,:,t)},dPrmPerm), [3,1,2])),3)); 
        colorbar;
        set(gca, 'ydir', 'normal'); 
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(mean(fiscPokeOutLat{ani}), [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(mean(mean(fiscPokeOutLat{ani})), [1,2]), '-k');
        plot(repmat(mean(fiscRwdDelivLat{ani}), [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(mean(mean(fiscRwdDelivLat{ani})), [1,2]), ':k', 'linewidth', 2);
        title(t);
    end
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', mlb.pathDir,...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    drawnow;
end    

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

%% Plot Group Data
dPrmPermByOdor = cell(1,mlb.seqLength, length(fileDirs));
dPrmPermByOdorChance = cell(1,mlb.seqLength, length(fileDirs));
dPrmPermByOdorTstat = cell(1,mlb.seqLength, length(fileDirs));
dPrmMeansByOdorDiff = cell(1,mlb.seqLength, length(fileDirs));
dPrmMeansByOdor = cell(1,mlb.seqLength, length(fileDirs));
for ani = 1:length(fileDirs)
    curAniPerms = grpDprmDecodePosts{ani};
    curAniChance = grpDprmDecodePostsChance{ani};
    for op = 1:mlb.seqLength
        dPrmPermByOdor{1,op,ani} = cell2mat(permute(cellfun(@(a){a(:,:,op)}, curAniPerms), [1,3,2]));
        dPrmPermByOdorChance{1,op,ani} = cell2mat(permute(cellfun(@(a){a(:,:,op)}, curAniChance), [1,3,2]));
        tempTstatMtx = nan(size(curAniPerms{1},1), size(curAniPerms{1},2));
        for t1 = 1:size(tempTstatMtx,1)
            for t2 = 1:size(tempTstatMtx,2)
                [~,~,~,stats] = ttest2(squeeze(dPrmPermByOdor{1,op,ani}(t1,t2,:)), squeeze(dPrmPermByOdorChance{1,op,ani}(t1,t2,:)));
                tempTstatMtx(t1,t2) = stats.tstat;
            end
        end
        dPrmMeansByOdorDiff{1,op,ani} = mean(dPrmPermByOdor{1,op,ani},3) - mean(dPrmPermByOdorChance{1,op,ani},3);
        dPrmMeansByOdor{1,op,ani} = mean(dPrmPermByOdor{1,op,ani},3);
        dPrmPermByOdorTstat{1,op,ani} = tempTstatMtx;
    end
end
figure;
odrMnVar = nan(mlb.seqLength,2);
for t = 1:mlb.seqLength    
    subplot(ceil(mlb.seqLength/2),floor(mlb.seqLength/2),t);
    curDprmMtx = mean(cell2mat(dPrmMeansByOdor(1,t,:)),3);
%     curDprmMtx = mean(cell2mat(dPrmMeansByOdorDiff(1,t,:)),3);
    curDprmTMtx = mean(cell2mat(dPrmPermByOdorTstat(1,t,:)),3);
    curDprmTMtx(curDprmTMtx<3) = 0;
    curDprmTMtx(curDprmTMtx>0) = 1;
    odrMnVar(t,1) = mean(curDprmMtx(:));
    odrMnVar(t,2) = std(curDprmMtx(:));
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, curDprmMtx);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(dPrmPermByOdorTstat(1,t,:)),3));
    hold on;
    bounds = bwboundaries(curDprmTMtx);
    for b = 1:length(bounds)
        if numel(bounds{b})>4
            plot(mlb.obsvTimeVect(bounds{b}(:,2)), mlb.obsvTimeVect(bounds{b}(:,1)), 'k', 'linewidth', 3);
        end
    end
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, zscore(curDprmMtx,0,'all'), [-3 3]);
    colorbar;
    set(gca, 'ydir', 'normal');
    hold on;
    plot([0 0], get(gca, 'ylim'), '-k');
    plot(get(gca, 'xlim'), [0 0], '-k');
    plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '-k');
    plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '-k');
    plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
    title(t);
end

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf('%s aligned; %s Ensemble, Trial Window = (%.0fms:%.0fms); %i Perms; %ims Window', alignment{1}, nsmblType, trlWindow{1}(1), trlWindow{1}(2), numPerms, binSize),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Process LFP Power
% First process individual animal beta power
if strcmp(alignment{1}, 'PokeIn')
    trlTimeLog = mlb.obsvTimeVect>0;
%     trlTimeLog = mlb.obsvTimeVect>nearestPOtime-0.5 & mlb.obsvTimeVect<nearestRWDtime;
else
    trlTimeLog = mlb.obsvTimeVect>-0.5;
end
aniTpwr = cell(length(fileDirs), size(grpLFP{1},2));
for ani = 1:length(fileDirs)
    aniLFP = grpLFP{ani};
    aniTrlNums = grpDecodeTrlNums{ani};
    trlMeans = squeeze(mean(aniLFP(trlTimeLog,:,aniTrlNums)));
    for band = 1:size(trlMeans,1)
        aniTpwr{ani,band} = trlMeans(band,:);
    end
end
for band = 1:size(aniTpwr,2)
    if band ==1
        bnm = 'Beta';
    elseif band == 2
        bnm = 'Theta';
    elseif band == 3
        bnm = 'Gamma';
    end
    tempPwr = cell2mat(aniTpwr(:,band)');
    sortedPower = sort(tempPwr);
    lfpThresh = [sortedPower(ceil(length(sortedPower)*0.4)), sortedPower(floor(length(sortedPower)*0.6))];
    tempDecodes = cell2mat(permute(grpDecodes, [1,3,2]));
    tempTrlIDs = cell2mat(grpDecodeTrlPosIDs);
    tempPermIDs = cell2mat(grpDecodePermIDs);
    hptDprmReal = cell(1,numPerms);
    lptDprmReal = cell(1,numPerms);
    hptDprmChance = cell(1,numPerms);
    lptDprmChance = cell(1,numPerms);
    for perm = 1:numPerms
        [hptDprmReal{perm},hptDprmChance{perm}] = mlb.CalcDprmMtxFromDecode(tempDecodes(:,:,tempPwr>lfpThresh(2) & tempPermIDs==perm), tempTrlIDs(tempPwr>lfpThresh(2) & tempPermIDs==perm));
        [lptDprmReal{perm},lptDprmChance{perm}] = mlb.CalcDprmMtxFromDecode(tempDecodes(:,:,tempPwr<lfpThresh(1) & tempPermIDs==perm), tempTrlIDs(tempPwr<lfpThresh(1) & tempPermIDs==perm));
    end
    figure;
    for op = 1:mlb.seqLength
        tempHPTdReal = cell2mat(permute(cellfun(@(a){a(:,:,op)},hptDprmReal), [1,3,2]));
        tempLPTdReal = cell2mat(permute(cellfun(@(a){a(:,:,op)},lptDprmReal), [1,3,2]));
        tempHPTdChance = cell2mat(permute(cellfun(@(a){a(:,:,op)},hptDprmChance), [1,3,2]));
        tempLPTdChance = cell2mat(permute(cellfun(@(a){a(:,:,op)},lptDprmChance), [1,3,2]));
        tHPT = nan(size(tempHPTdReal,1), size(tempHPTdReal,2));
        tLPT = nan(size(tempLPTdReal,1), size(tempLPTdReal,2));
        for t1 = 1:size(tHPT,1)
            for t2 = 1:size(tHPT,2)
                [~,~,~,hptStats] = ttest2(squeeze(tempHPTdReal(t1,t2,:)), squeeze(tempHPTdChance(t1,t2,:)));
                tHPT(t1,t2) = hptStats.tstat;
                [~,~,~,lptStats] = ttest2(squeeze(tempLPTdReal(t1,t2,:)), squeeze(tempLPTdChance(t1,t2,:)));
                tLPT(t1,t2) = lptStats.tstat;
            end
        end        
        
        subplot(mlb.seqLength,3,sub2ind([3,mlb.seqLength],1,op));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tHPT);
        set(gca, 'ydir', 'normal');
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '-k');
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
        title(sprintf('HPT Pos%i', op));
        colorbar;
        drawnow;
        
        subplot(mlb.seqLength,3,sub2ind([3,mlb.seqLength],2,op));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tLPT);
        set(gca, 'ydir', 'normal');
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '-k');
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
        title(sprintf('LPT Pos%i', op));
        colorbar;
        drawnow;
        
        subplot(mlb.seqLength,3,sub2ind([3,mlb.seqLength],3,op));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tHPT-tLPT);
        set(gca, 'ydir', 'normal');
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(nearestPOtime, [1,2]), '-k');
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(nearestRWDtime, [1,2]), ':k', 'linewidth', 2);
        title(sprintf('HPT-LPT Pos%i', op));        
        colorbar;
        drawnow;
    end
    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s Split; %s Ensemble, %s aligned; Trial Window = (%.0fms:%.0fms); %i Perms', bnm, nsmblType, alignment{1}, trlWindow{1}(1), trlWindow{1}(2),numPerms),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end