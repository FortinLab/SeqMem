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
    
end    

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

%% Plot individual animal dPrime Decodings
for ani = 1:length(fileDirs)
    aniPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscPokeOutLat{ani}),1,'last'));
    aniRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscRwdDelivLat{ani}),1,'last'));
    
    dPrmPerm = grpDprmDecodePosts{ani};
    figure;
    for o = 1:mlb.seqLength
        subplot(ceil(mlb.seqLength/2),floor(mlb.seqLength/2),o);
        imagesc(mlb.obsvTimeVect,mlb.obsvTimeVect,mean(cell2mat(permute(cellfun(@(a){a(:,:,o)},dPrmPerm), [1,3,2])),3));
        colorbar;
        set(gca, 'ydir', 'normal');
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(aniPOtime, [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(aniPOtime, [1,2]), '-k');
        plot(repmat(aniRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(aniRWDtime, [1,2]), ':k', 'linewidth', 2);
        title(o);
    end
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% Model analysis
persModSims = cell(3,length(fileDirs));
dynModSims = cell(3,length(fileDirs));
mixModSims = cell(3,length(fileDirs));
for ani = 1:length(fileDirs)
    dPrmPerm = grpDprmDecodePosts{ani};
    %% First create animal specific epoch logicals 
    aniPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscPokeOutLat{ani}),1,'last'));
    aniRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscRwdDelivLat{ani}),1,'last'));
    
    preTrlLog = mlb.obsvTimeVect<0;
    odrSampLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<aniPOtime;
    preRwdPrd = mlb.obsvTimeVect>=aniPOtime & mlb.obsvTimeVect<aniRWDtime;
    pstTrlPrd = mlb.obsvTimeVect>=aniRWDtime;
    
    %% Next find the best fit dynamic model for the animal's data
    % Create potential fully dynamic models based on 1/2 the width of the pre-trial or trial periods (whichever's shorter)
    dynCap = floor(length(mlb.obsvTimeVect)/2-1);
    fullDynMod = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), dynCap+1);
    fullDynMod(:,:,1) = eye(length(mlb.obsvTimeVect));
    for d = 1:dynCap
        fullDynMod(:,:,d+1) = triu(true(length(mlb.obsvTimeVect)), d*-1) & tril(true(length(mlb.obsvTimeVect)), d);
    end
            
    %% Next identify the best fit persistent model that fits the animal data
    nullModel = true(length(mlb.obsvTimeVect));
    pObound = false(length(mlb.obsvTimeVect)); 
        pObound(preTrlLog|odrSampLog,preTrlLog|odrSampLog)=true; 
        pObound(preRwdPrd|pstTrlPrd,preRwdPrd|pstTrlPrd)=true;
    pIbound = false(length(mlb.obsvTimeVect));
        pIbound(preTrlLog,preTrlLog)=true;
        pIbound(odrSampLog|preRwdPrd|pstTrlPrd,odrSampLog|preRwdPrd|pstTrlPrd) = true;
    pIpObound = false(length(mlb.obsvTimeVect)); 
        pIpObound(preTrlLog,preTrlLog)=true; 
        pIpObound(odrSampLog,odrSampLog)=true;
        pIpObound(preRwdPrd|pstTrlPrd,preRwdPrd|pstTrlPrd)=true;
    pIpOrwDbound = false(length(mlb.obsvTimeVect));
        pIpOrwDbound(preTrlLog,preTrlLog)=true; 
        pIpOrwDbound(odrSampLog,odrSampLog)=true;
        pIpOrwDbound(preRwdPrd,preRwdPrd)=true;
        pIpOrwDbound(pstTrlPrd,pstTrlPrd)=true;
%     persistMod = cell2mat(permute([{nullModel}, {pObound}, {pIbound}, {pIpObound}, {pIpOrwDbound}], [1,3,2])); persXtickLabels = [{'Null'}, {'PO'}, {'PI'}, {'PIPO'}, {'PIPORWD'}];
    persistMod = cell2mat(permute([{pObound}, {pIbound}, {pIpObound}, {pIpOrwDbound}], [1,3,2])); persXtickLabels = [{'PO'}, {'PI'}, {'PIPO'}, {'PIPORWD'}];
    
    %% Now create models with persistent epochs
    prePrst = false(length(mlb.obsvTimeVect));
        prePrst(preTrlLog,preTrlLog) = true;
    odrPrst = false(length(mlb.obsvTimeVect));
        odrPrst(odrSampLog,odrSampLog) = true;
    pstPrst = false(length(mlb.obsvTimeVect));
        pstPrst(preRwdPrd|pstTrlPrd,preRwdPrd|pstTrlPrd) = true;
    preOdrPrst = false(length(mlb.obsvTimeVect));
        preOdrPrst(preTrlLog|odrSampLog, preTrlLog|odrSampLog) = true;
    prePstPrst = prePrst | pstPrst;
    
    % Now create a fully mixed model hinging at port entry where persistence grows after poke in
    odrGrow = false(length(mlb.obsvTimeVect));
    odrGrow(preTrlLog,preTrlLog) = true;
    frstNdx = find(odrSampLog,1,'first');
    startNdx = frstNdx;
    for ndx = frstNdx:length(mlb.obsvTimeVect)
        tempLog = false(size(mlb.obsvTimeVect));
        maxNdx = ndx+(ndx-frstNdx);
        if maxNdx > length(mlb.obsvTimeVect)
            maxNdx = length(mlb.obsvTimeVect);
        end
        tempLog(ndx:maxNdx) = true;
        odrGrow(ndx,tempLog) = true;
        odrGrow(tempLog,ndx) = true;
    end
    
    mixModFrame = cell2mat(permute([{prePrst}, {odrPrst}, {pstPrst}, {prePstPrst}, {preOdrPrst}, {odrGrow}], [1,3,2])); mixXtickLabels = [{'Pre'}, {'Odr'}, {'Pst'}, {'PrePst'}, {'PreOdr'}, {'OdrGrow'}];
    
    %% Finally identify whether these models, or mixtures of the two, are best fits to the data.
    
    permDynSimsKL = nan(numPerms,4,size(fullDynMod,3));
    permDynSimsEUC = nan(numPerms,4,size(fullDynMod,3));
    permDynSimsCOS = nan(numPerms,4,size(fullDynMod,3));
    permPersSimsKL = nan(numPerms,4,size(persistMod,3));
    permPersSimsEUC = nan(numPerms,4,size(persistMod,3));
    permPersSimsCOS = nan(numPerms,4,size(persistMod,3));
    permMixModKL = nan(numPerms,4,size(fullDynMod,3),size(mixModFrame,3));
    permMixModEUC = nan(numPerms,4,size(fullDynMod,3),size(mixModFrame,3));
    permMixModCOS = nan(numPerms,4,size(fullDynMod,3),size(mixModFrame,3));
    for perm = 1:numPerms
        tempPerm = dPrmPerm{perm};
        for o = 1:mlb.seqLength
            tempPermZ = zscore(tempPerm(:,:,o),0,'all')+100;
            for d = 1:size(fullDynMod,3)
                tempDynModelZ = zscore(fullDynMod(:,:,d),0,'all')+100;
                tempKLmtx = tempPermZ.*log(tempPermZ./tempDynModelZ);
                permDynSimsKL(perm,o,d) = sum(tempKLmtx(:));
                permDynSimsEUC(perm,o,d) = pdist([tempPermZ(:)'; tempDynModelZ(:)']);
                permDynSimsCOS(perm,o,d) = pdist([tempPermZ(:)'; tempDynModelZ(:)']-100, 'cosine');
                for m = 1:size(mixModFrame,3)
                    tempMixModelZ = zscore(fullDynMod(:,:,d)|mixModFrame(:,:,m),0,'all')+100;
                    tempKLmtx = tempPermZ.*log(tempPermZ./tempMixModelZ);
                    permMixModKL(perm,o,d,m) = sum(tempKLmtx(:));
                    permMixModEUC(perm,o,d,m) = pdist([tempPermZ(:)'; tempMixModelZ(:)']);
                    permMixModCOS(perm,o,d,m) = pdist([tempPermZ(:)'; tempMixModelZ(:)']-100, 'cosine');
                end
            end
            for p = 1:size(persistMod,3)
                tempPersModelZ = zscore(persistMod(:,:,p),0,'all')+100;
                tempKLmtx = tempPermZ.*log(tempPermZ./tempPersModelZ);
                permPersSimsKL(perm,o,p) = sum(tempKLmtx(:));
                permPersSimsEUC(perm,o,p) = pdist([tempPermZ(:)'; tempPersModelZ(:)']);
                permPersSimsCOS(perm,o,p) = pdist([tempPermZ(:)'; tempPersModelZ(:)']-100, 'cosine');
            end
        end
    end    
    dynModSims(:,ani) = [{permDynSimsKL};{permDynSimsEUC};{permDynSimsCOS}];
    persModSims(:,ani) = [{permPersSimsKL};{permPersSimsEUC};{permPersSimsCOS}];
    mixModSims(:,ani) = [{permMixModKL};{permMixModEUC};{permMixModCOS}];
        
    %% Plot per animal
    tempDynKL = permute(mean(permDynSimsKL),[3,2,1]);
    tempPersKL = permute(mean(permPersSimsKL),[3,2,1]);
    tempMixKL = permute(mean(permMixModKL),[3,4,2,1]);
    
    figure;
    subplot(3,1,1)
    hold(gca, 'on');
    for o = 1:mlb.seqLength
        plot(tempDynKL(:,o), 'color', mlb.PositionColors(o,:));
        ndx = find(min(tempDynKL(:,o))==tempDynKL(:,o),1,'first');
        scatter(ndx,min(tempDynKL(:,o)),'x', 'markeredgecolor', mlb.PositionColors(o,:));
    end
    title(min(min(tempDynKL)));
    
    subplot(3,1,2)
    hold(gca,'on');
    for o = 1:mlb.seqLength
        plot(tempPersKL(:,o), 'color', mlb.PositionColors(o,:));
        ndx = find(min(tempPersKL(:,o))==tempPersKL(:,o),1,'first');
        scatter(ndx,min(tempPersKL(:,o)),'x', 'markeredgecolor', mlb.PositionColors(o,:));
    end
    set(gca,'xtick', 1:size(tempPersKL,1), 'xticklabels', persXtickLabels)
    title(min(min(tempPersKL)));
    
    for o = 1:mlb.seqLength
        subplot(3,mlb.seqLength, sub2ind([mlb.seqLength, 3],o,3));
        imagesc(tempMixKL(:,:,o));
        hold on;
        [d,p] = find(min(tempMixKL(:,:,o))==tempMixKL(:,:,o));
        scatter(p,d, 'wx');
        [d,p] = find(min(min(tempMixKL(:,:,o)))==tempMixKL(:,:,o));
        scatter(p,d, 'wo');
        set(gca, 'ydir', 'normal', 'xtick', 1:size(tempMixKL,2), 'xticklabels', mixXtickLabels, 'xticklabelrotation', 90);
        title(min(min(tempMixKL(:,:,o))));
    end
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    drawnow;
    
end