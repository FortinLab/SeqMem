tets = [];
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% setupSeqLength = 4; 

% CA1 Data
fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
tets = [1,22,17,18,17]; % Lateral/Proximal
% tets = [7,3,1,5,5]; % Medial/Distal
setupSeqLength = 5;

% % All well-trained files PFC
% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4; 

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
% Neural Variables
aniPosts = cell(size(fileDirs));
aniOdorDecodes = cell(size(fileDirs));
aniTimeDecodes = cell(size(fileDirs));
% aniPosAccLog = cell(size(fileDirs));
% aniTimeAccLog = cell(size(fileDirs));
aniTrlPosIDs = cell(size(fileDirs));
aniPermIDs = cell(size(fileDirs));
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
    [betaPhase, betaPower] = mlb.PP_TrialMatrix_LFP([20 40], trlWindow{1}, alignment{1});
    [thetaPhase, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
    %% Decode ISC via Sub-Sampling
%     mlb.bayesType = 1;  % Comment In to decode using Gaussian rather than Poisson
    mlb.SetLikes_SubSample;
%     mlb.bayesType = 3;  % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_Observes;
    tempPosts = cell2mat(reshape(mlb.post, [1,1,numel(mlb.post)]));
    postTrlIDs = cell2mat(reshape(mlb.postTrlIDs, [1,1,numel(mlb.postTrlIDs)]));
    tempOdorDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,3));
%     tempPosAccuracy = arrayfun(@(a,b)a==b,repmat(tempOdorDecode,[1,1,mlb.seqLength]), repmat(permute(1:mlb.seqLength,[3,1,2]), size(tempOdorDecode,1),size(tempOdorDecode,2),1));
    decodeTimeVect = round(mlb.decodeIDvects{1}(:,1)*1000)/1000;
    tempTimeDecode = mlb.DecodeBayesPost(tempPosts, decodeTimeVect);
    timePoints = unique(decodeTimeVect);
%     tempTimeAccuracy = arrayfun(@(a,b)a==b,repmat(tempTimeDecode,[1,1,length(timePoints)]), repmat(permute(timePoints,[3,2,1]), size(tempTimeDecode,1),size(tempTimeDecode,2),1));
    tempLFP = permute([betaPower(:,:,postTrlIDs), thetaPower(:,:,postTrlIDs)], [1,3,2]);
    permIDvect = cellfun(@(a,b){repmat(a(1,5), [1,size(b,3)])}, mlb.decodeIDvects, mlb.post);
    %% Outputs
    aniPosts{ani} = cell2mat(reshape(mlb.post, [1,1,mlb.numPerms]));
    aniOdorDecodes{ani} = tempOdorDecode;
    aniTimeDecodes{ani} = tempTimeDecode;
%     aniPosAccLog{ani} = tempPosAccuracy;
%     aniTimeAccLog{ani} = tempTimeAccuracy;
    aniTrlPosIDs{ani} = [mlb.trialInfo(postTrlIDs).Position];
    aniPermIDs{ani} = cell2mat(permIDvect);
    aniLFP{ani} = tempLFP;
end
pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

piNdx = find(abs(mlb.likeTimeVect)==min(abs(mlb.likeTimeVect)))+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
%% Create group trial data & ID vectors
if strcmp(alignment{1}, 'PokeIn')
    trlTimeLog = mlb.obsvTimeVect>0;
else
    trlTimeLog = mlb.obsvTimeVect>-0.5;
end
grpTrlPosts = cell2mat(reshape(aniPosts, [1,1,length(fileDirs)]));
grpTrlLFP = cell2mat(aniLFP);
grpTrlTimeDecodes = cell2mat(aniTimeDecodes);
grpTrlPosDecodes = cell2mat(aniOdorDecodes);

grpTrlPos = cell2mat(aniTrlPosIDs);
grpTrlPermIDs = cell2mat(aniPermIDs);
%% Determine LFP Power thresholds
for band = 1:size(grpTrlLFP,3)    
    tempTrlMeanPower = mean(grpTrlLFP(trlTimeLog,:,band));
    lfpThresh = nan(mlb.seqLength,2);
    for pos = 1:mlb.seqLength
        sortedPower = sort(tempTrlMeanPower(grpTrlPos==pos));
%         sortedPower = sort(tempTrlMeanPower);
%         betaThresh = [mean(betaTrlMean(grpTrlPos==p))-(std(betaTrlMean(grpTrlPos==p))*1), mean(betaTrlMean(grpTrlPos==p))+(std(betaTrlMean(grpTrlPos==p))*1)]; % INDIVIDUAL threshold (Mean+/-STD)
%         betaThresh = [mean(betaTrlMean(grpTrlPos))-(std(betaTrlMean(grpTrlPos))*1), mean(betaTrlMean(grpTrlPos))+(std(betaTrlMean(grpTrlPos))*1)]; % SESSION threshold (Mean +/-STD)
        lfpThresh(pos,:) = [sortedPower(ceil(length(sortedPower)*0.25)), sortedPower(floor(length(sortedPower)*0.75))];
    end
    % Split data based on LFP power
    splitPosts = cell(mlb.seqLength, 2, numPerms);
    splitTimeDecode = cell(mlb.seqLength, 2, numPerms);
    splitPosDecode = cell(mlb.seqLength, 2, numPerms);
    splitOPdPrmDecode = cell(1,2,numPerms);
    for perm = 1:numPerms
        curPermLog = grpTrlPermIDs==perm;
        for pos = 1:mlb.seqLength
            curPosLog = grpTrlPos==pos;
            curLowPwrTrlLog = tempTrlMeanPower<lfpThresh(pos,1);
            curHighPwrTrlLog = tempTrlMeanPower>lfpThresh(pos,2);
            
            splitPosts{pos,1,perm} = mean(grpTrlPosts(:,:,curPosLog & curPermLog & curLowPwrTrlLog),3,'omitnan');
            splitPosts{pos,2,perm} = mean(grpTrlPosts(:,:,curPosLog & curPermLog & curHighPwrTrlLog),3,'omitnan');
            
            tempLPTtimeDecode = nan(length(timePoints));
            tempHPTtimeDecode = nan(length(timePoints));
            for t1 = 1:length(timePoints)
                for t2 = 1:length(timePoints)
                    tempLPTtimeDecode(t1,t2) = mean(grpTrlTimeDecodes(t1,curPosLog & curPermLog & curLowPwrTrlLog)==timePoints(t2), 'omitnan');
                    tempHPTtimeDecode(t1,t2) = mean(grpTrlTimeDecodes(t1,curPosLog & curPermLog & curHighPwrTrlLog)==timePoints(t2), 'omitnan');
                end
            end                
            splitTimeDecode{pos,1,perm} = tempLPTtimeDecode;
            splitTimeDecode{pos,2,perm} = tempHPTtimeDecode;
            
            tempLPTodrDecode = nan(length(timePoints), mlb.seqLength);
            tempHPTodrDecode = nan(length(timePoints), mlb.seqLength);
            for odr = 1:mlb.seqLength
                tempLPTodrDecode(:,odr) = mean(grpTrlPosDecodes(:,curPosLog & curPermLog & curLowPwrTrlLog)==odr,2,'omitnan');
                tempHPTodrDecode(:,odr) = mean(grpTrlPosDecodes(:,curPosLog & curPermLog & curHighPwrTrlLog)==odr,2,'omitnan');
            end
            splitPosDecode{pos,1,perm} = tempLPTodrDecode;
            splitPosDecode{pos,2,perm} = tempHPTodrDecode;
        end
        splitOPdPrmDecode{1,1,perm} = mlb.CalcDprmVectFromDecode(grpTrlPosDecodes(:,curPermLog & curLowPwrTrlLog), grpTrlPos(curPermLog & curLowPwrTrlLog));
        splitOPdPrmDecode{1,2,perm} = mlb.CalcDprmVectFromDecode(grpTrlPosDecodes(:,curPermLog & curHighPwrTrlLog), grpTrlPos(curPermLog & curHighPwrTrlLog));
    end
    
    %% Now Plot things
%     h.XRuler.FirstCrossoverValue = 0;
%     h.XRuler.SecondCrossoverValue = 0;
    if band ==1
        bnm = 'Beta';
    elseif band == 2
        bnm = 'Theta';
    end
    curHighPowTrlPost = mean(cell2mat(splitPosts(:,2,:)), 3, 'omitnan'); cLim = [0 0.05];
%     curHighPowTrlPost = curHighPowTrlPost./max(curHighPowTrlPost(:));  cLim = [0 0.25];
%     curHighPowTrlPost = (curHighPowTrlPost-(mean(curHighPowTrlPost(:))))./std(curHighPowTrlPost(:));  cLim = [-6 6];    
    curLowPowTrlPost = mean(cell2mat(splitPosts(:,1,:)), 3, 'omitnan');
%     curLowPowTrlPost = curLowPowTrlPost./max(curLowPowTrlPost(:));
%     curLowPowTrlPost = (curLowPowTrlPost-(mean(curLowPowTrlPost(:))))./std(curLowPowTrlPost(:));
    imsp = nan(1,3);
    pasp = nan(1,3);
    tasp = nan(1,3);
    tdsp = nan(1,3);
    
    taClim = [0 0.15];
    
    figure;
    % High Power Trials
    % Time Accuracy HPT
    tasp(1) = subplot(3,9, [1,10]);
    highPowTimeAcc = mean(cell2mat(splitTimeDecode(:,2,:)), 3, 'omitnan');
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
    % Time Decoding HPT
    hptTimeAccDiag = mean(cell2mat(cellfun(@(a){diag(a)'}, splitTimeDecode(:,2,:))), 'omitnan');
    meanHPTtimeDiag = mean(hptTimeAccDiag,3);
%     varHPTtimeDiag = mlb.SEMcalc(hptTimeAccDiag,0,3);
    varHPTtimeDiag = std(hptTimeAccDiag,0,3);
    tdsp(1) = subplot(3,9,19);
    plot(meanHPTtimeDiag, 'k');
    hold on;    
    patch('XData', [1:length(timePoints), length(timePoints):-1:1],...
        'YData', [(meanHPTtimeDiag+varHPTtimeDiag), fliplr(meanHPTtimeDiag-varHPTtimeDiag)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '-');
    set(gca, 'ylim', [0 0.25], 'xticklabel', []);
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);    
    title('Mean Diag');
    % HPT Posteriors
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
    % Position Accuracy HPT
    hptPosAcc = cell2mat(splitPosDecode(:,2,:));
    meanHPTposAcc = mean(hptPosAcc,3);
%     varHPTposAcc = mlb.SEMcalc(hptPosAcc,0,3);
    varHPTposAcc = std(hptPosAcc,0,3);
    pasp(1) = subplot(3,9,20:21);
    hold on;
    for op = 1:mlb.seqLength
        plot(meanHPTposAcc(:,op), 'color', mlb.PositionColors(op,:));
        patch('XData', [1:size(meanHPTposAcc,1), size(meanHPTposAcc,1):-1:1],...
            'YData', [(meanHPTposAcc(:,op)+varHPTposAcc(:,op)); flipud(meanHPTposAcc(:,op)-varHPTposAcc(:,op))], 'edgecolor', mlb.PositionColors(op,:),...
            'facecolor', mlb.PositionColors(op,:), 'facealpha', 0.5, 'linestyle', '-');
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
    
        
    
    % Low Power Trials
    % Time Accuracy LPT
    tasp(2) = subplot(3,9, [4,13]);
    lowPowTimeAcc = mean(cell2mat(splitTimeDecode(:,1,:)), 3, 'omitnan');
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
    % Time Decoding HPT
    lptTimeAccDiag = mean(cell2mat(cellfun(@(a){diag(a)'}, splitTimeDecode(:,1,:))), 'omitnan');
    meanLPTtimeDiag = mean(lptTimeAccDiag,3);
%     varLPTtimeDiag = mlb.SEMcalc(lptTimeAccDiag,0,3);
    varLPTtimeDiag = std(lptTimeAccDiag,0,3);
    tdsp(2) = subplot(3,9,22);
    plot(meanLPTtimeDiag, 'k');
    hold on;    
    patch('XData', [1:length(timePoints), length(timePoints):-1:1],...
        'YData', [(meanLPTtimeDiag+varLPTtimeDiag), fliplr(meanLPTtimeDiag-varLPTtimeDiag)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '-');
    set(gca, 'ylim', [0 0.25], 'xticklabel', []);
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    title('Mean Diag');
    % LPT Posteriors
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
    % Position Accuracy LPT
    lptPosAcc = cell2mat(splitPosDecode(:,1,:));
    meanLPTposAcc = mean(lptPosAcc,3);
%     varLPTposAcc = mlb.SEMcalc(lptPosAcc,0,3);
    varLPTposAcc = std(lptPosAcc,0,3);
    pasp(2) = subplot(3,9,23:24);
    hold on;
    for op = 1:mlb.seqLength
        plot(meanLPTposAcc(:,op), 'color', mlb.PositionColors(op,:));
        patch('XData', [1:size(meanLPTposAcc,1), size(meanLPTposAcc,1):-1:1],...
            'YData', [(meanLPTposAcc(:,op)+varLPTposAcc(:,op)); flipud(meanLPTposAcc(:,op)-varLPTposAcc(:,op))], 'edgecolor', mlb.PositionColors(op,:),...
            'facecolor', mlb.PositionColors(op,:), 'facealpha', 0.5, 'linestyle', '-');
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
    
    
    % Power Difference
    % Time Accuracy HPT - LPT
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
    % Time Decoding HPT-LPT
    diffTimeAccDiag = hptTimeAccDiag-lptTimeAccDiag;
    meanDIFFtimeDiag = mean(diffTimeAccDiag,3);
%     varDIFFtimeDiag = mlb.SEMcalc(diffTimeAccDiag,0,3);
    varDIFFtimeDiag = std(diffTimeAccDiag,0,3);
    tdsp(3) = subplot(3,9,25);
    plot(meanDIFFtimeDiag, 'k');
    hold on;    
    patch('XData', [1:length(timePoints), length(timePoints):-1:1],...
        'YData', [(meanDIFFtimeDiag+varDIFFtimeDiag), fliplr(meanDIFFtimeDiag-varDIFFtimeDiag)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '-');
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(get(gca, 'xlim'), [0 0], '-k', 'linewidth', 1);
    set(gca, 'xticklabel', []);
    title('TimeDiag Diff');
    % HPT-LPT Posteriors
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
    % Position Accuracy HPT-LPT
    diffPosAcc = hptPosAcc-lptPosAcc;
    meanDIFFposAcc = mean(diffPosAcc,3);
%     varDIFFposAcc = mlb.SEMcalc(diffPosAcc,0,3);
    varDIFFposAcc = std(diffPosAcc,0,3);
    pasp(3) = subplot(3,9,26:27);
    hold on;
    for op = 1:mlb.seqLength
        plot(meanDIFFposAcc(:,op), 'color', mlb.PositionColors(op,:));
        patch('XData', [1:size(meanDIFFposAcc,1), size(meanDIFFposAcc,1):-1:1],...
            'YData', [(meanDIFFposAcc(:,op)+varDIFFposAcc(:,op)); flipud(meanDIFFposAcc(:,op)-varDIFFposAcc(:,op))], 'edgecolor', mlb.PositionColors(op,:),...
            'facecolor', mlb.PositionColors(op,:), 'facealpha', 0.25, 'linestyle', '-');
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
    
    %% Plot Summaries of Decoding Difference
    boundNdx = [find(mlb.likeTimeVect==min(mlb.likeTimeVect)); length(mlb.likeTimeVect)+1];
    avgAccDiffPerPerm = nan(length(timePoints),numPerms);
    avgDprmDiffPerPerm = nan(length(timePoints),numPerms);
    for perm = 1:numPerms
        tempMatch = nan(length(timePoints), mlb.seqLength);
        for p = 1:mlb.seqLength
            tempMatch(:,p) = diffPosAcc(boundNdx(p):boundNdx(p+1)-1,p,perm);
        end
        avgAccDiffPerPerm(:,perm) = mean(tempMatch, 2, 'omitnan');
        avgDprmDiffPerPerm(:,perm) = mean(splitOPdPrmDecode{1,2,perm}-splitOPdPrmDecode{1,1,perm}, 2, 'omitnan');
    end
    
    figure; 
    subplot(2,1,1);
    meanAccDiffPerPerm = mean(avgAccDiffPerPerm,2);
%     varAccDiffPerPerm = mlb.SEMcalc(avgAccDiffPerPerm,0,2);
    varAccDiffPerPerm = std(avgAccDiffPerPerm,0,2);
    plot(meanAccDiffPerPerm, 'k');
    hold on;    
    patch('XData', [1:length(timePoints), length(timePoints):-1:1],...
        'YData', [(meanAccDiffPerPerm+varAccDiffPerPerm); flipud(meanAccDiffPerPerm-varAccDiffPerPerm)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '-');
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(get(gca, 'xlim'), [0 0], '-k', 'linewidth', 1);
    set(gca, 'xticklabel', []);
    ylabel('Accuracy Diff');
    title('Average HPT-LPT Difference');
    
    subplot(2,1,2);
    meanDprmDiffPerPerm = mean(avgDprmDiffPerPerm,2);
%     varDprmDiffPerPerm = mlb.SEMcalc(avgDprmDiffPerPerm,0,2);
    varDprmDiffPerPerm = std(avgDprmDiffPerPerm,0,2);
    plot(meanDprmDiffPerPerm, 'k');
    hold on;    
    patch('XData', [1:length(timePoints), length(timePoints):-1:1],...
        'YData', [(meanDprmDiffPerPerm+varDprmDiffPerPerm); flipud(meanDprmDiffPerPerm-varDprmDiffPerPerm)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '-');
    plot(repmat(piNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(poNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(repmat(rwdNdx(1), [1,2]), get(gca, 'ylim'), '--k','linewidth', 1);
    plot(get(gca, 'xlim'), [0 0], '-k', 'linewidth', 1);
    set(gca, 'xticklabel', []);
    ylabel('d'' Diff');
    title('Average HPT-LPT Difference');
        
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s Split; %s aligned; Trial Window = (%.0fms:%.0fms); %i Perms', bnm, alignment{1}, trlWindow{1}(1), trlWindow{1}(2),numPerms),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
    