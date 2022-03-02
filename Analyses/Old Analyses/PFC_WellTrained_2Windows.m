clear all; %#ok<CLALL>
%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

binSize = 200;
dsRate = 50;
trlWindow = {[-800 500]; [-500 800]};
alignment = [{'PokeIn'}, {'PokeOut'}];
lfpWindow = [16 32];
numPerms = 10;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
% OSC Decoding
% Posteriors
tMatPosts = cell(1,1,length(fileDirs));
postSimPosOO = cell(1,1,length(fileDirs));
postSimOdrOO = cell(1,1,length(fileDirs));
postSimPosOI = cell(1,1,length(fileDirs));
postSimOdrOI = cell(1,1,length(fileDirs));
postSimLag = cell(1,1,length(fileDirs));
% Decodings
tMatTimeDecode = cell(1,1,length(fileDirs));
tMatPosDecode = cell(1,1,length(fileDirs));
tMatOdrDecode = cell(1,1,length(fileDirs));
oscPosAccuracy = cell(1,1,length(fileDirs));
oscOdrAccuracy = cell(1,1,length(fileDirs));

% TAO Decoding
taoMtxPosts = cell(1,1,length(fileDirs));
taoMtxTimeDecode = cell(1,1,length(fileDirs));
taoMtxItemDecode = cell(1,1,length(fileDirs));
taoPrevPosAccuracy = cell(1,1,length(fileDirs));
taoNextPosAccuracy = cell(1,1,length(fileDirs));
taoPrevOdrAccuracy = cell(1,1,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
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
    
    if ani == 1
        tMatPosts = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        tMatTimeDecode = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        tMatOdrDecode = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        tMatPosDecode = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        taoMtxPosts = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        taoMtxTimeDecode = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        taoMtxItemDecode = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
    end
    
    %% Decode OSC via FISC 
    mlb.SetLikes_FISC;
    mlb.SetObserves_Session;
    mlb.Process_Observes;
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Time')));
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Odor')));
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Position')));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    %% Organize decodings into TransMat Plots
    for o = 1:mlb.seqLength
        for p = 1:mlb.seqLength
            opLog = [mlb.trialInfo.Odor]==o & [mlb.trialInfo.Position]==p & [mlb.trialInfo.Performance]==1;
            if sum(opLog)>=1
                if o==p     %only works with well-trained sessions... need to adjust if running with Novel or Dual List sessions
                    %% Comment in to run with ONLY FISC trials
%                     opLog = false(size(opLog));
%                     opLog(mlb.fiscTrials(o,:)) = true;
                    %% Comment in to run with NON FISC OR TAO trials
                    % First remove the TAO trials
                    opTrls = find(opLog);
                    for opt = 1:length(opTrls)
                        if opTrls(opt)==1 || mlb.trialInfo(opTrls(opt)-1).TranspositionDistance==0
                            continue;
                        else
                            opLog(opTrls(opt)) = false;
                        end
                    end
                    % Then remove the FISC trials
                    opLog(mlb.fiscTrials(o,:)) = false;
                end
                tMatPosts{o,p,ani} = mean(mlb.post{1}(:,:,opLog),3);
                % Temporal Accuracy
                tempTimeDecode = timeDecode(:,opLog);
                tempTimeAccuracy = nan(size(timeDecode,1), length(timePoints));
                for t = 1:size(timePoints,1)
                    tempTimeAccuracy(:,t) = mean(tempTimeDecode==timePoints(t),2);
                end
                tMatTimeDecode{o,p,ani} = tempTimeAccuracy;
                % Ordinal Accuracy
                tempPosDecode = posDecode(:,opLog);
                tempPosAccuracy = nan(size(posDecode,1), mlb.seqLength);
                for pp = 1:mlb.seqLength
                    tempPosAccuracy(:,pp) = mean(tempPosDecode==pp,2);
                end
                tMatPosDecode{o,p,ani} = tempPosAccuracy;
                % Odor Accuracy
                tempOdrDecode = odrDecode(:,opLog);
                tempOdrAccuracy = nan(size(odrDecode,1), mlb.seqLength);
                for oo = 1:mlb.seqLength
                    tempOdrAccuracy(:,oo) = mean(tempOdrDecode==mlb.odrSeqs(oo),2);
                end
                tMatOdrDecode{o,p,ani} = tempOdrAccuracy;
            end
        end
    end       
    %% Evaluate Posteriors for similarity based on position, odor and lag
    tempPosts = tMatPosts(:,:,ani);
    tempPosts = cellfun(@(a)a(:),tempPosts, 'uniformoutput', 0);
    tempPosts(cellfun(@(a)isempty(a),tempPosts)) = {nan(size(tempPosts{1}))};
    lagSim = nan(1,length(mlb.lagVect));
    for lag = 1:length(mlb.lagVect)
        tempTMatLog = triu(true(mlb.seqLength),mlb.lagVect(lag)) & tril(true(mlb.seqLength),mlb.lagVect(lag));
        tempCorrMat = corr(cell2mat(tempPosts(tempTMatLog)'));
        lagSim(lag) = mean(tempCorrMat(triu(true(size(tempCorrMat)),1)), 'omitnan');
    end
    postSimLag{1,1,ani} = lagSim;
    posSim = nan(1,mlb.seqLength);
    isPosSim = nan(1,mlb.seqLength);
    odrSim = nan(1,mlb.seqLength);  
    isOdrSim = nan(1,mlb.seqLength);
    for op = 1:mlb.seqLength
        tempPosCorr = corr(cell2mat(tempPosts((1:mlb.seqLength)~=op,op)'));
        posSim(op) = mean(tempPosCorr(triu(true(size(tempPosCorr)),1)), 'omitnan');
        tempOdrCorr = corr(cell2mat(tempPosts(op,(1:mlb.seqLength)~=op)));
        odrSim(op) = mean(tempOdrCorr(triu(true(size(tempOdrCorr)),1)), 'omitnan');
        
        tempISposCorr = nan(1,mlb.seqLength);
        tempISodrCorr = nan(1,mlb.seqLength);
        for op2 = 1:mlb.seqLength
            if op2==op
                continue
            else
                 tempPosCorr = corr([tempPosts{op2,op}, tempPosts{op,op}]);
                 tempISposCorr(op2) = tempPosCorr(2);
                 tempOdrCorr = corr([tempPosts{op,op2}, tempPosts{op,op}]);
                 tempISodrCorr(op2) = tempOdrCorr(2);
            end
        end
        isPosSim(op) = mean(tempISposCorr, 'omitnan');
        isOdrSim(op) = mean(tempISodrCorr, 'omitnan');
    end    
    postSimPosOO{1,1,ani} = posSim;
    postSimOdrOO{1,1,ani} = odrSim;
    postSimPosOI{1,1,ani} = isPosSim;
    postSimOdrOI{1,1,ani} = isOdrSim;
    
    %% Evaluate OutSeq & Trial After OutSeq trials individually... a bit redundant to above but w/e all this is speedy meet me outside if you have a problem with it!
    osTrlNums = [mlb.trialInfo([mlb.trialInfo.TranspositionDistance]~=0).TrialNum];
    osTrlNums([mlb.trialInfo(osTrlNums).Performance]==0) = [];
    tempOSposAcc = nan(size(posDecode,1),length(osTrlNums));
    tempOSodrAcc = nan(size(odrDecode,1),length(osTrlNums));
    tempTAOprevPosAcc = nan(size(posDecode,1),length(osTrlNums));
    tempTAOnextPosAcc = nan(size(posDecode,1),length(osTrlNums));
    tempTAOprevOdrAcc = nan(size(odrDecode,1),length(osTrlNums));
    tempTAOposts = cell(mlb.seqLength, mlb.seqLength,length(osTrlNums));
    tempTAOmtxTimeAcc = cell(mlb.seqLength, mlb.seqLength,length(osTrlNums));
    tempTAOmtxItemAcc = cell(mlb.seqLength, mlb.seqLength,length(osTrlNums));
    for ost = 1:length(osTrlNums)
        tempOSposAcc(:,ost) = posDecode(:,osTrlNums(ost))==mlb.trialInfo(osTrlNums(ost)).Position;
        tempOSodrAcc(:,ost) = odrDecode(:,osTrlNums(ost))==mlb.trialInfo(osTrlNums(ost)).Odor;
        if mlb.trialInfo(osTrlNums(ost)).Position ~= mlb.seqLength
            if mlb.trialInfo(osTrlNums(ost)+1).Performance==1 && osTrlNums(ost)~=length(mlb.trialInfo)
                % Identify TAO Odor/Position... just so things are easier to read in the code below
                taoOdr = mlb.trialInfo(osTrlNums(ost)+1).Odor;
                taoPos = mlb.trialInfo(osTrlNums(ost)+1).Position;
                % Extract posteriors
                tempTAOposts{taoOdr,taoPos,ost} = mlb.post{1}(:,:,osTrlNums(ost)+1);
                % Calculate temporal decoding accuracy
                tempTimeAccuracy = nan(size(timeDecode,1), length(timePoints));
                for t = 1:size(timeDecode,1)
                    tempTimeAccuracy(:,t) = timeDecode(t,osTrlNums(ost)+1)==timePoints;
                end
                tempTAOmtxTimeAcc{taoOdr,taoPos,ost} = tempTimeAccuracy;
                % Calculate position/previous odor decoding accuracy
                if mlb.trialInfo(osTrlNums(ost)).Odor ~= mlb.trialInfo(osTrlNums(ost)+1).Position
                    tempTAOprevPosAcc(:,ost) = posDecode(:,osTrlNums(ost)+1)==mlb.trialInfo(osTrlNums(ost)).Position;
                    tempTAOnextPosAcc(:,ost) = posDecode(:,osTrlNums(ost)+1)==mlb.trialInfo(osTrlNums(ost)+1).Position;
                    tempTAOprevOdrAcc(:,ost) = odrDecode(:,osTrlNums(ost)+1)==mlb.trialInfo(osTrlNums(ost)).Odor;
                end
                tempItemDecode = nan(size(posDecode,1), mlb.seqLength);
                for o = 1:mlb.seqLength
                    tempItemDecode(:,o) = posDecode(:,osTrlNums(ost)+1)==o;
                end
                tempTAOmtxItemAcc{taoOdr,taoPos,ost} = tempItemDecode;
                    
            end
        end
    end
    for o = 1:mlb.seqLength
        for p = 1:mlb.seqLength
            tempTAOpostsTM = tempTAOposts(o,p,:);
            rmvLog = cellfun(@(a)isempty(a), tempTAOpostsTM);
            tempTAOpostsTM(rmvLog) = [];
            taoMtxPosts{o,p,ani} = mean(cell2mat(tempTAOpostsTM),3,'omitnan');
            tempTAOtimeTM = tempTAOmtxTimeAcc(o,p,:);
            tempTAOtimeTM(rmvLog) = [];
            taoMtxTimeDecode{o,p,ani} = mean(cell2mat(tempTAOtimeTM),3,'omitnan');
            tempTAOitemTM = tempTAOmtxItemAcc(o,p,:);
            tempTAOitemTM(rmvLog) = [];
            taoMtxItemDecode{o,p,ani} = mean(cell2mat(tempTAOitemTM),3,'omitnan');
        end
    end
    oscPosAccuracy{ani} = mean(tempOSposAcc,2,'omitnan');
    oscOdrAccuracy{ani} = mean(tempOSodrAcc,2,'omitnan');
    taoPrevPosAccuracy{ani} = mean(tempTAOprevPosAcc,2,'omitnan');
    taoNextPosAccuracy{ani} = mean(tempTAOnextPosAcc,2,'omitnan');
    taoPrevOdrAccuracy{ani} = mean(tempTAOprevOdrAcc,2,'omitnan');
end

%% Plot Posteriors
mlb.likeTimeVect = round(mlb.likeTimeVect*100)/100;
mlb.obsvTimeVect = round(mlb.obsvTimeVect*100)/100;
piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000))+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;

figure;
for sp = 1:mlb.seqLength^2
    [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
    tempPost = mean(cell2mat(tMatPosts(spR,spC,:)),3,'omitnan');
    if isempty(tempPost)
        continue;
    else
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
        imagesc(1:length(mlb.likeTimeVect), mlb.obsvTimeVect, tempPost, postCLim);
        hold on;
        plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]), '--k', 'linewidth', 2);
        for ndx = 1:length(piNdx)
            plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'),  '--k', 'linewidth', 2);
            plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'),  '--k', 'linewidth', 2);
            if ndx<length(piNdx)
                plot(repmat(posNdx(ndx),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
            end
        end
        set(gca, 'ydir', 'normal', 'xtick', []);
        title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));        
    end
end

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("FISC Decode: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
%% Plot Temporal Decoding
figure;
for sp = 1:mlb.seqLength*mlb.seqLength
    [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
    tempTime = mean(cell2mat(tMatTimeDecode(spR,spC,:)),3,'omitnan');
    if isempty(tempTime)
        continue;
    else
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempTime, decodeCLim);
        hold on;
        plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]), '--k', 'linewidth', 2);
        plot([0 0], get(gca, 'xlim'), '--k', 'linewidth', 2);
        plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);        
        set(gca, 'ydir', 'normal', 'xtick', []);
        title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));        
    end
end

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Temporal Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot Odor & Position Decoding
trlTypeDecodes = [{tMatOdrDecode}, {tMatPosDecode}];
for cog = 1:length(trlTypeDecodes)
    curCog = trlTypeDecodes{cog};
    figure;
    for sp = 1:mlb.seqLength*mlb.seqLength
        [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
        if isempty(curCog{spR,spC,1})
            continue;
        else
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
            for plt = 1:mlb.seqLength
                tempPlt = cell2mat(cellfun(@(a)a(:,plt), curCog(spR,spC,:), 'uniformoutput',0));
                pltMean = mean(tempPlt,3,'omitnan');
                pltVar = mlb.SEMcalc(tempPlt,0,3);
                plot(mlb.obsvTimeVect, pltMean, 'color', mlb.PositionColors(plt,:));
                hold on;
                patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                    'YData', [(pltMean+pltVar); flipud(pltMean-pltVar)], 'edgecolor', mlb.PositionColors(plt,:),...
                    'facecolor', mlb.PositionColors(plt,:), 'facealpha', 0.5);
                axis tight
                set(gca, 'ylim', [0 1]);
                box off
                plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
                plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);
                title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));
            end
        end
    end
    if cog == 1
        annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
            'String', sprintf("Odor Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    else        
        annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
            'String', sprintf("Position Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    end
end

%% Plot TAO Posteriors
figure; 
for sp = 1:mlb.seqLength^2
    [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
    tempPost = mean(cell2mat(taoMtxPosts(spR,spC,:)),3,'omitnan');
    if isempty(tempPost)
        continue;
    else
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
        imagesc(1:length(mlb.likeTimeVect), mlb.obsvTimeVect, tempPost, postCLim);
        hold on;
        plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]), '--k', 'linewidth', 2);
        for ndx = 1:length(piNdx)
            plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'),  '--k', 'linewidth', 2);
            plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'),  '--k', 'linewidth', 2);
            if ndx<length(piNdx)
                plot(repmat(posNdx(ndx),[1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
            end
        end
        set(gca, 'ydir', 'normal', 'xtick', []);
        title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));        
    end
end

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("FISC Decode TAO: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot TAO Temporal Decoding
figure;
for sp = 1:mlb.seqLength*mlb.seqLength
    [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
    tempTime = mean(cell2mat(taoMtxTimeDecode(spR,spC,:)),3,'omitnan');
    if isempty(tempTime)
        continue;
    else
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempTime, decodeCLim);
        hold on;
        plot(get(gca, 'xlim'), [0 0], '--k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]), '--k', 'linewidth', 2);
        plot([0 0], get(gca, 'xlim'), '--k', 'linewidth', 2);
        plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);        
        set(gca, 'ydir', 'normal', 'xtick', []);
        title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));        
    end
end

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("TAO Temporal Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot TAO Item Decoding
figure;
for sp = 1:mlb.seqLength*mlb.seqLength
    [spR, spC] = ind2sub([mlb.seqLength, mlb.seqLength], sp);
    if isempty(taoMtxItemDecode{spR,spC,1})
        continue;
    else
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], spC, spR));
        for plt = 1:mlb.seqLength
            tempPlt = cell2mat(cellfun(@(a)a(:,plt), taoMtxItemDecode(spR,spC,:), 'uniformoutput',0));
            pltMean = mean(tempPlt,3,'omitnan');
            pltVar = mlb.SEMcalc(tempPlt,0,3);
            plot(mlb.obsvTimeVect, pltMean, 'color', mlb.PositionColors(plt,:));
            hold on;
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(pltMean+pltVar); flipud(pltMean-pltVar)], 'edgecolor', mlb.PositionColors(plt,:),...
                'facecolor', mlb.PositionColors(plt,:), 'facealpha', 0.5);
        end
        axis tight
        set(gca, 'ylim', [0 1]);
        box off
        plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
        plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);
        title(sprintf('%s/%i', mlb.Rosetta{spR}, spC));
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("TAO Item Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot OSC Summary Figure
oscPosMean = mean(cell2mat(oscPosAccuracy),3,'omitnan');
oscPosVar = mlb.SEMcalc(cell2mat(oscPosAccuracy),0,3);

oscOdrMean = mean(cell2mat(oscOdrAccuracy),3,'omitnan');
oscOdrVar = mlb.SEMcalc(cell2mat(oscOdrAccuracy),0,3);

figure;
subplot(1,4,1:3)
plot(mlb.obsvTimeVect, oscPosMean, '-k');
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(oscPosMean+oscPosVar); flipud(oscPosMean-oscPosVar)], 'edgecolor', 'k',...
    'facecolor', 'k', 'facealpha', 0.5);
plot(mlb.obsvTimeVect, oscOdrMean, '-r');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(oscOdrMean+oscOdrVar); flipud(oscOdrMean-oscOdrVar)], 'edgecolor', 'r',...
    'facecolor', 'r', 'facealpha', 0.5);
axis tight
set(gca, 'ylim', [0 1]);
box off
plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);
title('OSC Decoding Accuracy');
patch('XData', [0,0,max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000),max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000)],...
    'YData', [get(gca, 'ylim'), fliplr(get(gca, 'ylim'))], 'edgecolor', 'none',...
    'facecolor', 'k', 'facealpha', 0.15);

% Plot bar graphs of Trial Period decodings
trlLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<=max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000);
oscPosAcc = cell2mat(oscPosAccuracy);
barMeanOSCpos = mean(oscPosAcc(trlLog,:,:));

oscOdrAcc = cell2mat(oscOdrAccuracy);
barMeanOSCodr = mean(oscOdrAcc(trlLog,:,:));

subplot(1,4,4)
b = bar([mean(barMeanOSCpos), mean(barMeanOSCodr)]);
b.FaceColor = 'flat';
b.CData(1,:) = 0;
b.CData(2,:) = [1 0 0];
hold on;
errorbar([mean(barMeanOSCpos), mean(barMeanOSCodr)],...
    [mlb.SEMcalc(barMeanOSCpos,0,3), mlb.SEMcalc(barMeanOSCodr,0,3)],...
    'linestyle', 'none', 'color', 'k', 'CapSize', 0);
title('OSC Trial Accuracy');
set(gca, 'xticklabel', [{'Position'}, {'Odor'}], 'XTickLabelRotation', 45);

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("OSC Trial Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
%% Plot TAO Summary Figure
taoPrevPosMean = mean(cell2mat(taoPrevPosAccuracy),3,'omitnan');
taoPrevPosVar = mlb.SEMcalc(cell2mat(taoPrevPosAccuracy),0,3);

taoNextPosMean = mean(cell2mat(taoNextPosAccuracy),3,'omitnan');
taoNextPosVar = mlb.SEMcalc(cell2mat(taoNextPosAccuracy),0,3);

taoPrevOdrMean = mean(cell2mat(taoPrevOdrAccuracy),3,'omitnan');
taoPrevOdrVar = mlb.SEMcalc(cell2mat(taoPrevOdrAccuracy),0,3);

figure;
subplot(1,4,1:3)
plot(mlb.obsvTimeVect, taoPrevPosMean, '--k');
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(taoPrevPosMean+taoPrevPosVar); flipud(taoPrevPosMean-taoPrevPosVar)], 'edgecolor', 'k',...
    'facecolor', 'k', 'facealpha', 0.5, 'linestyle', '--');
plot(mlb.obsvTimeVect, taoPrevOdrMean, '--r');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(taoPrevOdrMean+taoPrevOdrVar); flipud(taoPrevOdrMean-taoPrevOdrVar)], 'edgecolor', [0.9 0.025 0.025],...
    'facecolor', [0.9 0.025 0.025], 'facealpha', 0.5, 'linestyle', '--');
plot(mlb.obsvTimeVect, taoNextPosMean, '-k');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(taoNextPosMean+taoNextPosVar); flipud(taoNextPosMean-taoNextPosVar)], 'edgecolor', 'k',...
    'facecolor', 'k', 'facealpha', 0.5);
axis tight
set(gca, 'ylim', [0 1]);
box off
plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
plot(repmat(max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000), [1,2]),get(gca, 'ylim'), '--k', 'linewidth', 2);
title('TAO Decoding Accuracy');
patch('XData', [repmat(mlb.obsvTimeVect(1), [1,2]), 0,0,],...
    'YData', [get(gca, 'ylim'), fliplr(get(gca, 'ylim'))], 'edgecolor', 'none',...
    'facecolor', 'k', 'facealpha', 0.15);

% Plot bar graphs of pre-trial decodings
preTrlLog = mlb.obsvTimeVect<=0;
taoPrevPosAcc = cell2mat(taoPrevPosAccuracy);
barMeanTAOpp = mean(taoPrevPosAcc(preTrlLog,:,:));

taoPrevOdrAcc = cell2mat(taoPrevOdrAccuracy);
barMeanTAOpo = mean(taoPrevOdrAcc(preTrlLog,:,:));

taoNextPosAcc = cell2mat(taoNextPosAccuracy);
barMeanTAOnp = mean(taoNextPosAcc(preTrlLog,:,:));

subplot(1,4,4)
b = bar([mean(barMeanTAOnp), mean(barMeanTAOpp), mean(barMeanTAOpo)]);
b.FaceColor = 'flat';
b.CData(1:2,:) = 0;
b.CData(3,:) = [1 0 0];
hold on;
errorbar([mean(barMeanTAOnp), mean(barMeanTAOpp), mean(barMeanTAOpo)], [mlb.SEMcalc(barMeanTAOnp,0,3),...
    mlb.SEMcalc(barMeanTAOpp,0,3), mlb.SEMcalc(barMeanTAOpo,0,3)],...
    'linestyle', 'none', 'color', 'k', 'CapSize', 0);
title('Pre-Trial Accuracy');
set(gca, 'xticklabel', [{'Position'}, {'Prev Position'}, {'Prev Odor'}], 'XTickLabelRotation', 45);

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("TAO Pre-Trial Decoding: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Evaluate Posterior Similarities 
figure;
% Plot Lag Similarities
postSimLags = cell2mat(postSimLag);
subplot(3,2,1)
errorbar(mlb.lagVect, mean(postSimLags,3), mlb.SEMcalc(postSimLags,0,3));
title('Lag Similarity');

% Plot OutSeq-OutSeq Similarities
subplot(3,2,3)
psPosOO = cell2mat(postSimPosOO);
psOdrOO = cell2mat(postSimOdrOO);
bar([mean(psPosOO,3)',mean(psOdrOO,3)']);
title('Similarity To OutSeqs');

% Plot OutSeq-InSeq Similarities
subplot(3,2,5)
psPosOI = cell2mat(postSimPosOI);
psOdrOI = cell2mat(postSimOdrOI);
bar([mean(psPosOI,3)', mean(psOdrOI,3)']);
title('Similarity To InSeqs');

% Compare Across Different Siliarity Contrasts
subplot(3,2,2:2:6)
aniLagSim = reshape(mean(postSimLags, 'omitnan'), [1,length(fileDirs)]);
aniOdrSim = reshape(mean(psOdrOO, 'omitnan'), [1,length(fileDirs)]);
aniPosSim = reshape(mean(psPosOO, 'omitnan'), [1,length(fileDirs)]);
bar([mean(aniLagSim), mean(aniOdrSim), mean(aniPosSim)]);
hold on;
errorbar([mean(aniLagSim), mean(aniOdrSim), mean(aniPosSim)], [mlb.SEMcalc(aniLagSim,0,2), mlb.SEMcalc(aniOdrSim,0,2), mlb.SEMcalc(aniPosSim,0,2)], 'linestyle', 'none', 'color', 'k', 'capsize', 0);
plot([aniLagSim;aniOdrSim;aniPosSim], 'k');
set(gca, 'xtick', 1:3, 'xticklabels', [{'Lag'},{'Odor'}, {'Position'}], 'XTickLabelRotation', 45)
title('OutSeq Posterior Similarity');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("OSC Posterior Similarity: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
