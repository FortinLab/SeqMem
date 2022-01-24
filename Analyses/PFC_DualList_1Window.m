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
ssProportion = 0.5;
ssType = 0; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Data Outputs
% Behavior Variables
fiscPokeOutLat = cell(length(fileDirs),1);
fiscRwdDelivLat = cell(length(fileDirs),1);
smi = nan(length(fileDirs),2);
dPrm = nan(length(fileDirs),2);
ri = nan(length(fileDirs),2);
smiByOP = nan(length(fileDirs),8,2);
dPrmByOP = nan(length(fileDirs),8,2);
riByOP = nan(length(fileDirs),8,2);
% Decoding Variables
grpPosts = cell(1,1,length(fileDirs));
grpTimeDecode = cell(1,1,length(fileDirs));
grpPosDecode = cell(1,1,length(fileDirs));
grpOdrDecode = cell(1,1,length(fileDirs));

%%
for ani = 1:length(fileDirs)
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
% %     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
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
    smi(ani,:) = mlb.smi;
    dPrm(ani,:) = mlb.dPrime;
    ri(ani,:) = mlb.ri;
    for op = 1:mlb.seqLength
        smiByOP(ani,:,1) = reshape(mlb.smiByPos', [1,numel(mlb.smiByPos)]);
        smiByOP(ani,:,2) = reshape(mlb.smiByOdr', [1,numel(mlb.smiByOdr)]);
        dPrmByOP(ani,:,1) = reshape(mlb.dPrimeByPos', [1,numel(mlb.dPrimeByPos)]);
        dPrmByOP(ani,:,2) = reshape(mlb.dPrimeByOdr', [1,numel(mlb.dPrimeByOdr)]);
        riByOP(ani,:,1) = reshape(mlb.riByPos', [1,numel(mlb.riByPos)]);
        riByOP(ani,:,2) = reshape(mlb.riByOdr', [1,numel(mlb.riByOdr)]);
    end
    %% Process Neural Data
    mlb.SetLikes_SubSample;
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
    timeAccuracy = nan(size(timeDecode,1), length(timePoints), length(trlOdrs));
    timePosts = nan(size(tempPosts,1), size(tempPosts,2), length(trlOdrs));
    for odr = 1:length(trlOdrs)
        timePosts(:,:,odr) = mean(tempPosts(:,:,trlOdrVect==trlOdrs(odr)),3,'omitnan'); 
        for t = 1:length(timePoints)
            timeAccuracy(:,t,odr) = mean(timeDecode(:,trlOdrVect==trlOdrs(odr))==timePoints(t),2, 'omitnan');
        end
    end
    grpPosts{ani} = timePosts;
    grpTimeDecode{ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions), length(positions));
    for pos = 1:length(positions)
        for p = 1:length(positions)
            posAccuracy(:,p,pos) = mean(posDecode(:,trlPosVect==positions(pos))==positions(p),2, 'omitnan');
        end
    end
    grpPosDecode{ani} = posAccuracy; 
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(tempPosts, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors), length(odors));
    for odr = 1:length(odors)
        for o = 1:length(odors)
            odrAccuracy(:,o,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(o),2, 'omitnan');
        end
    end
    grpOdrDecode{ani} = odrAccuracy;     
    
    if ani==1
        iscBS_DecodeTime = mlb.likeTimeVect;
        iscBS_ObserveTime = mlb.obsvTimeVect;
    end
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%%
figure;
sp = nan(size(trlOdrs));
spp = nan(size(trlOdrs));
for o = 1:length(trlOdrs)
    subplot(4,length(trlOdrs),[sub2ind([length(trlOdrs),4],o,1), sub2ind([length(trlOdrs),4],o,2)]);
    imagesc(mlb.obsvTimeVect, 1:length(mlb.likeTimeVect), mean(cell2mat(cellfun(@(a)a(:,:,o)', grpPosts, 'uniformoutput', 0)),3), [0 0.015]);
    set(gca, 'yticklabel', [], 'ydir', 'normal');
    title(trlOdrs(o));
    sp(o) = subplot(4,length(trlOdrs),sub2ind([length(trlOdrs),4],o,3));
    hold on;
    grpDecode = cell2mat(cellfun(@(a)a(:,:,o), grpOdrDecode, 'uniformoutput',0));
    for op = 1:length(trlOdrs)
        [~,odrPos] = find(mlb.odrSeqs==trlOdrs(op));
        meanDec = mean(grpDecode(:,op,:),3);
        semDec = mlb.SEMcalc(grpDecode(:,op,:),0,3);
        if trlOdrs(op)>=10
            plot(mlb.obsvTimeVect, meanDec, 'color', mlb.PositionColors(odrPos,:), 'linestyle', '--');
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(meanDec+semDec)', flipud(meanDec-semDec)'], 'edgecolor', mlb.PositionColors(odrPos,:),...
            'facecolor', mlb.PositionColors(odrPos,:), 'facealpha', 0.5, 'linestyle', '--');
        else
            plot(mlb.obsvTimeVect, meanDec, 'color', mlb.PositionColors(odrPos,:));
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(meanDec+semDec)', flipud(meanDec-semDec)'], 'edgecolor', mlb.PositionColors(odrPos,:),...
            'facecolor', mlb.PositionColors(odrPos,:), 'facealpha', 0.5);
        end
    end
    axis tight;
    spp(o) = subplot(4,length(trlOdrs), sub2ind([length(trlOdrs),4],o,4));
    [~,odrPos] = find(mlb.odrSeqs==trlOdrs(o));
    grpDecodePos = cell2mat(cellfun(@(a)a(:,:,odrPos), grpPosDecode, 'uniformoutput', 0));
    for p = 1:mlb.seqLength
        meanDec = mean(grpDecodePos(:,p,:),3);
        semDec = mlb.SEMcalc(grpDecodePos(:,p,:),0,3);
        plot(mlb.obsvTimeVect, meanDec, 'color', mlb.PositionColors(p,:));
        hold on;
        patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
            'YData', [(meanDec+semDec)', flipud(meanDec-semDec)'], 'edgecolor', mlb.PositionColors(p,:),...
            'facecolor', mlb.PositionColors(p,:), 'facealpha', 0.5);
    end
    axis tight;
end
linkaxes(sp, 'xy');
linkaxes(spp, 'xy');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("Dual-list 1window: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%%
posOdrDiff = cell(size(fileDirs));
aniOdrDiff = cell(size(fileDirs));
for ani = 1:length(fileDirs)
    tempPosOdrDiff = nan(size(grpOdrDecode{1},1),mlb.seqLength, 2);
    for pos = 1:mlb.seqLength
        posOdrs = mlb.odrSeqs(:,pos);
        tempOdrDecode(:,1,1) = grpOdrDecode{ani}(:,trlOdrs==posOdrs(1), trlOdrs==posOdrs(1));
        tempOdrDecode(:,2,1) = grpOdrDecode{ani}(:,trlOdrs==posOdrs(1), trlOdrs==posOdrs(2));
        tempOdrDecode(:,1,2) = grpOdrDecode{ani}(:,trlOdrs==posOdrs(2), trlOdrs==posOdrs(1));
        tempOdrDecode(:,2,2) = grpOdrDecode{ani}(:,trlOdrs==posOdrs(2), trlOdrs==posOdrs(2));
        tempPosOdrDiff(:,pos,:) = diff(tempOdrDecode,[],2);
    end
    posOdrDiff{ani} = tempPosOdrDiff;
    aniOdrDiff{ani} = mean(tempPosOdrDiff,2);
end

figure; 
sps = nan(1,mlb.seqLength+1);
for op = 1:mlb.seqLength
    sps(op) = subplot(1,mlb.seqLength+1,op);
    hold(sps(op),'on');
    tempOPs = cell2mat(cellfun(@(a)a(:,op,:), posOdrDiff, 'uniformoutput',0));
    tempMean = mean(tempOPs,2);
    tempVar = mlb.SEMcalc(tempOPs,0,2);
    
    plot(sps(op), mlb.obsvTimeVect, tempMean(:,:,1), 'color', mlb.PositionColors(op,:));
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(tempMean(:,:,1)+tempVar(:,:,1))', flipud(tempMean(:,:,1)-tempVar(:,:,1))'],...
        'edgecolor', mlb.PositionColors(op,:), 'facecolor', mlb.PositionColors(op,:), 'facealpha', 0.5);
    plot(sps(op), mlb.obsvTimeVect, tempMean(:,:,2), 'color', mlb.PositionColors(op,:), 'linestyle', '--');
    patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
        'YData', [(tempMean(:,:,2)+tempVar(:,:,2))', flipud(tempMean(:,:,2)-tempVar(:,:,2))'],...
        'edgecolor', mlb.PositionColors(op,:), 'facecolor', mlb.PositionColors(op,:), 'facealpha', 0.5,...
        'linestyle','--');    
end
grpPosOdrDiff = cell2mat(posOdrDiff);
% grpPosOdrDiff = cell2mat(aniOdrDiff);
gpodMean = median(grpPosOdrDiff,2);
gpodVar = mlb.SEMcalc(grpPosOdrDiff,0,2);
sps(end) = subplot(1,mlb.seqLength+1, mlb.seqLength+1);
hold(sps(end), 'on');
plot(sps(end), mlb.obsvTimeVect, gpodMean(:,:,1), '-k');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(gpodMean(:,:,1)+gpodVar(:,:,1))', flipud(gpodMean(:,:,1)-gpodVar(:,:,1))'],...
    'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.5);
plot(sps(end), mlb.obsvTimeVect, gpodMean(:,:,2), '--k');
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(gpodMean(:,:,2)+gpodVar(:,:,2))', flipud(gpodMean(:,:,2)-gpodVar(:,:,2))'],...
    'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.5,...
    'linestyle','--');
linkaxes(sps, 'xy');
for sp = 1:length(sps)
    plot(sps(sp), [0 0], get(gca,'ylim'), '--k', 'linewidth', 2);
    plot(sps(sp), repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
    plot(sps(sp), repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
    plot(sps(sp), [min(mlb.obsvTimeVect) max(mlb.obsvTimeVect)], [0 0], '-k', 'linewidth', 2);
end
annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', sprintf("Dual-list 1window: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');