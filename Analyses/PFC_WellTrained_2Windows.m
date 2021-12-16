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
trlWindow = {[-500 500]; [-200 500]};
alignment = [{'PokeIn'}, {'PokeOut'}];
lfpWindow = [16 32];
numPerms = 100;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
% OSC Decoding
tMatPosts = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
tMatTimeDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
tMatOdrDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
tMatPosDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
% TAO Decoding
taoPosts = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
taoTimeDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
taoOdrDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
taoPosDecode = cell(numel(mlb.odrSeqs), mlb.seqLength, length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.numPerms = numPerms;
    mlb.ssProportion = ssProportion;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    
    %% Decode OSC via FISC 
    mlb.SetLikes_FISC;
    mlb.SetObserves_Session;
    mlb.Process_Observes;
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Time')));
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Odor')));
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,strcmp(mlb.trialIDs, 'Position')));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    for o = 1:numel(mlb.odrSeqs)
        for p = 1:mlb.seqLength
            opLog = [mlb.trialInfo.Odor]==o & [mlb.trialInfo.Position]==p & [mlb.trialInfo.Performance]==1;
            if sum(opLog)>=1
                if o==p     %only works with well-trained sessions... need to adjust if running with Novel or Dual List sessions
%                     opLog = false(size(opLog));
%                     opLog(mlb.fiscTrials(o,:)) = true;
                    opLog(mlb.fiscTrials(o,:)) = false;
                end
                tMatPosts{o,p,ani} = mean(mlb.post{1}(:,:,opLog),3);
                % Temporal Accuracy
                tempTimeDecode = timeDecode(:,opLog);
                tempTimeAccuracy = nan(size(timeDecode,1), length(timePoints));
                for t = 1:size(timePoints,1)
                    tempTimeAccuracy(:,t) = sum(tempTimeDecode==timePoints(t),2);
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
                tempOdrAccuracy = nan(size(odrDecode,1), numel(mlb.odrSeqs));
                for oo = 1:numel(mlb.odrSeqs)
                    tempOdrAccuracy(:,oo) = mean(tempOdrDecode==mlb.odrSeqs(oo),2);
                end
                tMatOdrDecode{o,p,ani} = tempOdrAccuracy;
            end
        end
    end                    
end

%%
piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==max(mlb.likeTimeVect)-(mlb.windows{end}(2)/1000));
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;

figure;
for sp = 1:numel(mlb.odrSeqs)*mlb.seqLength
    [spR, spC] = ind2sub([numel(mlb.odrSeqs), mlb.seqLength], sp);
    tempPost = mean(cell2mat(tMatPosts(spR,spC,:)),3,'omitnan');
    if isempty(tempPost)
        continue;
    else
        subplot(numel(mlb.odrSeqs), mlb.seqLength, sub2ind([numel(mlb.odrSeqs), mlb.seqLength], spC, spR));
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
    'String', sprintf("FISC Decode ISC Group: binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s; %.0fms:%.0fms to %s), NumPerms = %.0f, Subsample Proportion = %.01f, Subsample Type = %.0f, BayesType = %.0f",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, trlWindow{2}(1), trlWindow{2}(2), alignment{2}, numPerms, ssProportion, ssType, bayesType),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');