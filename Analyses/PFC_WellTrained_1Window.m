% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Session096'}];

binSize = 100;
dsRate = 20;
trlWindow = {[-800 1200]};
alignment = {'PokeIn'};
lfpWindow = [16 32];
numPerms = 10;
ssProportion = 0.7;
ssType = 0; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

%% Create Output Variables
% fiscLOO analysis
fiscLOO_Posts = cell(size(fileDirs));
fiscLOO_Decodes = cell(3,1,length(fileDirs));
% fiscISC analysis
fiscISC_Posts = cell(size(fileDirs));
fiscISC_Decodes = cell(3,1,length(fileDirs));
% iscBS analysis
iscBS_Posts = cell(size(fileDirs));
iscBS_Decodes = cell(3,1,length(fileDirs));
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
    
    %% Decode FISC via Leave-1-Out
%     trialIDs = [{'Time'}, {'Window'}, {'Position'}, {'Odor'}];
    mlb.SetLikes_FISC;
    mlb.Process_LikelyLOO;
    fiscLOO_Posts{ani} = mlb.post;
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints));
    for t = 1:length(timePoints)
        timeAccuracy(:,t) = mean(timeDecode==timePoints(t),2, 'omitnan');
    end
    fiscLOO_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions));
    for p = 1:length(positions)
        posAccuracy(:,p) = mean(posDecode==positions(p),2, 'omitnan');
    end
    fiscLOO_Decodes{2,1,ani} = posAccuracy;
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors));
    for o = 1:length(odors)
        odrAccuracy(:,o) = mean(odrDecode==odors(o),2, 'omitnan');
    end
    fiscLOO_Decodes{3,1,ani} = odrAccuracy;    
    
    %% Decode ISC via FISC 
    mlb.Process_Observes;
    fiscISC_Posts{ani} = mlb.post;
    mlb.trialInfo(mlb.postTrlIDs{1})
    trlOdrVect = [mlb.trialInfo(reshape(mlb.postTrlIDs{1}, [1,numel(mlb.postTrlIDs{1})])).Odor];
    trlOdrs = unique(trlOdrVect);
    trlPosVect = [mlb.trialInfo(reshape(mlb.postTrlIDs{1}, [1,numel(mlb.postTrlIDs{1})])).Position];
    trlPoss = unique(trlPosVect);
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints), length(trlPoss));
    for pos = 1:length(trlPoss)
        for t = 1:length(timePoints)
            timeAccuracy(:,t,pos) = mean(timeDecode(:,trlPosVect==trlPoss(pos))==timePoints(t),2, 'omitnan');
        end
    end
    fiscISC_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions), length(positions));
    for pos = 1:length(positions)
        for p = 1:length(positions)
            posAccuracy(:,p,pos) = mean(posDecode(:,trlPosVect==positions(pos))==positions(p),2, 'omitnan');
        end
    end
    fiscISC_Decodes{2,1,ani} = posAccuracy; 
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(mlb.post{1}, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors), length(odors));
    for odr = 1:length(odors)
        for o = 1:length(odors)
            odrAccuracy(:,o,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(o),2, 'omitnan');
        end
    end
    fiscISC_Decodes{3,1,ani} = odrAccuracy;    
    
    %% Decode ISC via Sub-Sampling
    mlb.SetLikes_SubSample;
    mlb.Process_Observes;
    iscBS_Posts{ani} = cell2mat(reshape(mlb.post, [1,1,numel(mlb.post)]));
    postTrlIDs = cell2mat(reshape(mlb.postTrlIDs, [1,1,numel(mlb.postTrlIDs)]));
    
    trlOdrVect = [mlb.trialInfo(postTrlIDs).Odor];
    trlOdrs = unique(trlOdrVect);
    trlPosVect = [mlb.trialInfo(postTrlIDs).Position];
    trlPoss = unique(trlPosVect);
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(iscBS_Posts{ani}, mlb.decodeIDvects{1}(:,1));
    timePoints = unique(mlb.decodeIDvects{1}(:,1));
    timeAccuracy = nan(size(timeDecode,1), length(timePoints), length(trlPoss));
    for pos = 1:length(trlPoss)
        for t = 1:length(timePoints)
            timeAccuracy(:,t,pos) = mean(timeDecode(:,trlPosVect==trlPoss(pos))==timePoints(t),2, 'omitnan');
        end
    end
    iscBS_Decodes{1,1,ani} = timeAccuracy;
    % Decode Position
    posDecode = mlb.DecodeBayesPost(iscBS_Posts{ani}, mlb.decodeIDvects{1}(:,3));
    positions = unique(mlb.decodeIDvects{1}(:,3));
    posAccuracy = nan(size(posDecode,1), length(positions), length(positions));
    for pos = 1:length(positions)
        for p = 1:length(positions)
            posAccuracy(:,p,pos) = mean(posDecode(:,trlPosVect==positions(pos))==positions(p),2, 'omitnan');
        end
    end
    iscBS_Decodes{2,1,ani} = posAccuracy; 
    % Decode Odor
    odrDecode = mlb.DecodeBayesPost(iscBS_Posts{ani}, mlb.decodeIDvects{1}(:,4));
    odors = unique(mlb.decodeIDvects{1}(:,4));
    odrAccuracy = nan(size(odrDecode,1),length(odors), length(odors));
    for odr = 1:length(odors)
        for o = 1:length(odors)
            odrAccuracy(:,o,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(o),2, 'omitnan');
        end
    end
    iscBS_Decodes{3,1,ani} = odrAccuracy;    
end

figure;
annotation(gcf,'textbox', [0.1 0.95 0.8 0.05],...
    'String', sprintf('Dual List Group:\n binSize = %.0fms, dsRate = %.0fms, SartWindow = (%.0fms:%.0fms to pokeIn), EndWindow = (%.0fms:%.0fms to pokeOut), Likelihood Porportion = %.02f, NumPerms = %.0f',...
    binSize, dsRate, stWin(1), stWin(2), ndWin(1), ndWin(2), templateProportion, numPerms),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

