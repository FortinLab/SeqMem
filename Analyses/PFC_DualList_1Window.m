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
numPerms = 10;
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
        for p = 1:length(odors)
            odrAccuracy(:,p,odr) = mean(odrDecode(:,trlOdrVect==odors(odr))==odors(p),2, 'omitnan');
        end
    end
    grpOdrDecode{ani} = odrAccuracy;     
    
    if ani==1
        iscBS_DecodeTime = mlb.likeTimeVect;
        iscBS_ObserveTime = mlb.obsvTimeVect;
    end
end

%%
figure;
for o = 1:length(trlOdrs)
    subplot(length(trlOdrs),1,o);
    imagesc(1:length(mlb.likeTimeVect), mlb.obsvTimeVect, mean(cell2mat(cellfun(@(a)a(:,:,o), grpPosts, 'uniformoutput', 0)),3), [0 0.025]);
end
