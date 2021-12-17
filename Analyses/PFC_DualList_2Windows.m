clear all; %#ok<CLALL>
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

binSize = 200;
dsRate = 50;
trlWindow = {[-500 500]; [-500 500]};
alignment = [{'PokeIn'}, {'PokeOut'}];
lfpWindow = [16 32];
numPerms = 10;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
grpPosts = cell(1,1,length(fileDirs));
grpTimeDecode = cell(1,1,length(fileDirs));
grpOdrDecode = cell(1,1,length(fileDirs));
grpPosDecode = cell(1,1,length(fileDirs));

%% Analyze Data
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
    mlb.SetLikes_SubSample;
    mlb.Process_Observes;
    %%
    timePoints = sort(unique(mlb.likeTimeVect));
    odors = sort(mlb.odrSeqs(:));
    timeAccuracy = cell(1,mlb.numPerms);
    odorAccuracy = cell(1,mlb.numPerms);
    posAccuracy = cell(1,mlb.numPerms);
    for rep = 1:mlb.numPerms
        tempTimeDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Time')));
        tempTimeAccuracy = nan(length(mlb.likeTimeVect), length(mlb.obsvTimeVect), length(odors));
        for t = 1:length(timePoints)
        end
        
        tempOdrDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Odor')));
        tempOdrAccuracy = nan(length(mlb.obsvTimeVect), length(odors), length(odors));
        for o = 1:length(odors)
            for oo = 1:length(odors)
            end
        end
        
        tempPosDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Position')));
        tempPosAccuracy = nan(length(mlb.obsvTimeVect), mlb.seqLength, length(odors));
        for o = 1:length(odors)
            for p = 1:mlb.seqLength
            end
        end
    end
end

    