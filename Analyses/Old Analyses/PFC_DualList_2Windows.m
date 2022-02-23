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
ssType = 0; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
grpPosts = cell(1,1,length(fileDirs));
grpTimeAccuracy = cell(1,1,length(fileDirs));
grpOdrAccuracy = cell(1,1,length(fileDirs));
grpPosAccuracy = cell(1,1,length(fileDirs));

%% Analyze Data
for ani = 1:length(fileDirs)
    %%
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
    %% Organize posteriors and decodings based on Odor-In-Position
    timePoints = sort(unique(mlb.likeTimeVect));
    odors = sort(mlb.odrSeqs(:));
    posts = cell(1,mlb.numPerms);
    timeAccuracy = cell(1,mlb.numPerms);
    odorAccuracy = cell(1,mlb.numPerms);
    posAccuracy = cell(1,mlb.numPerms);
    for rep = 1:mlb.numPerms
        tempTimeDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Time')));
        tempOdrDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Odor')));
        tempPosDecode = mlb.DecodeBayesPost(mlb.post{rep}, mlb.decodeIDvects{rep}(:,strcmp(mlb.trialIDs, 'Position')));
        
        tempPost = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), length(odors));
        tempTimeAccuracy = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), length(odors));
        tempOdrAccuracy = nan(length(mlb.obsvTimeVect), length(odors), length(odors));
        tempPosAccuracy = nan(length(mlb.obsvTimeVect), mlb.seqLength, length(odors));
        for o = 1:length(odors)
            odrLog = [mlb.trialInfo(mlb.postTrlIDs{rep}).Odor]==odors(o);
            tempPost(:,:,o) = mean(mlb.post{rep}(:,:,odrLog),3,'omitnan');
            for t = 1:length(timePoints)
                tempTimeAccuracy(:,t,o) = mean(tempTimeDecode(:,odrLog)==timePoints(t),2);                
            end
            for oo = 1:length(odors)
                tempOdrAccuracy(:,oo,o) = mean(tempOdrDecode(:,odrLog)==odors(oo),2);
            end
            for p = 1:mlb.seqLength
                tempPosAccuracy(:,p,o) = mean(tempPosDecode(:,odrLog)==p,2);
            end
        end
        posts{rep} = tempPost;
        timeAccuracy{rep} = tempTimeAccuracy;
        odorAccuracy{rep} = tempOdrAccuracy;
        posAccuracy{rep} = tempPosAccuracy;
    end
    %%
    aniPost = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), length(odors));
    aniTimeAccuracy = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), length(odors));
    aniOdrAccuracy = nan(length(mlb.obsvTimeVect), length(odors), length(odors));
    aniPosAccuracy = nan(length(mlb.obsvTimeVect), mlb.seqLength, length(odors));
    for o = 1:length(odors)
        aniPost(:,:,o) = mean(cell2mat(reshape(cellfun(@(a)a(:,:,o),posts, 'uniformoutput', 0), [1,1,mlb.numPerms])),3,'omitnan');
        aniTimeAccuracy(:,:,o) = mean(cell2mat(reshape(cellfun(@(a)a(:,:,o),timeAccuracy, 'uniformoutput', 0), [1,1,mlb.numPerms])),3,'omitnan');
        aniOdrAccuracy(:,:,o) = mean(cell2mat(reshape(cellfun(@(a)a(:,:,o),odorAccuracy, 'uniformoutput', 0), [1,1,mlb.numPerms])),3,'omitnan');
        aniPosAccuracy(:,:,o) = mean(cell2mat(reshape(cellfun(@(a)a(:,:,o),posAccuracy, 'uniformoutput', 0), [1,1,mlb.numPerms])),3,'omitnan');
    end
    grpPosts{ani} = aniPost;
    grpTimeAccuracy{ani} = aniTimeAccuracy;
    grpOdrAccuracy{ani} = aniOdrAccuracy;
    grpPosAccuracy{ani} = aniPosAccuracy;        
end
%% Plot everything!
posts = cell(length(odors),1);
for o = 1:length(odors)
    posts{o} = mean(cell2mat(cellfun(@(a)a(:,:,o),grpPosts, 'uniformoutput', 0)),3,'omitnan');
end
figure;
subplot(5,5,[2:
imagesc(cell2mat(posts));
set(gca,'ydir', 'normal');

    