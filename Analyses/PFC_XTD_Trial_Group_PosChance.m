% PFC_XTD_Trial_Group_PosChance
%%
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
% trlWindow = {[-1000 2000]}; %pi1
% trlWindow = {[-1500 2000]}; %pi2
% alignment = {'PokeIn'};
% trlWindow = {[-1800 2000]}; %po1
trlWindow = {[-2000 1500]}; %po2
alignment = {'PokeOut'};
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

postCLim = [0 0.05];
decodeCLim = [0 0.2];


%%
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb = MLB_SM(fileDirs{ani});
    % Create Analysis Variables
    if ani == 1 
        % Behavior Variables
        fiscPokeOutLat = cell(length(fileDirs),1);
        fiscRwdDelivLat = cell(length(fileDirs),1);
        fiscNxtTrlStrtLag = cell(length(fileDirs),1);
        smi = nan(length(fileDirs),1);
        dPrm = nan(length(fileDirs),1);
        ri = nan(length(fileDirs),1);
        smiByOP = nan(length(fileDirs),mlb.seqLength,2);
        dPrmByOP = nan(length(fileDirs),mlb.seqLength,2);
        riByOP = nan(length(fileDirs),mlb.seqLength,2);
        % Posteriors
        realPost = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
        realTrlD = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
        realSsnD = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
        realDecode = cell(mlb.seqLength,1,length(fileDirs));
        realDecodeD = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
    end
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.bayesType = bayesType;
    
    %% Extract Behavioral Variables
    if strcmp(alignment{1}, 'PokeIn')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
    elseif strcmp(alignment{1}, 'PokeOut')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeInIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
        fiscNxtTrlStrtLag{ani} = ([mlb.trialInfo(mlb.fiscTrials(1:end-1,:)+1).PokeInIndex]-[mlb.trialInfo(mlb.fiscTrials(1:end-1,:)).PokeOutIndex])'/1000;
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
%     mlb.SetLikes_FISC;
    mlb.SetLikes_ISC;
    %% Comment in if running L1O
    mlb.Process_IterativeLikelyL1O;
    trlIDs = mlb.postTrlIDs(:);
    trlIDlog = ~isnan(trlIDs);
    posts = permute(mlb.post(:),[3,2,1]);
    posts = posts(:,:,trlIDlog);
    unq = unique(cell2mat(posts));
    minProb = unq(2);
    maxProb = unq(end-1);
    trlIDs = trlIDs(trlIDlog);
    trlPosVect = [mlb.trialInfo(trlIDs).Position];
    %% Decode & calculate d' for session
    tempDecode = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), size(posts,3));
    for trl = 1:size(posts,3)
        [~,tempDecode(:,:,trl)] = max(posts{trl},[],3);
    end
    for pos = 1:mlb.seqLength
        realDecode{pos,1,ani} = tempDecode(:,:,trlPosVect==pos);
        tempNonPos = tempDecode(:,:,trlPosVect~=pos);
        for odr = 1:mlb.seqLength
            tempHits = sum(realDecode{pos,1,ani}==odr,3);
            tempFNs = sum(realDecode{pos,1,ani}~=odr,3);
            tempFAs = sum(tempNonPos==odr,3);
            tempCRs = sum(tempNonPos~=odr,3);
            realDecodeD{pos,odr,ani} = arrayfun(@(a,b,c,d)mlb.CalculateDprime([a,b;c,d]), tempHits, tempFNs, tempFAs, tempCRs);
        end
    end
%     figure; for pos = 1:mlb.seqLength; for p = 1:mlb.seqLength; subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], pos,p)); imagesc(zscore(realDecodeD{pos,p,ani}, 0, 'all'), [-2 2]); set(gca,'ydir','normal'); end; end
%     drawnow;           
            
    %% Comment in if running ISC-FISC or vice versa
%     mlb.Process_IterativeObserves;
%     unq = unique(mlb.post);
%     minProb = unq(2);
%     maxProb = unq(end-1);
%     trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
%     trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
%     trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
%     iscLog = trlOdrVect==trlPosVect & trlPerfVect==1;
%     posts = cell(1,1,size(mlb.post,3));
%     for trl = 1:size(mlb.post,4)
%         posts{trl} = mlb.post(:,:,:,trl);
%     end
%     posts = posts(:,:,iscLog);
%     trlIDs = mlb.postTrlIDs(iscLog);
%     trlPosVect = [mlb.trialInfo(mlb.postTrlIDs(iscLog)).Position];    
    %% Calculate d' values from posteriors
    for pos = 1:mlb.seqLength
        tempPosts = posts(:,:,trlPosVect==pos);
        for p = 1:mlb.seqLength
            realPost{pos,p,ani} = cell2mat(cellfun(@(a){a(:,:,p)}, tempPosts));
            hits = cell2mat(cellfun(@(a){a(:,:,p)}, tempPosts));
            fas = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect~=pos)));
            meanHR = mean(hits,3);
            meanHR(meanHR==0) = minProb;
            meanHR(meanHR==1) = maxProb;
            meanFAR = mean(fas,3);
            meanFAR(meanFAR==0) = minProb;
            meanFAR(meanFAR==1) = maxProb;
            realSsnD{pos,p,ani} = arrayfun(@(a,b)norminv(a)-norminv(b), meanHR, meanFAR);
            tempTrlD = nan(size(hits));
            for trl = 1:size(hits,3)
                tempHR = hits(:,:,trl);
                tempHR(tempHR==0) = minProb;
                tempHR(tempHR==1) = maxProb;
                tempTrlD(:,:,trl) = arrayfun(@(a,b)norminv(a)-norminv(b), tempHR, meanFAR);
            end
            realTrlD{pos,p,ani} = tempTrlD;
        end
    end    
end

%% Collapse trial data across animals
grpTrlPost = cell(mlb.seqLength);
grpTrlD = cell(mlb.seqLength);
grpAniD = cell(mlb.seqLength);
for pos = 1:mlb.seqLength
    for p = 1:mlb.seqLength
        grpTrlPost{pos,p} = cell2mat(realPost(pos,p,:));
        grpTrlD{pos,p} = cell2mat(realTrlD(pos,p,:));
        grpAniD{pos,p} = cell2mat(realSsnD(pos,p,:));
    end
end
    
pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%% Process Chance Data    
fprintf('Processing chance now...\n');
chancePostPerms = cell(mlb.seqLength, mlb.seqLength, numChancePerms);
chanceSsnDPerms = cell(mlb.seqLength, mlb.seqLength, numChancePerms);
chanceTrlDPerms = cell(mlb.seqLength, mlb.seqLength, numChancePerms);
chancePostLagPerms = cell(1, length(mlb.lagVect), numChancePerms);
chancePostLag23Perms = cell(1, length(mlb.lagVect), numChancePerms);
chancePostLag123Perms = cell(1, length(mlb.lagVect), numChancePerms);
chancePostLag234Perms = cell(1, length(mlb.lagVect), numChancePerms);
for perm = 1:numChancePerms
    tempChancePost = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
    tempChanceSsnD = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
    tempChanceTrlD = cell(mlb.seqLength,mlb.seqLength,length(fileDirs));
    for ani = 1:length(fileDirs)
        % Load and setup the MLB object
        mlb = MLB_SM(fileDirs{ani});
        mlb.binSize = binSize;
        mlb.dsRate = dsRate;
        mlb.windows = trlWindow;
        mlb.alignments = alignment;
        mlb.bayesType = bayesType;
        % Set the likelihoods
        mlb.SetLikes_ISC;
        % Shuffle the trials around to destroy position information for the likelihoods
        tempTrlIDmtx = permute(cellfun(@(a)a(1,end),mlb.likeTrlIDs), [1,3,2]);
        chancePerm = sortrows([randperm(numel(tempTrlIDmtx(~isnan(tempTrlIDmtx))))',find(~isnan(tempTrlIDmtx))]);
        mlb.likeTrlSpikes(~isnan(tempTrlIDmtx)) = mlb.likeTrlSpikes(chancePerm(:,2));
        % Process the trials using L1O
        mlb.Process_IterativeLikelyL1O;
        trlIDs = mlb.postTrlIDs(:);
        trlIDlog = ~isnan(trlIDs);
        posts = permute(mlb.post(:),[3,2,1]);
        unq = unique(cell2mat(posts));
        minProb = unq(2);
        maxProb = unq(end-1);
        posts = posts(:,:,trlIDlog);
        trlIDs = trlIDs(trlIDlog);
        trlPosVect = [mlb.trialInfo(trlIDs).Odor];
        for pos = 1:mlb.seqLength
            tempPosts = posts(:,:,trlPosVect==pos);
            for p = 1:mlb.seqLength
                tempChancePost{pos,p,ani} = cell2mat(cellfun(@(a){a(:,:,p)}, tempPosts));
                hits = cell2mat(cellfun(@(a){a(:,:,p)}, tempPosts));
                fas = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect~=pos)));
                meanHR = mean(hits,3);
                meanHR(meanHR==0) = minProb;
                meanHR(meanHR==1) = maxProb;
                meanFAR = mean(fas,3);
                meanFAR(meanFAR==0) = minProb;
                meanFAR(meanFAR==1) = maxProb;
                tempChanceSsnD{pos,p,ani} = arrayfun(@(a,b)norminv(a)-norminv(b), meanHR, meanFAR);
                tempTrlD = nan(size(hits));
                for trl = 1:size(hits,3)
                    tempHR = hits(:,:,trl);
                    tempHR(tempHR==0) = minProb;
                    tempHR(tempHR==1) = maxProb;
                    tempTrlD(:,:,trl) = arrayfun(@(a,b)norminv(a)-norminv(b), tempHR, meanFAR);
                end
                tempChanceTrlD{pos,p,ani} = tempTrlD;
            end
        end
    end
    tempGrpTrlD = cell(mlb.seqLength, mlb.seqLength);
    for pos = 1:mlb.seqLength
        for p = 1:mlb.seqLength
            chancePostPerms{pos,p,perm} = mean(cell2mat(tempChancePost(pos,p,:)),3);
            chanceSsnDPerms{pos,p,perm} = mean(cell2mat(tempChanceSsnD(pos,p,:)),3);
            chanceTrlDPerms{pos,p,perm} = mean(cell2mat(tempChanceTrlD(pos,p,:)),3);
            tempGrpTrlD{pos,p} = cell2mat(tempChanceTrlD(pos,p,:));
        end
    end    
    tempChancePostLag = cell(size(mlb.lagVect));
    tempChancePost23Lag = cell(size(mlb.lagVect));
    tempChancePost123Lag = cell(size(mlb.lagVect));
    tempChancePost234Lag = cell(size(mlb.lagVect));
    for pos = 1:mlb.seqLength
        for odr = 1:mlb.seqLength
            tempLagLog = (pos-odr)==mlb.lagVect;
            if pos~=4 && odr~=4
                tempChancePost123Lag{tempLagLog} = cat(3,tempChancePost123Lag{tempLagLog}, tempGrpTrlD{pos,odr});
            end
            if pos~=1 && odr~=1
                tempChancePost234Lag{tempLagLog} = cat(3,tempChancePost234Lag{tempLagLog}, tempGrpTrlD{pos,odr});
            end
            if (pos~=1 && odr~=1) && (pos~=4 && odr~=4)
                tempChancePost23Lag{tempLagLog} = cat(3,tempChancePost23Lag{tempLagLog}, tempGrpTrlD{pos,odr});
            end
            tempChancePostLag{tempLagLog} = cat(3,tempChancePostLag{tempLagLog}, tempGrpTrlD{pos,odr});
        end
    end
    for lag = 1:length(mlb.lagVect)
        chancePostLagPerms{1,lag,perm} = mean(tempChancePostLag{lag},3);
        chancePostLag23Perms{1,lag,perm} = mean(tempChancePost23Lag{lag},3);
        chancePostLag123Perms{1,lag,perm} = mean(tempChancePost123Lag{lag},3);
        chancePostLag234Perms{1,lag,perm} = mean(tempChancePost234Lag{lag},3);
    end
    fprintf('Iteration %i complete\n', perm);
end
chancePost = cell(mlb.seqLength, mlb.seqLength, 2);
chanceSsnD = cell(mlb.seqLength, mlb.seqLength, 2);
chanceTrlD = cell(mlb.seqLength, mlb.seqLength, 2);
for pos = 1:mlb.seqLength
    for p = 1:mlb.seqLength
        chancePost{pos,p,1} = mean(cell2mat(chancePostPerms(pos,p,:)),3);
        chancePost{pos,p,2} = std(cell2mat(chancePostPerms(pos,p,:)),0,3);
        
        chanceSsnD{pos,p,1} = mean(cell2mat(chanceSsnDPerms(pos,p,:)),3);
        chanceSsnD{pos,p,2} = std(cell2mat(chanceSsnDPerms(pos,p,:)),0,3);
        
        chanceTrlD{pos,p,1} = mean(cell2mat(chanceTrlDPerms(pos,p,:)),3);
        chanceTrlD{pos,p,2} = std(cell2mat(chanceTrlDPerms(pos,p,:)),0,3);
    end
end
chanceLag = cell(2,length(mlb.lagVect));
chance23Lag = cell(2,length(mlb.lagVect));
chance123Lag = cell(2,length(mlb.lagVect));
chance234Lag = cell(2,length(mlb.lagVect));
for lag = 1:length(mlb.lagVect)
    chanceLag{1,lag} = mean(cell2mat(chancePostLagPerms(1,lag,:)),3);
    chanceLag{2,lag} = std(cell2mat(chancePostLagPerms(1,lag,:)),0,3);
    chance23Lag{1,lag} = mean(cell2mat(chancePostLag23Perms(1,lag,:)),3);
    chance23Lag{2,lag} = std(cell2mat(chancePostLag23Perms(1,lag,:)),0,3);
    chance123Lag{1,lag} = mean(cell2mat(chancePostLag123Perms(1,lag,:)),3);
    chance123Lag{2,lag} = std(cell2mat(chancePostLag123Perms(1,lag,:)),0,3);
    chance234Lag{1,lag} = mean(cell2mat(chancePostLag234Perms(1,lag,:)),3);
    chance234Lag{2,lag} = std(cell2mat(chancePostLag234Perms(1,lag,:)),0,3);    
end



%% Calculate optimal dynamic model fit
% Create pre-trial dynamic models
if strcmp(alignment{1}, 'PokeIn')
    preTrlPrdLog = mlb.obsvTimeVect<0;
    trlPrdLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<nearestPOtime;
    rlyTrlPrdLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<(nearestPOtime/2);
    latTrlPrdLog = mlb.obsvTimeVect>=(nearestPOtime/2) & mlb.obsvTimeVect<nearestPOtime;
elseif strcmp(alignment{1}, 'PokeOut')
    preTrlPrdLog = mlb.obsvTimeVect>0;
    trlPrdLog = mlb.obsvTimeVect<0 & mlb.obsvTimeVect>=nearestPOtime;
    rlyTrlPrdLog = mlb.obsvTimeVect>=nearestPOtime & mlb.obsvTimeVect<(nearestPOtime/2);
    latTrlPrdLog = mlb.obsvTimeVect>=(nearestPOtime/2) & mlb.obsvTimeVect<0;    
end
preTrlLogMtx = false(length(mlb.obsvTimeVect));
preTrlLogMtx(preTrlPrdLog, preTrlPrdLog) = true;
preTrlDynCap = sum(preTrlPrdLog)-1;
preTrlDynMod = nan(sum(preTrlPrdLog), sum(preTrlPrdLog), preTrlDynCap-1);
preTrlDynMod(:,:,1) = eye(sum(preTrlPrdLog));
for d = 1:preTrlDynCap-1
    preTrlDynMod(:,:,d+1) = triu(true(sum(preTrlPrdLog)), d*-1) & tril(true(sum(preTrlPrdLog)), d);
end
% Create trial dynamic models
trlLogMtx = false(length(mlb.obsvTimeVect));
trlLogMtx(trlPrdLog,trlPrdLog) = true;
trlDynCap = sum(trlPrdLog)-1;
trlDynMod = nan(sum(trlPrdLog), sum(trlPrdLog), trlDynCap-1);
trlDynMod(:,:,1) = eye(sum(trlPrdLog));
for d = 1:trlDynCap-1
    trlDynMod(:,:,d+1) = triu(true(sum(trlPrdLog)), d*-1) & tril(true(sum(trlPrdLog)), d);
end
% Create Early Trial Dynamic Models
rlyTrlLogMtx = false(length(mlb.obsvTimeVect));
rlyTrlLogMtx(rlyTrlPrdLog,rlyTrlPrdLog) = true;
rlyTrlDynCap = sum(rlyTrlPrdLog)-1;
rlyTrlDynMod = nan(sum(rlyTrlPrdLog), sum(rlyTrlPrdLog), rlyTrlDynCap-1);
rlyTrlDynMod(:,:,1) = eye(sum(rlyTrlPrdLog));
for d = 1:rlyTrlDynCap-1
    rlyTrlDynMod(:,:,d+1) = triu(true(sum(rlyTrlPrdLog)), d*-1) & tril(true(sum(rlyTrlPrdLog)), d);
end
% Create Late Trial Dynamic Models
latTrlLogMtx = false(length(mlb.obsvTimeVect));
latTrlLogMtx(latTrlPrdLog,latTrlPrdLog) = true;
latTrlDynCap = sum(latTrlPrdLog)-1;
latTrlDynMod = nan(sum(latTrlPrdLog), sum(latTrlPrdLog), latTrlDynCap-1);
latTrlDynMod(:,:,1) = eye(sum(latTrlPrdLog));
for d = 1:latTrlDynCap-1
    latTrlDynMod(:,:,d+1) = triu(true(sum(latTrlPrdLog)), d*-1) & tril(true(sum(latTrlPrdLog)), d);
end


% Calculate fit for time lag invariance across trials
%%% maybe?
trialDynFit_PreTrial = cell(mlb.seqLength,1);
trialDynFit_Trial = cell(mlb.seqLength,1);
trialDynFit_RlyTrial = cell(mlb.seqLength,1);
trialDynFit_LatTrial = cell(mlb.seqLength,1);
trialDynFit_Latencies = cell(mlb.seqLength,1);

aniDynFit_PreTrial = cell(mlb.seqLength,1);
aniDynFit_Trial = cell(mlb.seqLength,1);
for pos = 1:mlb.seqLength
    % Calculate ideal dynamic fit across trials
    curTrlDecode = grpTrlD{pos,pos};
%     curTrlDecode(curTrlDecode<0) = 0;                                               % Comment out if negative decodability is desired
    tempPreTrlFit = nan(size(curTrlDecode,3), size(preTrlDynMod,3));
    tempTrlFit = nan(size(curTrlDecode,3), size(trlDynMod,3));
    tempRlyTrlFit = nan(size(curTrlDecode,3), size(rlyTrlDynMod,3));
    tempLatTrlFit = nan(size(curTrlDecode,3), size(rlyTrlDynMod,3));
    for trl = 1:size(curTrlDecode,3)
        curPreTrl = zscore(curTrlDecode(preTrlPrdLog,preTrlPrdLog,trl),0,'all');
        for pd = 1:size(preTrlDynMod,3)
            tempDynZ = zscore(preTrlDynMod(:,:,pd),0,'all');
            tempPreTrlFit(trl,pd) = pdist([curPreTrl(:)';tempDynZ(:)'], 'cosine');
        end
        curTrl = zscore(curTrlDecode(trlPrdLog,trlPrdLog,trl),0,'all');
        for td = 1:size(trlDynMod,3)
            tempDynZ = zscore(trlDynMod(:,:,td),0,'all');
            tempTrlFit(trl,td) = pdist([curTrl(:)'; tempDynZ(:)'], 'cosine');
        end
        curRlyTrl = zscore(curTrlDecode(rlyTrlPrdLog,rlyTrlPrdLog,trl),0,'all');
        for ed = 1:size(rlyTrlDynMod,3)
            tempDynZ = zscore(rlyTrlDynMod(:,:,ed),0,'all');
            tempRlyTrlFit(trl,ed) = pdist([curRlyTrl(:)'; tempDynZ(:)'], 'cosine');
        end
        curLatTrl = zscore(curTrlDecode(latTrlPrdLog,latTrlPrdLog,trl),0,'all');
        for ld = 1:size(latTrlDynMod,3)
            tempDynZ = zscore(latTrlDynMod(:,:,ld),0,'all');
            tempLatTrlFit(trl,ld) = pdist([curLatTrl(:)'; tempDynZ(:)'], 'cosine');
        end
    end
    trialDynFit_PreTrial{pos} = tempPreTrlFit;
    trialDynFit_Trial{pos} = tempTrlFit;
    trialDynFit_RlyTrial{pos} = tempRlyTrlFit;
    trialDynFit_LatTrial{pos} = tempLatTrlFit;
    [~,preLat] = min(tempPreTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,1) = preLat/(1/dsRate);
%     trialDynFit_Latencies{pos}(:,1) = preLat./preTrlDynCap;
    [~,trlLat] = min(tempTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,2) = trlLat/(1/dsRate);
%     trialDynFit_Latencies{pos}(:,2) = trlLat./trlDynCap;
    [~,rlyLat] = min(tempRlyTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,3) = rlyLat/(1/dsRate);
%     trialDynFit_Latencies{pos}(:,3) = rlyLat./rlyTrlDynCap;
    [~,latLat] = min(tempLatTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,4) = latLat/(1/dsRate);
%     trialDynFit_Latencies{pos}(:,4) = latLat./latTrlDynCap;
end

%% Plot Fits
pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
figure;
for pos = 1:mlb.seqLength
    sp1(pos) = subplot(mlb.seqLength,4,sub2ind([4, mlb.seqLength], 1:2,ones(1,2)*pos));
    tempPreTrl = trialDynFit_PreTrial{pos};
    preMean = mean(tempPreTrl,1,'omitnan');
    preSEM = mlb.SEMcalc(tempPreTrl,0,1);
    preCI = tinv(0.975, size(tempPreTrl,1)-1).*preSEM;
    plot((1:length(preMean))/(1/dsRate),preMean, 'color', 'k', 'linewidth', 1.5);
    hold on;
    patch('XData', [(1:length(preMean))/(1/dsRate), fliplr((1:length(preMean))/(1/dsRate))],...
        'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
        'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
    patch('XData', [(1:length(preMean))/(1/dsRate), fliplr((1:length(preMean))/(1/dsRate))],...
        'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
    
    tempTrl = trialDynFit_Trial{pos};
    trlMean = mean(tempTrl,1,'omitnan');
    trlSEM = mlb.SEMcalc(tempTrl,0,1);
    trlCI = tinv(0.975, size(tempTrl,1)-1).*trlSEM;
    plot((1:length(trlMean))/(1/dsRate),trlMean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
    hold on;    
    patch('XData', [(1:length(trlMean))/(1/dsRate), fliplr((1:length(trlMean))/(1/dsRate))],...
        'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [(1:length(trlMean))/(1/dsRate), fliplr((1:length(trlMean))/(1/dsRate))],...
        'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    xlabel('Persistence duration (ms)');
    ylabel('Distance (cosine similarity)');
    
    subplot(mlb.seqLength,4,sub2ind([4, mlb.seqLength], 3,pos));
    tempLats = trialDynFit_Latencies{pos};%./dsRate;
    bar([mean(tempLats(:,1)), mean(tempLats(:,2))]);
    hold on;
    swarmchart(ones(size(tempLats,1),1),tempLats(:,1), 'ok', 'filled', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.2);
    swarmchart(ones(size(tempLats,1),1)+1, tempLats(:,2), 'ok', 'filled', 'markerfacecolor', mlb.PositionColors(pos,:), 'markeredgecolor', 'none', 'markerfacealpha', 0.2);
    [h,p,ci,stats] = ttest(tempLats(:,1), tempLats(:,2));
    if p<0.05
        title(sprintf('t = %.02f; p = %.02i', stats.tstat, p), 'color','r');
    else
        title(sprintf('t = %.02f; p = %.02f', stats.tstat, p));
    end
        set(gca, 'xtick', [1 2], 'xticklabel', [{'Pre-Trial'}, {'Trial'}]);

    
    subplot(mlb.seqLength,4,sub2ind([4, mlb.seqLength], 4,pos));
    corrScatPlot(tempLats(:,1), tempLats(:,2), 'Pre-Trial', 'Trial', []);    
end
linkaxes(sp1, 'xy');

figure;
subplot(1,4,sub2ind([4, 1], 1:2,ones(1,2)));
tempPreTrl = cell2mat(trialDynFit_PreTrial);
% tempPreTrl = cell2mat(trialDynFit_PreTrial(1:3));
% tempPreTrl = cell2mat(trialDynFit_PreTrial(2:end));
preMean = mean(tempPreTrl,1,'omitnan');
preSEM = mlb.SEMcalc(tempPreTrl,0,1);
preCI = tinv(0.975, size(tempPreTrl,1)-1).*preSEM;
plot((1:length(preMean))/(1/dsRate),preMean, 'color', 'k', 'linewidth', 1.5);
hold on;
patch('XData', [(1:length(preMean))/(1/dsRate), fliplr((1:length(preMean))/(1/dsRate))],...
    'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
    'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
patch('XData', [(1:length(preMean))/(1/dsRate), fliplr((1:length(preMean))/(1/dsRate))],...
    'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);

tempTrl = cell2mat(trialDynFit_Trial);
% tempTrl = cell2mat(trialDynFit_Trial(1:3));
% tempTrl = cell2mat(trialDynFit_Trial(2:end));
trlMean = mean(tempTrl,1,'omitnan');
trlSEM = mlb.SEMcalc(tempTrl,0,1);
trlCI = tinv(0.975, size(tempTrl,1)-1).*trlSEM;
plot((1:length(trlMean))/(1/dsRate),trlMean, 'color', 'r', 'linewidth', 1.5);
hold on;
patch('XData', [(1:length(trlMean))/(1/dsRate), fliplr((1:length(trlMean))/(1/dsRate))],...
    'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
    'linestyle', 'none', 'edgecolor', 'r', 'facecolor', 'r', 'facealpha', 0.25);
patch('XData', [(1:length(trlMean))/(1/dsRate), fliplr((1:length(trlMean))/(1/dsRate))],...
    'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'r', 'facealpha', 0);
xlabel('Persistence duration (ms)');
ylabel('Distance (cosine similarity)');

subplot(1,4,sub2ind([4, 1], 3,1));
tempLats = cell2mat(trialDynFit_Latencies);
% tempLats = cell2mat(trialDynFit_Latencies(1:3));
% tempLats = cell2mat(trialDynFit_Latencies(2:end));
bar([mean(tempLats(:,1)), mean(tempLats(:,2))]);
hold on;
errorbar([mean(tempLats(:,1)), mean(tempLats(:,2))], [mlb.SEMcalc(tempLats(:,1)), mlb.SEMcalc(tempLats(:,2))],...
    'linestyle', 'none', 'color', 'k', 'capsize', 0);
swarmchart(ones(size(tempLats,1),1), tempLats(:,1), 'ok', 'filled', 'markerfacecolor', 'k', 'markeredgecolor','none', 'markerfacealpha', 0.2);
swarmchart(ones(size(tempLats,1),1)+1, tempLats(:,2), 'ok', 'filled', 'markerfacecolor', 'r', 'markeredgecolor','none', 'markerfacealpha', 0.2);
[h,p,ci,stats] = ttest(tempLats(:,1), tempLats(:,2));
if p<0.05
    title(sprintf('t = %.02f; p = %.02i', stats.tstat, p), 'color','r');
else
    title(sprintf('t = %.02f; p = %.02f', stats.tstat, p));
end
set(gca, 'xtick', [1 2], 'xticklabel', [{'Pre-Trial'}, {'Trial'}]);


subplot(1,4,sub2ind([4, 1], 4,1));
corrScatPlot(tempLats(:,1), tempLats(:,2), 'Pre-Trial', 'Trial', []);
%% Plot things
pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);

piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
figure;
colormap(cMap);
for a = 1:length(fileDirs)
    for p = 1:mlb.seqLength
        subplot(length(fileDirs),mlb.seqLength,sub2ind([mlb.seqLength, length(fileDirs)], p,a));
        imagesc(grpAniD{p,p}(:,:,a)');     
        colorbar;    
        set(gca,'ydir', 'normal', 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
        hold on;
        plot(get(gca, 'xlim'),repmat(piNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(poNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(rwdNdx(1), [1,2]), ':k','linewidth', 2);

        plot(repmat(piNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(poNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(rwdNdx(1), [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
        if p==1
            title(fileDirs{a}, 'interpreter', 'none', 'horizontalalignment', 'left');
            ylabel('Train Time');
        end
        if a == length(fileDirs)
            xlabel('Test Time');
        end
            
    end
end

figure;
colormap(cMap);
for p = 1:mlb.seqLength
    % Trial Averages
    subplot(2,mlb.seqLength,sub2ind([mlb.seqLength, 2],p,1));
    tempGroupD = mean(grpTrlD{p,p},3);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempGroupD');    colorbar;    set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
    hold on;
    tempDthresh = chanceTrlD{p,p,1}+(tinv(0.99,numChancePerms-1).*(chanceTrlD{p,p,2}./sqrt(numChancePerms-1)));
    abvThresh = tempGroupD-(tinv(0.99,sum(~isnan(grpTrlD{p,p}(1,1,:)),3)-1).*mlb.SEMcalc(grpTrlD{p,p},0,3))>tempDthresh;
    bounds = bwboundaries(abvThresh);
    for b = 1:length(bounds)
        if numel(bounds{b})>4
            plot(mlb.obsvTimeVect(bounds{b}(:,1)), mlb.obsvTimeVect(bounds{b}(:,2)), 'k', 'linewidth', 2);
        end
    end
    colorbar;
    set(gca,'ydir', 'normal');
    plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
    
    plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    title(sprintf('Trial Average %i',p));
    xlabel('Test Time');
    if p==1
        ylabel('Train Time');
    end
    
    % Animal Averages
    subplot(2,mlb.seqLength,sub2ind([mlb.seqLength, 2],p,2));
    tempGroupD = mean(grpAniD{p,p},3);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempGroupD'); colorbar;    set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
    hold on;
    tempDthresh = chanceSsnD{p,p,1}+(tinv(0.95,numChancePerms-1).*(chanceSsnD{p,p,2}./sqrt(numChancePerms-1)));
%     abvThresh = tempGroupD>tempDthresh;
    abvThresh = tempGroupD-(tinv(0.95,sum(~isnan(grpAniD{p,p}(1,1,:)),3)-1).*mlb.SEMcalc(grpAniD{p,p},0,3))>tempDthresh;
    bounds = bwboundaries(abvThresh);
    for b = 1:length(bounds)
        if numel(bounds{b})>4
            plot(mlb.obsvTimeVect(bounds{b}(:,1)),mlb.obsvTimeVect(bounds{b}(:,2)), 'k', 'linewidth', 2);
        end
    end
    colorbar;
    set(gca,'ydir', 'normal');
    plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
    
    plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    title(sprintf('Animal Average %i',p));
    xlabel('Test Time');
    if p==1
        ylabel('Train Time');
    end
end
%% Plot TransMat
figure;
colormap(cMap);
for p = 1:mlb.seqLength
    for o = 1:mlb.seqLength
        % Trial Averages
        subplot(mlb.seqLength,mlb.seqLength,sub2ind([mlb.seqLength, mlb.seqLength],p,o));
        tempGroupD = mean(grpTrlD{p,o},3);
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempGroupD');    
        colorbar;
        if o==p
            set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
            tempClim = get(gca,'clim');
        end
        hold on;
        tempDthresh = chanceTrlD{p,o,1}+(tinv(0.99,numChancePerms-1).*(chanceTrlD{p,o,2}./sqrt(numChancePerms-1)));
        abvThresh = tempGroupD-(tinv(0.99,sum(~isnan(grpTrlD{p,o}(1,1,:)),3)-1).*mlb.SEMcalc(grpTrlD{p,o},0,3))>tempDthresh;
        bounds = bwboundaries(abvThresh);
        for b = 1:length(bounds)
            if numel(bounds{b})>4
                plot(mlb.obsvTimeVect(bounds{b}(:,1)), mlb.obsvTimeVect(bounds{b}(:,2)), 'k', 'linewidth', 2);
            end
        end
        colorbar;
        set(gca,'ydir', 'normal');
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestPOtime, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestRWDtime, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
        title(sprintf('Pos %i; Decode %i',p,o));
        if o==mlb.seqLength
            xlabel('Test Time');
        end
        if p==1
            ylabel('Train Time');
        end
    end
    for o = 1:mlb.seqLength
        % Trial Averages
        subplot(mlb.seqLength,mlb.seqLength,sub2ind([mlb.seqLength, mlb.seqLength],p,o));
        set(gca,'clim', tempClim);
    end
end

%% Plot Lag Averages
lagTrlVect = cell(size(mlb.lagVect));
for pos = 1:mlb.seqLength
    for odr = 1:mlb.seqLength       
%         if (odr~=1 && pos~=1) && (odr~=4 && pos~=4)
%         if odr~=4 && pos~=4
%         if odr~=1 && pos~=1
            tempLagLog = (pos-odr)==mlb.lagVect;
            lagTrlVect{tempLagLog} = cat(3,lagTrlVect{tempLagLog},grpTrlD{pos,odr});
%         end
    end
end
figure;
colormap(cMap);
for lag = 1:length(mlb.lagVect)
    subplot(1,length(mlb.lagVect), lag);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(lagTrlVect{lag},3)', [-1.5 1.5]); 
    hold on;
    tempDthresh = chanceLag{1,lag}+(tinv(0.99,numChancePerms-1).*(chanceLag{2,lag}./sqrt(numChancePerms-1)));
%     tempDthresh = chance23Lag{1,lag}+(tinv(0.99,numChancePerms-1).*(chance23Lag{2,lag}./sqrt(numChancePerms-1)));
%     tempDthresh = chance123Lag{1,lag}+(tinv(0.99,numChancePerms-1).*(chance123Lag{2,lag}./sqrt(numChancePerms-1)));
%     tempDthresh = chance234Lag{1,lag}+(tinv(0.99,numChancePerms-1).*(chance234Lag{2,lag}./sqrt(numChancePerms-1)));
    abvThresh = mean(lagTrlVect{lag},3)-(tinv(0.99,sum(~isnan(lagTrlVect{lag}),3)-1).*mlb.SEMcalc(lagTrlVect{lag},0,3))>tempDthresh;
    bounds = bwboundaries(abvThresh);
    for b = 1:length(bounds)
        if numel(bounds{b})>4
            plot(mlb.obsvTimeVect(bounds{b}(:,1)), mlb.obsvTimeVect(bounds{b}(:,2)), 'k', 'linewidth', 2);
        end
    end
    colorbar;
    set(gca,'ydir', 'normal');
    plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
    
    plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    title(mlb.lagVect(lag));
    xlabel('Test Time');
    ylabel('Train Time');
end

%% Quantify Off Diagonal Quadrents in decodings
% Create Logical Vectors
if strcmp(alignment{1}, 'PokeIn')
    preTrlPrdLog = mlb.obsvTimeVect<0;
    trlPrdLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<nearestPOtime;
    pstTrlPrdLog = mlb.obsvTimeVect>=nearestPOtime;
elseif strcmp(alignment{1}, 'PokeOut')
    preTrlPrdLog = mlb.obsvTimeVect<nearestPOtime;
    trlPrdLog = mlb.obsvTimeVect<0 & mlb.obsvTimeVect>=nearestPOtime;
    pstTrlPrdLog = mlb.obsvTimeVect>=0;
end
% Pre-Trial & Trial
ttPrTrMtxLog = false(length(mlb.obsvTimeVect));
ttTrPrMtxLog = false(length(mlb.obsvTimeVect));
ttTrPrMtxLog(preTrlPrdLog,trlPrdLog) = true;
ttPrTrMtxLog(trlPrdLog, preTrlPrdLog) = true;
% Pre-Trial & Post-Trial
ttPrPoMtxLog = false(length(mlb.obsvTimeVect));
ttPoPrMtxLog = false(length(mlb.obsvTimeVect));
ttPoPrMtxLog(preTrlPrdLog, pstTrlPrdLog) = true;
ttPrPoMtxLog(pstTrlPrdLog, preTrlPrdLog) = true;
% Trial & Post-Trial
ttPoTrMtxLog = false(length(mlb.obsvTimeVect));
ttTrPoMtxLog = false(length(mlb.obsvTimeVect));
ttTrPoMtxLog(pstTrlPrdLog, trlPrdLog) = true;
ttPoTrMtxLog(trlPrdLog, pstTrlPrdLog) = true;

prTrSVM = cell(size(grpTrlPost));
trPrSVM = cell(size(grpTrlPost));
prPoSVM = cell(size(grpTrlPost));
poPrSVM = cell(size(grpTrlPost));
poTrSVM = cell(size(grpTrlPost));
trPoSVM = cell(size(grpTrlPost));

prTrTrPrSVM = cell(size(grpTrlPost));
trPoPoTrSVM = cell(size(grpTrlPost));
poPrPrPoSVM = cell(size(grpTrlPost));
for op = 1:numel(grpTrlPost)
%     tempPosts = grpTrlPost{op};
    tempPosts = grpTrlD{op};
    tempPrTrSVM = nan(size(tempPosts,3),1);
    tempTrPrSVM = nan(size(tempPosts,3),1);
    tempPrPoSVM = nan(size(tempPosts,3),1);
    tempPoPrSVM = nan(size(tempPosts,3),1);
    tempPoTrSVM = nan(size(tempPosts,3),1);
    tempTrPoSVM = nan(size(tempPosts,3),1);
    tempPrTrTrPrSVM = nan(size(tempPosts,3),1);
    tempTrPoPoTrSVM = nan(size(tempPosts,3),1);
    tempPoPrPrPoSVM = nan(size(tempPosts,3),1);
    for trl = 1:size(tempPosts,3)
        tempPostTrl = tempPosts(:,:,trl);
        tempPrTrSVM(trl) = mean(tempPostTrl(ttPrTrMtxLog));
        tempTrPrSVM(trl) = mean(tempPostTrl(ttTrPrMtxLog));
        tempPrPoSVM(trl) = mean(tempPostTrl(ttPrPoMtxLog));
        tempPoPrSVM(trl) = mean(tempPostTrl(ttPoPrMtxLog));
        tempPoTrSVM(trl) = mean(tempPostTrl(ttPoTrMtxLog));
        tempTrPoSVM(trl) = mean(tempPostTrl(ttTrPoMtxLog));
        
        tempPrTrTrPrSVM(trl) = mean(tempPostTrl(ttPrTrMtxLog | ttTrPrMtxLog));
        tempTrPoPoTrSVM(trl) = mean(tempPostTrl(ttPrPoMtxLog | ttPoPrMtxLog));
        tempPoPrPrPoSVM(trl) = mean(tempPostTrl(ttPoTrMtxLog | ttTrPoMtxLog));
    end
    prTrSVM{op} = tempPrTrSVM;
    trPrSVM{op} = tempTrPrSVM;
    prPoSVM{op} = tempPrPoSVM;
    poPrSVM{op} = tempPoPrSVM;
    poTrSVM{op} = tempPoTrSVM;
    trPoSVM{op} = tempTrPoSVM;
    
    prTrTrPrSVM{op} = tempPrTrTrPrSVM;
    trPoPoTrSVM{op} = tempTrPoPoTrSVM;
    poPrPrPoSVM{op} = tempPoPrPrPoSVM;
end
figure;
% Pre-Trial vs Trial 
subplot(3,6,1)
imagesc(ttPrTrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,2)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(prTrSVM{pos,pos}))+pos-1, prTrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Trial vs Pre-Trial
subplot(3,6,3)
imagesc(ttTrPrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,4)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(trPrSVM{pos,pos}))+pos-1, trPrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Both Pre-Trl & Trial
subplot(3,6,5)
imagesc(ttPrTrMtxLog | ttTrPrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,6)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(prTrTrPrSVM{pos,pos}))+pos-1, prTrTrPrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Pre-Trial vs Post-Trial
subplot(3,6,7)
imagesc(ttPrPoMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,8)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(prPoSVM{pos,pos}))+pos-1, prPoSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Post-Trial vs Pre-Trial
subplot(3,6,9)
imagesc(ttPoPrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,10)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(poPrSVM{pos,pos}))+pos-1, poPrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Both Post-Trial & Pre-Trial
subplot(3,6,11)
imagesc(ttPrPoMtxLog | ttPoPrMtxLog);
set(gca,'ydir', 'normal');
subplot(3,6,12)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(trPoPoTrSVM{pos,pos}))+pos-1, trPoPoTrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Trial vs Post-Trial
subplot(3,6,13)
imagesc(ttTrPoMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,14)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(trPoSVM{pos,pos}))+pos-1, trPoSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Post-Trial vs Trial
subplot(3,6,15)
imagesc(ttPoTrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,16)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(poTrSVM{pos,pos}))+pos-1, poTrSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;
% Both Post-Trial & Trial
subplot(3,6,17)
imagesc(ttTrPoMtxLog | ttPoTrMtxLog);
set(gca, 'ydir', 'normal');
subplot(3,6,18)
for pos = 1:mlb.seqLength
    swarmchart(ones(size(poPrPrPoSVM{pos,pos}))+pos-1, poPrPrPoSVM{pos,pos}, 10, mlb.PositionColors(pos,:));
    hold on;
end
% set(gca, 'ylim', [0 1]);
% set(gca, 'ylim', [-10 10]); grid on;

% Set up big ANOVA
anovaData = [prTrSVM(logical(eye(size(prTrSVM))));...
    trPrSVM(logical(eye(size(trPrSVM))));...
    prPoSVM(logical(eye(size(prPoSVM))));...
    poPrSVM(logical(eye(size(poPrSVM))));...
    trPoSVM(logical(eye(size(trPoSVM))));...
    poTrSVM(logical(eye(size(poTrSVM))))];

windowLog = [cellfun(@(a){ones(size(a))},prTrSVM(logical(eye(size(prTrSVM)))));...
    cellfun(@(a){ones(size(a))+1},trPrSVM(logical(eye(size(trPrSVM)))));...
    cellfun(@(a){ones(size(a))+2},prPoSVM(logical(eye(size(prPoSVM)))));...
    cellfun(@(a){ones(size(a))+3},poPrSVM(logical(eye(size(poPrSVM)))));...
    cellfun(@(a){ones(size(a))+5},trPoSVM(logical(eye(size(trPoSVM)))));...
    cellfun(@(a){ones(size(a))+4},poTrSVM(logical(eye(size(poTrSVM)))))];

symLog = [cellfun(@(a){ones(size(a))},prTrSVM(logical(eye(size(prTrSVM)))));...
    cellfun(@(a){ones(size(a))},trPrSVM(logical(eye(size(trPrSVM)))));...
    cellfun(@(a){ones(size(a))+1},prPoSVM(logical(eye(size(prPoSVM)))));...
    cellfun(@(a){ones(size(a))+1},poPrSVM(logical(eye(size(poPrSVM)))));...
    cellfun(@(a){ones(size(a))+2},trPoSVM(logical(eye(size(trPoSVM)))));...
    cellfun(@(a){ones(size(a))+2},poTrSVM(logical(eye(size(poTrSVM)))))];

posLog = [cellfun(@(a,b){zeros(size(a))+b},prTrSVM(logical(eye(size(prTrSVM)))), num2cell(1:mlb.seqLength)');...
    cellfun(@(a,b){zeros(size(a))+b},trPrSVM(logical(eye(size(trPrSVM)))), num2cell(1:mlb.seqLength)');...
    cellfun(@(a,b){zeros(size(a))+b},prPoSVM(logical(eye(size(prPoSVM)))), num2cell(1:mlb.seqLength)');...
    cellfun(@(a,b){zeros(size(a))+b},poPrSVM(logical(eye(size(poPrSVM)))), num2cell(1:mlb.seqLength)');...
    cellfun(@(a,b){zeros(size(a))+b},trPoSVM(logical(eye(size(trPoSVM)))), num2cell(1:mlb.seqLength)');...
    cellfun(@(a,b){zeros(size(a))+b},poTrSVM(logical(eye(size(poTrSVM)))), num2cell(1:mlb.seqLength)')];

% [p,tbl,stats] = anovan(cell2mat(anovaData')', [cell2mat(windowLog')', cell2mat(symLog')', cell2mat(posLog')'],...
%     'model', 'full', 'sstype', 2,'varnames', [{'Window'}, {'WindowGroup'}, {'Position'}]);

[p,tbl,stats] = anova1(cell2mat(anovaData), cell2mat(windowLog));
multcompare(stats, 'CType', 'bonferroni');

[p,tbl,stats] = anova1(cell2mat([prTrTrPrSVM(logical(eye(size(prTrTrPrSVM)))); trPoPoTrSVM(logical(eye(size(trPoPoTrSVM)))); poPrPrPoSVM(logical(eye(size(poPrPrPoSVM))))]),...
    cell2mat([cellfun(@(a){ones(size(a))},prTrTrPrSVM(logical(eye(size(prTrTrPrSVM))))); cellfun(@(a){ones(size(a))+1},trPoPoTrSVM(logical(eye(size(trPoPoTrSVM))))); cellfun(@(a){ones(size(a))+2},poPrPrPoSVM(logical(eye(size(poPrPrPoSVM)))))]));
multcompare(stats, 'CType', 'bonferroni');



%% Save output
save('PFC_XTD_PosChance_-2to1p5_PO.mat', '-v7.3');
