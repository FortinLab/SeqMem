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
% trlWindow = {[-1000 2000]};
% alignment = {'PokeIn'};
trlWindow = {[-1500 2000]};
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
        grpLFP = cell(mlb.seqLength,length(fileDirs),2);
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
        grpLFP{pos,ani,1} = permute(betaPower(:,:,trlIDs(trlPosVect==pos)), [1,3,2]);
        grpLFP{pos,ani,2} = permute(thetaPower(:,:,trlIDs(trlPosVect==pos)), [1,3,2]);
    end    
end

%% Collapse trial data across animals
grpTrlPost = cell(mlb.seqLength);
grpTrlD = cell(mlb.seqLength);
grpAniD = cell(mlb.seqLength);
grpTrlLFP = cell(mlb.seqLength,2);
for pos = 1:mlb.seqLength
    for p = 1:mlb.seqLength
        grpTrlPost{pos,p} = cell2mat(realPost(pos,p,:));
        grpTrlD{pos,p} = cell2mat(realTrlD(pos,p,:));
        grpAniD{pos,p} = cell2mat(realSsnD(pos,p,:));
    end
    grpTrlLFP{pos,1} = cell2mat(grpLFP(pos,:,1));
    grpTrlLFP{pos,2} = cell2mat(grpLFP(pos,:,2));
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
for perm = 53:numChancePerms
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
    for pos = 1:mlb.seqLength
        for p = 1:mlb.seqLength
            chancePostPerms{pos,p,perm} = mean(cell2mat(tempChancePost(pos,p,:)),3);
            chanceSsnDPerms{pos,p,perm} = mean(cell2mat(tempChanceSsnD(pos,p,:)),3);
            chanceTrlDPerms{pos,p,perm} = mean(cell2mat(tempChanceTrlD(pos,p,:)),3);
        end
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



%% Calculate optimal dynamic model fit
% Create pre-trial dynamic models
preTrlPrdLog = mlb.obsvTimeVect<0;
preTrlLogMtx = false(length(mlb.obsvTimeVect));
preTrlLogMtx(preTrlPrdLog, preTrlPrdLog) = true;
preTrlDynCap = sum(preTrlPrdLog)-1;
preTrlDynMod = nan(sum(preTrlPrdLog), sum(preTrlPrdLog), preTrlDynCap-1);
preTrlDynMod(:,:,1) = eye(sum(preTrlPrdLog));
for d = 1:preTrlDynCap-1
    preTrlDynMod(:,:,d+1) = triu(true(sum(preTrlPrdLog)), d*-1) & tril(true(sum(preTrlPrdLog)), d);
end
% Create trial dynamic models
trlPrdLog = mlb.obsvTimeVect>0 & mlb.obsvTimeVect<nearestPOtime;
trlLogMtx = false(length(mlb.obsvTimeVect));
trlLogMtx(trlPrdLog,trlPrdLog) = true;
trlDynCap = sum(trlPrdLog)-1;
trlDynMod = nan(sum(trlPrdLog), sum(trlPrdLog), trlDynCap-1);
trlDynMod(:,:,1) = eye(sum(trlPrdLog));
for d = 1:trlDynCap-1
    trlDynMod(:,:,d+1) = triu(true(sum(trlPrdLog)), d*-1) & tril(true(sum(trlPrdLog)), d);
end
% Create Early Trial Dynamic Models
rlyTrlPrdLog = mlb.obsvTimeVect>0 & mlb.obsvTimeVect<(nearestPOtime/2);
rlyTrlLogMtx = false(length(mlb.obsvTimeVect));
rlyTrlLogMtx(rlyTrlPrdLog,rlyTrlPrdLog) = true;
rlyTrlDynCap = sum(rlyTrlPrdLog)-1;
rlyTrlDynMod = nan(sum(rlyTrlPrdLog), sum(rlyTrlPrdLog), rlyTrlDynCap-1);
rlyTrlDynMod(:,:,1) = eye(sum(rlyTrlPrdLog));
for d = 1:rlyTrlDynCap-1
    rlyTrlDynMod(:,:,d+1) = triu(true(sum(rlyTrlPrdLog)), d*-1) & tril(true(sum(rlyTrlPrdLog)), d);
end
% Create Late Trial Dynamic Models
latTrlPrdLog = mlb.obsvTimeVect>(nearestPOtime/2) & mlb.obsvTimeVect<nearestPOtime;
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
    trialDynFit_Latencies{pos}(:,1) = preLat./preTrlDynCap;
    [~,trlLat] = min(tempTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,2) = trlLat./trlDynCap;
    [~,rlyLat] = min(tempRlyTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,3) = rlyLat./rlyTrlDynCap;
    [~,latLat] = min(tempLatTrlFit,[],2);
    trialDynFit_Latencies{pos}(:,4) = latLat./latTrlDynCap;
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
    plot((1:length(preMean))./length(preMean),preMean, 'color', 'k', 'linewidth', 1.5);
    hold on;
%     temp = plot(repmat(((1:length(preMean))./length(preMean))', [1,size(tempPreTrl,1)]),tempPreTrl');
%     for t = 1:length(temp)
%         temp(t).Color = [0,0,0,0.2];
%     end
    patch('XData', [(1:length(preMean))./length(preMean), fliplr((1:length(preMean))./length(preMean))],...
        'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
        'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
    patch('XData', [(1:length(preMean))./length(preMean), fliplr((1:length(preMean))./length(preMean))],...
        'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
%     plot(1:length(preMean),preMean, 'color', 'k', 'linewidth', 1.5);
%     hold on;
%     patch('XData', [1:length(preMean), length(preMean):-1:1],...
%         'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
%         'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
%     patch('XData', [1:length(preMean), length(preMean):-1:1],...
%         'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
%         'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
    
    tempTrl = trialDynFit_Trial{pos};
    trlMean = mean(tempTrl,1,'omitnan');
    trlSEM = mlb.SEMcalc(tempTrl,0,1);
    trlCI = tinv(0.975, size(tempTrl,1)-1).*trlSEM;
    plot((1:length(trlMean))./length(trlMean),trlMean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
    hold on;
%     temp = plot(repmat(((1:length(trlMean))./length(trlMean))', [1,size(tempTrl,1)]),tempTrl');
%     for t = 1:length(temp)
%         temp(t).Color = [mlb.PositionColors(pos,:),0.2];
%     end
    patch('XData', [(1:length(trlMean))./length(trlMean), fliplr((1:length(trlMean))./length(trlMean))],...
        'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [(1:length(trlMean))./length(trlMean), fliplr((1:length(trlMean))./length(trlMean))],...
        'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
%     plot(1:length(trlMean),trlMean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
%     hold on;    
%     patch('XData', [1:length(trlMean), length(trlMean):-1:1],...
%         'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
%         'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
%     patch('XData', [1:length(trlMean), length(trlMean):-1:1],...
%         'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
%         'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    xlabel('Proportion Persistent');
    ylabel('Distance (cosine similarity)');
    
    subplot(mlb.seqLength,4,sub2ind([4, mlb.seqLength], 3,pos));
    tempLats = trialDynFit_Latencies{pos};%./dsRate;
    bar([mean(tempLats(:,1)), mean(tempLats(:,2))]);
    hold on;
    scatter(normrnd(0,0.1, [size(tempLats,1),1])+1,tempLats(:,1), 'ok', 'filled', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);
    scatter(normrnd(0,0.1, [size(tempLats,1),1])+2, tempLats(:,2), 'ok', 'filled', 'markerfacecolor', mlb.PositionColors(pos,:), 'markeredgecolor', 'k', 'markerfacealpha', 0.2);
    [h,p,ci,stats] = ttest(tempLats(:,1), tempLats(:,2));
%     annotation(gcf, 'textbox', get(gca, 'position'), 'string', sprintf('t = %.02f; p = %.02i', stats.tstat, p), 'linestyle', 'none');
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
    preMean = mean(tempPreTrl,1,'omitnan');
    preSEM = mlb.SEMcalc(tempPreTrl,0,1);
    preCI = tinv(0.975, size(tempPreTrl,1)-1).*preSEM;
    plot((1:length(preMean))./length(preMean),preMean, 'color', 'k', 'linewidth', 1.5);
    hold on;
%     temp = plot(repmat(((1:length(preMean))./length(preMean))', [1,size(tempPreTrl,1)]),tempPreTrl');
%     for t = 1:length(temp)
%         temp(t).Color = [0,0,0,0.2];
%     end
    patch('XData', [(1:length(preMean))./length(preMean), fliplr((1:length(preMean))./length(preMean))],...
        'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
        'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
    patch('XData', [(1:length(preMean))./length(preMean), fliplr((1:length(preMean))./length(preMean))],...
        'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
%     plot(1:length(preMean),preMean, 'color', 'k', 'linewidth', 1.5);
%     hold on;
%     patch('XData', [1:length(preMean), length(preMean):-1:1],...
%         'YData', [(preMean+preSEM), fliplr(preMean-preSEM)],...
%         'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
%     patch('XData', [1:length(preMean), length(preMean):-1:1],...
%         'YData', [(preMean+preCI), fliplr(preMean-preCI)],...
%         'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
    
    tempTrl = cell2mat(trialDynFit_Trial);
    trlMean = mean(tempTrl,1,'omitnan');
    trlSEM = mlb.SEMcalc(tempTrl,0,1);
    trlCI = tinv(0.975, size(tempTrl,1)-1).*trlSEM;
    plot((1:length(trlMean))./length(trlMean),trlMean, 'color', 'r', 'linewidth', 1.5);
    hold on;
%     temp = plot(repmat(((1:length(trlMean))./length(trlMean))', [1,size(tempTrl,1)]),tempTrl');
%     for t = 1:length(temp)
%         temp(t).Color = [mlb.PositionColors(pos,:),0.2];
%     end
    patch('XData', [(1:length(trlMean))./length(trlMean), fliplr((1:length(trlMean))./length(trlMean))],...
        'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
        'linestyle', 'none', 'edgecolor', 'r', 'facecolor', 'r', 'facealpha', 0.25);
    patch('XData', [(1:length(trlMean))./length(trlMean), fliplr((1:length(trlMean))./length(trlMean))],...
        'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'r', 'facealpha', 0);
%     plot(1:length(trlMean),trlMean, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
%     hold on;    
%     patch('XData', [1:length(trlMean), length(trlMean):-1:1],...
%         'YData', [(trlMean+trlSEM), fliplr(trlMean-trlSEM)],...
%         'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
%     patch('XData', [1:length(trlMean), length(trlMean):-1:1],...
%         'YData', [(trlMean+trlCI), fliplr(trlMean-trlCI)],...
%         'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    xlabel('Proportion Persistent');
    ylabel('Distance (cosine similarity)');
    
    subplot(1,4,sub2ind([4, 1], 3,1));
    tempLats = cell2mat(trialDynFit_Latencies);%./dsRate;
    bar([mean(tempLats(:,1)), mean(tempLats(:,2))]);
    hold on;
    scatter(normrnd(0,0.1, [size(tempLats,1),1])+1,tempLats(:,1), 'ok', 'filled', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);
    scatter(normrnd(0,0.1, [size(tempLats,1),1])+2, tempLats(:,2), 'ok', 'filled', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);
    [h,p,ci,stats] = ttest(tempLats(:,1), tempLats(:,2));
%     annotation(gcf, 'textbox', get(gca, 'position'), 'string', sprintf('t = %.02f; p = %.02i', stats.tstat, p), 'linestyle', 'none');
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
        imagesc(grpAniD{p,p}(:,:,a)', [-2 2]);
%         imagesc(zscore(grpAniD{p,p}(:,:,a)',0, 'all'), [-2 2]);
        colorbar;
        set(gca,'ydir', 'normal');
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
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, zscore(tempGroupD',0,'all'), [-2 2]);    
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
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, zscore(tempGroupD',0,'all'), [-2 2]);    
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
%%
figure;
colormap(cMap);
for p = 1:mlb.seqLength
    for o = 1:mlb.seqLength
        % Trial Averages
        subplot(mlb.seqLength,mlb.seqLength,sub2ind([mlb.seqLength, mlb.seqLength],p,o));
        tempGroupD = mean(grpTrlD{p,o},3);
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, tempGroupD');    colorbar;    %set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
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
    for odr = 2:3
%         if odr~=4 && pos~=4
            tempLagLog = (pos-odr)==mlb.lagVect;
            lagTrlVect{tempLagLog} = cat(3,lagTrlVect{tempLagLog},grpTrlD{pos,odr});
%         end
    end
end
figure;
colormap(cMap);
for lag = 1:length(mlb.lagVect)
    subplot(1,length(mlb.lagVect), lag);
    imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, zscore(mean(lagTrlVect{lag},3),0,'all')', [-2 2]); 
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(lagTrlVect{lag},3)', [-1 1]); 
%     imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(lagTrlVect{lag},3)'); 
%     if lag==0
%         set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
%         tempClim = get(gca,'clim');
%     end
    colorbar;
    set(gca,'ydir', 'normal');
    hold on;
    plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
    
    plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestPOtime, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(nearestRWDtime, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    title(mlb.lagVect(lag));
end
% 
% for lag = 1:length(mlb.lagVect)
%     subplot(1,length(mlb.lagVect), lag);
%     set(gca, 'clim', tempClim);
% end
        
    

%% Quiiiiick look at beta and xtd
% figure;
% trl = 2;
% for t = 1:size(grpTrlLFP{1,1},2)
%     
%     sp1 = subplot(3,1,1:2);
%     imagesc(grpTrlD{trl,trl}(:,:,t)',[-2 2]);
%     hold on;
%     plot(get(gca, 'xlim'),repmat(piNdx(1), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(poNdx(1), [1,2]), '--k','linewidth', 2);
%     plot(get(gca, 'xlim'),repmat(rwdNdx(1), [1,2]), ':k','linewidth', 2);
%     
%     plot(repmat(piNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(1), [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     hold off;
%     set(gca, 'ydir', 'normal');
%     sp2 = subplot(3,1,3);
%     plot(grpTrlLFP{trl,1}(:,t), 'k');
%     hold on;
%     plot(grpTrlLFP{trl,2}(:,t), 'r');
%     plot(repmat(piNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(1), [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     hold off;
%     axis tight
%     linkaxes([sp1, sp2], 'x');
%     
%     pause(2);
% end
% 
% %% Quiiiiick look at beta and xtd (XTDd diagonal across positions)
% figure;
% trl = 2;
% for t = 1:size(grpTrlLFP{1,1},2)
%     
%     sp1 = subplot(2,1,1);
%     tempDecode = cell2mat(cellfun(@(a){diag(a)},cellfun(@(a){a(:,:,t)},grpTrlD(trl,:))));
%     for pos = 1:mlb.seqLength
%         plot(tempDecode(:,pos),'color',mlb.PositionColors(pos,:), 'linewidth', 2);
%         hold on;
%     end    
%     axis tight;
%     set(gca,'ylim', [0 10]);
%     plot(repmat(piNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(1), [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     hold off;
%     set(gca, 'ydir', 'normal');
%     sp2 = subplot(2,1,2);
%     plot(grpTrlLFP{trl,1}(:,t), 'k');
%     hold on;
%     plot(grpTrlLFP{trl,2}(:,t), 'r');
%     plot(repmat(piNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(poNdx(1), [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%     plot(repmat(rwdNdx(1), [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     hold off;
%     axis tight
%     set(gca, 'ylim', [-2 8]);
%     linkaxes([sp1, sp2], 'x');
%     
%     pause(2);
% end

clear chancePostPerms chanceSsnDPerms chanceTrlDPerms
save('PFC_XTD_PosChance_PO.mat', '-v7.3');
