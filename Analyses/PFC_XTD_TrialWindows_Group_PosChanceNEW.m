% PFC_XTD_TrialWindows_Group_PosChance

%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'}
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];
% % CA1 Data
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% tets = [1,22,17,18,17]; % Lateral/Distal
% % tets = [7,3,1,5,5]; % Medial/Proximal

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
binSize = 200;
dsRate = 50;
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts
% 
% alignments = [{'PokeIn'}, {'PokeOut'}];
% trlWindows = [{[-1200 2000]}, {[-2000 1200]}];
% trlWindows = [{[-700 2000]}, {[-2000 700]}]; %% Non-Overlapping Pre/Post ITD

alignments = {'PokeIn'};
% trlWindows = {[-1200 2000]}; %% Overlapping Pre/Post ITD... used for persistence evaluation
trlWindows = {[-700 2000]}; %% Non-Overlapping Pre/Post ITD

% alignments = {'PokeOut'};
% trlWindows = {[-2000 1200]};  %% Overlapping Pre/Post ITD... used for persistence evaluation
% trlWindows = {[-2000 700]}; %% Non-Overlapping Pre/Post ITD
    
lfpWindow = [16 32];
numChancePerms = 100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Chance Perm Nums

cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);
%% Create the mlb objects
realTic = tic;
mlb = cell(size(fileDirs));
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb{ani} = MLB_SM(fileDirs{ani});
end

%% Extract Behavioral Variables
piPokeOutLat = cell(size(fileDirs));
piRwdLat = cell(size(fileDirs));
poPokeInLat = cell(size(fileDirs));
poRwdLat = cell(size(fileDirs));
for ani = 1:length(fileDirs)
    piPokeOutLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex])'/1000;
    piRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex])'/1000;
    poPokeInLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex])'/1000;
    poRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex])'/1000;
end
grpPiPoLat = median(cell2mat(piPokeOutLat'));
grpPiRwdLat = median(cell2mat(piRwdLat'));
grpPoPiLat = median(cell2mat(poPokeInLat'));
grpPoRwdLat = median(cell2mat(poRwdLat'));

%% Run XTD
xtdPostRaw = cell(length(alignments),length(fileDirs));
xtdISClog = cell(1,length(fileDirs));
xtdTrlIDs = cell(1,length(fileDirs));
obsvTimeVect = cell(size(alignments));
for ani = 1:length(fileDirs)
    fprintf('Processing XTD for Anim #%i.....',ani);
    mlb{ani}.binSize = binSize;
    mlb{ani}.dsRate = dsRate;
    mlb{ani}.bayesType = bayesType;    
    seqs = unique([mlb{ani}.trialInfo.SequenceNum]);
    % Make a blank trialInfo struct
    blankStruct = mlb{ani}.trialInfo(1);
    tempFields = fieldnames(blankStruct);
    for f = 1:length(tempFields)
        blankStruct.(tempFields{f}) = nan;
    end
    tempSeqInfo = repmat(blankStruct,[mlb{ani}.seqLength, max(seqs)]);
    for s = 1:length(seqs)
        tempSeqInfo(1:sum([mlb{ani}.trialInfo.SequenceNum]==s),s) = mlb{ani}.trialInfo([mlb{ani}.trialInfo.SequenceNum]==s);
    end
    xtdTrlIDs{ani} = tempSeqInfo;
    tempSeqPosts = cell(numel(mlb{ani}.odrSeqs),max(seqs));
    for al = 1:length(alignments)
        mlb{ani}.alignments = alignments(al);
        mlb{ani}.windows = trlWindows(al);
        mlb{ani}.SetLikes_ISC;
        mlb{ani}.Process_IterativeLikelyL1O;
        if ani == 1
            obsvTimeVect{al} = mlb{ani}.obsvTimeVect;
        end
        tempPost = permute(mlb{ani}.post, [1,3,2]);
        trlIDs = unique(mlb{ani}.postTrlIDs(~isnan(mlb{ani}.postTrlIDs)));
        for trl = 1:length(trlIDs)
            trlInfo = mlb{ani}.trialInfo(trlIDs(trl));
            tempSeqPosts{trlInfo.Position, trlInfo.SequenceNum} = tempPost{mlb{ani}.postTrlIDs==trlIDs(trl)};
        end
        if al==1
            xtdISClog{ani} = cellfun(@(a)~isempty(a),tempSeqPosts);
        end
        % Run XTD on all non-ISC trials
        nonISCtrls = xtdTrlIDs{ani}(~xtdISClog{ani});
        nonISCtrls(arrayfun(@(a)isnan(a.TrialNum),nonISCtrls)) = [];
        [ssnSpikes, ~] = mlb{ani}.PP_ConcatTrialData;
        tempLike = cell2mat(mlb{ani}.likeTrlSpikes);
        tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));
        for nit = 1:length(nonISCtrls)
            if ~isempty(tempSeqPosts{nonISCtrls(nit).Position, nonISCtrls(nit).SequenceNum})
                error('Something wrong this should no have data dummy code monkey!');
            end
            tempSeqPosts{nonISCtrls(nit).Position, nonISCtrls(nit).SequenceNum} = mlb{ani}.CalcIterativeBayesPost_Poisson(mean(tempLike,3,'omitnan'), ssnSpikes(:,:,nonISCtrls(nit).TrialNum), mlb{ani}.decodeIDvects(:,1), mlb{ani}.decodeIDvects(:,4), tempProb);
        end
        unqPostVals = unique(cell2mat(reshape(cellfun(@(a){a(:)},tempSeqPosts),[numel(tempSeqPosts),1])));
        for t = 1:numel(tempSeqPosts)
            if ~isempty(tempSeqPosts{t})
                tempSeqPosts{t}(tempSeqPosts{t}==0) = unqPostVals(find(unqPostVals>0,1,'first'));
                tempSeqPosts{t}(tempSeqPosts{t}==1) = unqPostVals(find(unqPostVals<1,1,'last'));
            end
        end
        xtdPostRaw{al,ani} = tempSeqPosts;
    end
    fprintf('Completed\n');
end

%% Calculate D-Prime values (per-animal)
aniFAR = cell(length(fileDirs), length(alignments));
aniFARtime = cell(length(fileDirs), length(alignments));
aniTrlD = cell(size(xtdPostRaw));       % Decodability with session level estimates of FAR
aniTrlDtime = cell(size(xtdPostRaw));   % Decodability with local level (+/- 5 sequences) estimates of FAR
seqTimeWin = 3;
for ani = 1:length(fileDirs)
    fprintf('Processing Decodability for Anim #%i.....',ani);
    seqNums = unique([mlb{ani}.trialInfo.SequenceNum]);
    for al = 1:length(alignments)
        iscLog = xtdISClog{ani};
        trlIDs = xtdTrlIDs{ani};
        iscTrlIDs = xtdTrlIDs{ani}(iscLog);
        posIDs = [iscTrlIDs.Position];
        iscPosts = permute(xtdPostRaw{al,ani}(iscLog), [2,3,1]);
        tempFAR = cell(mlb{ani}.seqLength);
        tempFARtime = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, max(seqNums));
        for trlPos = 1:mlb{ani}.seqLength
            tpLog = [iscTrlIDs.Position]==trlPos;
            tempFAposts = iscPosts(~tpLog);
            for decPos = 1:mlb{ani}.seqLength
                tempFAR{trlPos,decPos} = mean(cell2mat(cellfun(@(a){a(:,:,decPos)},tempFAposts)),3,'omitnan');
                for seq = 1:max(seqNums)
                    if seq < (seqTimeWin*2)
                        seqLog = [iscTrlIDs.SequenceNum]<=seqTimeWin*2;
                    elseif seq > max(seqNums)-(seqTimeWin*2)
                        seqLog = [iscTrlIDs.SequenceNum]>=max(seqNums)-(seqTimeWin*2);
                    else
                        seqLog = [iscTrlIDs.SequenceNum]<=seq+(seqTimeWin-1) & [iscTrlIDs.SequenceNum]>=seq-seqTimeWin;
                    end
                    tempFARtime{trlPos,decPos,seq} = mean(cell2mat(cellfun(@(a){a(:,:,decPos)},iscPosts(~tpLog & seqLog))),3,'omitnan');
                end
            end     
        end
        aniFAR{ani,al} = tempFAR;
        aniFARtime{ani,al} = tempFARtime;
        
        posts = xtdPostRaw{al,ani};
        trlDs = cell(size(posts));
        trlDsTime = cell(size(posts));
        for pos = 1:size(posts,1)
            curFAR = tempFAR(pos,:);
            curFARtime = tempFARtime(pos,:,:);
            for seq = 1:size(posts,2)
                trlPost = posts{pos,seq};
                if ~isempty(trlPost)
                    tempTrlD = nan(size(trlPost));
                    tempTrlDtime = nan(size(trlPost));
                    for posD = 1:size(trlPost,3)
                        %%%%% Straight Difference HR-FAR
%                         tempTrlD(:,:,posD) = trlPost(:,:,posD)-curFAR{posD};
%                         if ~isempty(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum})
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)-curFARtime{1,posD,trlIDs(pos,seq).SequenceNum};
%                         elseif trlIDs(pos,seq).SequenceNum == max(seqNums)
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)-curFARtime{1,posD,trlIDs(pos,seq).SequenceNum-1};
%                         else
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)-curFARtime{1,posD,trlIDs(pos,seq).SequenceNum+1};
%                         end
%                         
                        %%%%% Ratio Difference HR-FAR
%                         tempTrlD(:,:,posD) = trlPost(:,:,posD)./curFAR{posD};
%                         if ~isempty(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum})
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)./curFARtime{1,posD,trlIDs(pos,seq).SequenceNum};
%                         elseif trlIDs(pos,seq).SequenceNum == max(seqNums)
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)./curFARtime{1,posD,trlIDs(pos,seq).SequenceNum-1};
%                         else
%                             tempTrlDtime(:,:,posD) = trlPost(:,:,posD)./curFARtime{1,posD,trlIDs(pos,seq).SequenceNum+1};
%                         end
                        
                        %%%%% D' Difference HR-FAR
                        tempTrlD(:,:,posD) = norminv(trlPost(:,:,posD))-norminv(curFAR{posD});
                        if ~isempty(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum})
                            tempTrlDtime(:,:,posD) = norminv(trlPost(:,:,posD))-norminv(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum});
                        elseif trlIDs(pos,seq).SequenceNum == max(seqNums)
                            tempTrlDtime(:,:,posD) = norminv(trlPost(:,:,posD))-norminv(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum-1});
                        else
                            tempTrlDtime(:,:,posD) = norminv(trlPost(:,:,posD))-norminv(curFARtime{1,posD,trlIDs(pos,seq).SequenceNum+1});
                        end
                    end
                    trlDs{pos,seq} = tempTrlD;
                    trlDsTime{pos,seq} = tempTrlDtime;
                end
            end
        end
        aniTrlD{al,ani} = trlDs;    
        aniTrlDtime{al,ani} = trlDsTime;
    end
    fprintf('Completed\n');
end
            
%% Trial Specific Analyses: Epoched Values & Model Fits
trialEpochs = cell(size(aniTrlD));
trialEpochComps = cell(size(aniTrlD));
trialEpochFits = cell(size(aniTrlD));
curD = aniTrlD; % <----- Session level FAR d' value
% curD = aniTrlDtime; % <----- Local level FAR d' value
for al = 1:length(alignments)
    for ani = 1:length(fileDirs)
        tempD = curD{al,ani};
        tempTrlInfo = xtdTrlIDs{ani};
        tempTmatEpochs = cell(size(tempD));
        tempTmatEpochComp = cell(size(tempD));
        tempEpochFits = cell(size(tempD));
        for trl = 1:numel(tempD)
            if ~isnan(tempTrlInfo(trl).TrialNum)
                % Identify Trial Periods
                curTrlInfo = tempTrlInfo(trl);
                if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                    nearestPIval = 0;
                    nearestPOval = mlb{ani}.obsvTimeVect(find(curTrlInfo.PokeDuration>mlb{ani}.obsvTimeVect,1,'last'));
                elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                    nearestPIval = mlb{ani}.obsvTimeVect(find(curTrlInfo.PokeDuration>mlb{ani}.obsvTimeVect,1,'last'))*-1;
                    nearestPOval = 0;
                end
                % Define Trial Period Logical Masks
                preTrialLog = obsvTimeVect{al}<nearestPIval;
                trialLog = obsvTimeVect{al}>=nearestPIval & obsvTimeVect{al}<nearestPOval;
                pstTrialLog = obsvTimeVect{al}>=nearestPOval;
                if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                    itiLog = preTrialLog;
                elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                    itiLog = pstTrialLog;
                end
                % Define Epoch Masks
                % NOTE: In the data, resulting from the initial analysis designs rows refer to trial time, columns refer to model time. 
                % For visualizations this is flipped... it doesn't need be, that's just how we've been looking at things. All the diagrams below refer
                % to this real data organization with trial time as rows and training time as columns.
                blankMask = false(size(tempD{1}(:,:,1)));
                %   Train/Test Pre/Pre
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       ----------------
                prPrMsk = blankMask;
                prPrMsk(preTrialLog, preTrialLog) = true;
                %   Train/Test Pre/Trial
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       |%%%%|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                prTrMsk = blankMask;    
                prTrMsk(trialLog, preTrialLog) = true;
                %   Train/Test Pre/Post
                %       ----------------
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       |%%%%|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                prPoMsk = blankMask;    
                prPoMsk(pstTrialLog, preTrialLog) = true;
                %   Train/Test Trial/Pre
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       ----------------
                trPrMsk = blankMask;    
                trPrMsk(preTrialLog, trialLog) = true;
                %   Train/Test Trial/Trial
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |____|%%%%|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                trTrMsk = blankMask;    
                trTrMsk(trialLog, trialLog) = true;
                %   Train/Test Trial/Post
                %       ----------------
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |____|%%%%|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                trPoMsk = blankMask;    
                trPoMsk(pstTrialLog, trialLog) = true;
                %   Train/Test Post/Pre
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |%%%%|
                %       |    |    |%%%%|
                %       |    |    |%%%%|
                %       ----------------
                poPrMsk = blankMask;    
                poPrMsk(preTrialLog, pstTrialLog) = true;
                %   Train/Test Post/Trial
                %       ----------------
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |%%%%|
                %       |    |    |%%%%|
                %       |____|____|%%%%|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                poTrMsk = blankMask;   
                poTrMsk(trialLog, pstTrialLog) = true;
                %   Train/Test Post/Post
                %       ----------------
                %       |    |    |%%%%|
                %       |    |    |%%%%|
                %       |____|____|%%%%|
                %       |    |    |    |
                %       |    |    |    |
                %       |____|____|____|
                %       |    |    |    |
                %       |    |    |    |
                %       |    |    |    |
                %       ----------------
                poPoMsk = blankMask;    
                poPoMsk(pstTrialLog, pstTrialLog) = true;
                %   Within Epoch
                %       ----------------
                %       |    |    |%%%%|
                %       |    |    |%%%%|
                %       |____|____|%%%%|
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |____|%%%%|____|
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       |%%%%|    |    |
                %       ----------------
                weMsk = blankMask;      
                weMsk(prPrMsk | trTrMsk | poPoMsk) = true;
                % Across Epoch
                %       ----------------
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |____|%%%%|____|
                %       |%%%%|    |%%%%|
                %       |%%%%|    |%%%%|
                %       |%%%%|____|%%%%|
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       |    |%%%%|    |
                %       ----------------
                xeMsk = blankMask;          
                xeMsk(prTrMsk | trPrMsk | poTrMsk | trPoMsk) = true;
                % Extract Epoch Values
                curTrlD = tempD{trl};
                tempEpochMeans = nan(3,3,mlb{ani}.seqLength);
                tempEpochComp = nan(2,1,mlb{ani}.seqLength);
                % Trial Dynamics Analysis Parameters
                trialDynCap = sum(trialLog)-1;
                itiDynCap = sum(itiLog)-1;
                tempTrialEpochFit = nan(2,length(obsvTimeVect{al}),mlb{ani}.seqLength);
                for pos = 1:mlb{ani}.seqLength
                    curTrlDpos = curTrlD(:,:,pos);
                    % NOTE: organization for these summary values is transposed relative to the organization outlined above to be more in line with
                    % the visuals. Rows here are training time (Pre = 1; Trial = 2; Post = 3) with columns reflecting testing time.
                    tempEpochMeans(1,1,pos) = mean(curTrlDpos(prPrMsk));    
                    tempEpochMeans(1,2,pos) = mean(curTrlDpos(prTrMsk));
                    tempEpochMeans(1,3,pos) = mean(curTrlDpos(prPoMsk));
                    tempEpochMeans(2,1,pos) = mean(curTrlDpos(trPrMsk));
                    tempEpochMeans(2,2,pos) = mean(curTrlDpos(trTrMsk));
                    tempEpochMeans(2,3,pos) = mean(curTrlDpos(trPoMsk));
                    tempEpochMeans(3,1,pos) = mean(curTrlDpos(poPrMsk));
                    tempEpochMeans(3,2,pos) = mean(curTrlDpos(poTrMsk));
                    tempEpochMeans(3,3,pos) = mean(curTrlDpos(poPoMsk));
                    tempEpochComp(1,1,pos) = mean(curTrlDpos(weMsk));
                    tempEpochComp(2,1,pos) = mean(curTrlDpos(xeMsk));
                    
                    trialPeriod = zscore(curTrlDpos(trialLog,trialLog), 0, 'all');
                    for d = 1:trialDynCap-1
                        tempDyn = zscore(triu(true(sum(trialLog)), d*-1) & tril(true(sum(trialLog)), d), 0, 'all');
                        tempTrialEpochFit(1,d,pos) = pdist([trialPeriod(:)';tempDyn(:)'], 'cosine');
                    end
                    itiPeriod = zscore(curTrlDpos(itiLog,itiLog), 0 , 'all');
                    for d = 1:itiDynCap-1
                        tempDyn = zscore(triu(true(sum(itiLog)), d*-1) & tril(true(sum(itiLog)), d), 0, 'all');
                        tempTrialEpochFit(2,d,pos) = pdist([itiPeriod(:)';tempDyn(:)'], 'cosine');
                    end
                end
                tempTmatEpochs{trl} = tempEpochMeans;
                tempTmatEpochComp{trl} = tempEpochComp;
                tempEpochFits{trl} = tempTrialEpochFit; 
            end
        end
        trialEpochs{al,ani} = tempTmatEpochs;
        trialEpochComps{al,ani} = tempTmatEpochComp;
        trialEpochFits{al,ani} = tempEpochFits;
    end
end
fprintf('Trial analyses complete\n');

%% Combine Trial Decodability Across Animals into Transposition Matrix Organization
%%%% Can rework this via [variable{:}] to skip the initial temporary reshaping step... do this tomorrow.
trlD_tMat = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
trlDtime_tMat = trlD_tMat;
trlEpochs_tMat = trlD_tMat;
trlEpochComp_tMat = trlD_tMat;
trlFits_tMat = trlD_tMat;
trlEpochs_tMat_OSC = trlD_tMat;
trlFits_tMat_OSC = trlD_tMat;
trlEpochs_tMat_ERR = trlD_tMat;
trlFits_tMat_ERR = trlD_tMat;
for al = 1:length(alignments)
    % Concatenate logical/selection matrices
    tempISClog = [xtdISClog{:}];
    tempTrlInfo = [xtdTrlIDs{:}];
    tempOSClog = reshape([tempTrlInfo.TranspositionDistance], size(tempTrlInfo))~=0 & reshape([tempTrlInfo.Performance], size(tempTrlInfo))==1;
    tempErrLog = reshape([tempTrlInfo.Performance], size(tempTrlInfo))==0;
    % Concatenate data matrices
    tempD = [aniTrlD{al,:}];
    tempDtime = [aniTrlDtime{al,:}];
    tempTrlEpochs = [trialEpochs{al,:}];
    tempTrlEpochComps = [trialEpochComps{al,:}];
    tempTrlFits = [trialEpochFits{al,:}];
    % Select Data
    %   ISC
    tempD_ISC = tempD;
    tempD_ISC(~tempISClog) = {[]};
    tempDtime_ISC = tempDtime;
    tempDtime_ISC(~tempISClog) = {[]};
    tempTrlEpochs_ISC = tempTrlEpochs;
    tempTrlEpochs_ISC(~tempISClog) = {[]};
    tempTrlEpochComps_ISC = tempTrlEpochComps;
    tempTrlEpochComps_ISC(~tempISClog) = {[]};
    tempTrlFits_ISC = tempTrlFits;
    tempTrlFits_ISC(~tempISClog) = {[]};
    %   OSC
    tempTrlEpochs_OSC = tempTrlEpochs;
    tempTrlEpochs_OSC(~tempOSClog) = {[]};
    tempTrlFits_OSC = tempTrlFits;
    tempTrlFits_OSC(~tempOSClog) = {[]};
    %   ERR
    tempTrlEpochs_ERR = tempTrlEpochs;
    tempTrlEpochs_ERR(~tempErrLog) = {[]};
    tempTrlFits_ERR = tempTrlFits;
    tempTrlFits_ERR(~tempErrLog) = {[]};
    
    % Organize data into tMat shape
    for posT = 1:mlb{ani}.seqLength
        tempD_ISCpos = permute(tempD_ISC(posT,:), [1,3,2]);
        tempD_ISCpos(cellfun(@(a)isempty(a), tempD_ISCpos)) = [];
        tempDtime_ISCpos = permute(tempDtime_ISC(posT,:), [1,3,2]);
        tempDtime_ISCpos(cellfun(@(a)isempty(a), tempDtime_ISCpos)) = [];
        tempTrlEpochs_ISCpos = permute(tempTrlEpochs_ISC(posT,:), [1,3,2]);
        tempTrlEpochs_ISCpos(cellfun(@(a)isempty(a), tempTrlEpochs_ISCpos)) = [];
        tempTrlEpochComps_ISCpos = permute(tempTrlEpochComps_ISC(posT,:), [1,3,2]);
        tempTrlEpochComps_ISCpos(cellfun(@(a)isempty(a), tempTrlEpochComps_ISCpos)) = [];
        tempTrlFits_ISCpos = permute(tempTrlFits_ISC(posT,:), [1,3,2]);
        tempTrlFits_ISCpos(cellfun(@(a)isempty(a), tempTrlFits_ISCpos)) = [];
        tempTrlEpochs_OSCpos = permute(tempTrlEpochs_OSC(posT,:), [1,3,2]);
        tempTrlEpochs_OSCpos(cellfun(@(a)isempty(a), tempTrlEpochs_OSCpos)) = [];
        tempTrlFits_OSCpos = permute(tempTrlFits_OSC(posT,:), [1,3,2]);
        tempTrlFits_OSCpos(cellfun(@(a)isempty(a), tempTrlFits_OSCpos)) = [];
        tempTrlEpochs_ERRpos = permute(tempTrlEpochs_ERR(posT,:), [1,3,2]);
        tempTrlEpochs_ERRpos(cellfun(@(a)isempty(a), tempTrlEpochs_ERRpos)) = [];
        tempTrlFits_ERRpos = permute(tempTrlFits_ERR(posT,:), [1,3,2]);
        tempTrlFits_ERRpos(cellfun(@(a)isempty(a), tempTrlFits_ERRpos)) = [];
        for posD = 1:mlb{ani}.seqLength            
            trlD_tMat{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempD_ISCpos));
            trlDtime_tMat{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempDtime_ISCpos));
            trlEpochs_tMat{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlEpochs_ISCpos));
            trlEpochComp_tMat{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlEpochComps_ISCpos));
            trlFits_tMat{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlFits_ISCpos));
            trlEpochs_tMat_OSC{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlEpochs_OSCpos));
            trlFits_tMat_OSC{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlFits_OSCpos));
            trlEpochs_tMat_ERR{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlEpochs_ERRpos));
            trlFits_tMat_ERR{posT,posD,al} = cell2mat(cellfun(@(a){a(:,:,posD)}, tempTrlFits_ERRpos));
        end
    end
end
fprintf('Trial data complied\n')

%% Z-Normalize Relevant Data to True Correct
trlEpochs_tMatZ = cell(size(trlEpochs_tMat));
trlEpochComp_tMatZ = cell(size(trlEpochComp_tMat));
trlEpochs_tMat_OSCz = cell(size(trlEpochs_tMat_OSC));
trlEpochs_tMat_ERRz = cell(size(trlEpochs_tMat_ERR));
for pos = 1:mlb{ani}.seqLength
    tempEpochMean = mean(reshape(trlEpochs_tMat{pos,pos}, [numel(trlEpochs_tMat{pos,pos}), 1]));
    tempEpochVar = std(reshape(trlEpochs_tMat{pos,pos}, [numel(trlEpochs_tMat{pos,pos}), 1]));
    trlEpochs_tMatZ(pos,:) = cellfun(@(a){(a-tempEpochMean)./tempEpochVar}, trlEpochs_tMat(pos,:));
    tempEpochCompMean = mean(reshape(trlEpochComp_tMat{pos,pos}, [numel(trlEpochComp_tMat{pos,pos}), 1]));
    tempEpochCompVar = std(reshape(trlEpochComp_tMat{pos,pos}, [numel(trlEpochComp_tMat{pos,pos}), 1]));
    trlEpochComp_tMatZ(pos,:) = cellfun(@(a){(a-tempEpochCompMean)./tempEpochCompVar}, trlEpochComp_tMat(pos,:));
    tempEpochOSCmean = mean(reshape(trlEpochs_tMat_OSC{pos,pos}, [numel(trlEpochs_tMat_OSC{pos,pos}), 1]));
    tempEpochOSCvar = std(reshape(trlEpochs_tMat_OSC{pos,pos}, [numel(trlEpochs_tMat_OSC{pos,pos}), 1]));
%     trlEpochs_tMat_OSCz(pos,:) = cellfun(@(a){(a-tempEpochOSCmean)./tempEpochOSCvar}, tempEpochs_tMat_OSC(pos,:));
end
fprintf('Epoch data normalized\n');

%% Within vs Across Epoch Analyses
% posVect = 1:mlb{ani}.seqLength; tit = 'All';
% posVect = 1:mlb{ani}.seqLength-1; tit = 'All but Final';
% posVect = 2:3;    tit = 'Pos 2&3';
posVect = 2:mlb{ani}.seqLength;   tit = 'SANSA';
% posVect = 3; tit = sprintf('Pos=%i',posVect);

% epochCompData = trlEpochComp_tMat;
epochCompData = trlEpochComp_tMatZ;

binNum = 100;
for al = 1:length(alignments)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Pre-Processing
    curLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
    curLog(:,:,al) = logical(eye(mlb{ani}.seqLength));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Analyses
    curEpochComps = permute(epochCompData(curLog), [2,3,1]);
    curEpochComps = cell2mat(curEpochComps(posVect));

    figure; 
    mlb{ani}.PlotMeanVarViolin(1,curEpochComps(1,1,:),3,0.05,'r','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(2,curEpochComps(2,1,:),3,0.05,'k','error', 'CI', 'binNum', binNum);  
    [~,p,~,stats] = ttest(curEpochComps(1,1,:),curEpochComps(2,1,:));
    if p<0.001
        text(0.25, max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02i', p)}], 'verticalalignment', 'top');
    else
        text(0.25, max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02f', p)}], 'verticalalignment', 'top');
    end
    title('Epoched XTD');
    set(gca,'xtick', 1:2, 'xticklabels', [{'Within'}, {'Across'}]);
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s: %s', alignments{al}, tit),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% Lag Analyses
posVect = 1:mlb{ani}.seqLength; tit = 'All';
% posVect = 1:mlb{ani}.seqLength-1; tit = 'All but Final';
% posVect = 2:3;    tit = 'Pos 2&3';
% posVect = 2:mlb{ani}.seqLength;   tit = 'SANSA';
% posVect = 3; tit = sprintf('Pos=%i',posVect);

epochData = trlEpochs_tMat;
% epochData = trlEpochs_tMatZ;

binNum = 100;

% Logical Vectors for tMat space
for al = 1:length(alignments)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Pre-Processing
    prvLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
    prvLog(:,:,al) = logical(triu(ones(mlb{ani}.seqLength),-1) & tril(ones(mlb{ani}.seqLength),-1));
    curLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
    curLog(:,:,al) = logical(eye(mlb{ani}.seqLength));
    nxtLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
    nxtLog(:,:,al) = logical(triu(ones(mlb{ani}.seqLength),1) & tril(ones(mlb{ani}.seqLength),1));
    % Extract Lag Data
    nxtEpochs = permute(epochData(nxtLog), [2,3,1]);
    prvEpochs = permute(epochData(prvLog), [2,3,1]);
    curEpochs = permute(epochData(curLog), [2,3,1]);
    if strcmp(tit, 'SANSA')
        nxtEpochComps = nxtEpochComps(2:end);
        nxtEpochs = nxtEpochs(2:end);
    elseif strcmp(tit, 'Pos 2&3')
        nxtEpochComps = nxtEpochComps(2:end);
        nxtEpochs = nxtEpochs(2:end);
        prvEpochComps = prvEpochComps(1:end-1);
        prvEpochs = prvEpochs(1:end-1);
    elseif contains(tit, 'Pos=')
        if sum(posVect~=4)==1
            nxtEpochComps = nxtEpochComps(posVect);
            nxtEpochs = nxtEpochs(posVect);
        else
            nxtEpochComps = [];
            nxtEpochs = [];
        end
        if sum(posVect~=1)==1
            prvEpochComps = prvEpochComps(posVect-1);
            prvEpochs = prvEpochs(posVect-1);
        else
            prvEpochComps = [];
            prvEpochs = [];
        end
    end
    prvEpochs = cell2mat(prvEpochs);
    curEpochs = cell2mat(curEpochs(posVect));
    nxtEpochs = cell2mat(nxtEpochs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Analyses
    figure;
    % Train/Test Within vs Across Epochs
    % Train/Test Post/Pre:
    %       N-1 vs N
    subplot(4,1,1);
    mlb{ani}.PlotMeanVarViolin(2,prvEpochs(3,2,:),3,0.05,'r','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(3,curEpochs(3,1,:),3,0.05,'k','error', 'CI', 'binNum', binNum);
    title([{'Train: Post-Trial'};{'Test: Pre-Trial'}]);
    set(gca, 'xlim', [0 4], 'xtick', 2:3, 'xticklabels', [{'Same ITI'},{'Next ITI'}], 'xticklabelrotation', 45);
    [~,p,~,stats] = ttest2(prvEpochs(3,2,:),curEpochs(3,1,:));
    if p<0.001
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02i', p)}]);
    else
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02f', p)}]);
    end
    % Train/Test Pre/Post:
    %       N+1 vs N
    subplot(4,1,2);
    mlb{ani}.PlotMeanVarViolin(1,curEpochs(1,3,:),3,0.05,'k','error', 'CI', 'binNum', binNum);
%     mlb{ani}.PlotMeanVarViolin(1,prvEpochs(1,3,:),3,0.05,'k','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(2,nxtEpochs(1,3,:),3,0.05,'r','error', 'CI', 'binNum', binNum);
    title([{'Train: Pre-Trial'};{'Test: Post-Trial'}]);
    set(gca, 'xlim', [0 4], 'xtick', 1:2, 'xticklabels', [{'Prev ITI'},{'Same ITI'}], 'xticklabelrotation', 45);
    [~,p,~,stats] = ttest2(curEpochs(1,3,:),nxtEpochs(1,3,:));
    if p<0.001
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02i', p)}]);
    else
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02f', p)}]);
    end
    % Train/Test Pre/Trial:
    %       N-1 vs N vs N+1
    subplot(4,1,3);
    mlb{ani}.PlotMeanVarViolin(1,prvEpochs(1,2,:),3,0.05,'b','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(2,curEpochs(1,2,:),3,0.05,'r','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(3,nxtEpochs(1,2,:),3,0.05,'k','error','CI', 'binNum', binNum);
    anovaDta = [squeeze(prvEpochs(1,2,:)); squeeze(curEpochs(1,2,:)); squeeze(nxtEpochs(1,2,:))];
    grpLog = [ones(size(prvEpochs,3),1); ones(size(curEpochs,3),1)+1; ones(size(nxtEpochs,3),1)+2];
    [~,prTrTbl,prTrStats] = anovan(anovaDta, grpLog);
    prTrMC = multcompare(prTrStats, 'CType', 'bonferroni', 'display', 'off');
    title([{'Train: Pre-Trial'};{'Test:Trial'}]);
    set(gca,'xtick', 1:3, 'xticklabels', [{'Prev'}, {'Upcoming'}, {'Next'}], 'xticklabelrotation', 45);
    % Train/Test Post/Trial:
    %       N-1 vs N vs N+1
    subplot(4,1,4);
    mlb{ani}.PlotMeanVarViolin(1,prvEpochs(3,2,:),3,0.05,'b','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(2,curEpochs(3,2,:),3,0.05,'r','error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(3,nxtEpochs(3,2,:),3,0.05,'k','error','CI', 'binNum', binNum);
    anovaDta = [squeeze(prvEpochs(3,2,:)); squeeze(curEpochs(3,2,:)); squeeze(nxtEpochs(3,2,:))];
    grpLog = [ones(size(prvEpochs,3),1); ones(size(curEpochs,3),1)+1; ones(size(nxtEpochs,3),1)+2];
    [p,poTrTbl,poTrStats] = anovan(anovaDta, grpLog);
    poTrMC = multcompare(poTrStats, 'CType', 'bonferroni', 'display', 'off');
    title([{'Train: Post-Trial'};{'Test:Trial'}]);
    set(gca,'xtick', 1:3, 'xticklabels', [{'Prev'}, {'Recent'}, {'Next'}], 'xticklabelrotation', 45);
        
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s: %s', alignments{al}, tit),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    linkaxes
end

%% Performance Analyses
binNum = 20;

for al = 1:length(alignments)
    curLog = false(size(prvLog,1), size(prvLog,2), length(alignments));
    curLog(:,:,al) = logical(eye(mlb{ani}.seqLength));
    
    tempISCepochs = permute(trlEpochs_tMat(curLog), [2,3,1]);
    tempISCepochs(1) = [];
    tempOSCepochs = permute(trlEpochs_tMat_OSC(curLog), [2,3,1]);
    tempOSCepochs(1) = [];
    tempERRepochs = permute(trlEpochs_tMat_ERR(curLog), [2,3,1]);
    tempERRepochs(1) = [];
    
    preTrialISC = cellfun(@(a){a(1,1,:)}, tempISCepochs);
    preTrialOSC = cellfun(@(a){a(1,1,:)}, tempOSCepochs);
    preTrialERR = cellfun(@(a){a(1,1,:)}, tempERRepochs);
    
    figure;
    subplot(4,1,1);
    mlb{ani}.PlotMeanVarViolin(1, squeeze(cell2mat(preTrialISC)), 1, 0.05, 'k', 'error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(2, squeeze(cell2mat(preTrialOSC)), 1, 0.05, 'r', 'error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(3, squeeze(cell2mat(preTrialERR)), 1, 0.05, 'b', 'error', 'CI', 'binNum', binNum);
    set(gca, 'xtick', 1:3, 'xticklabel', [{'InSeq-Corr'}, {'OutSeq-Corr'}, {'Errors'}], 'xticklabelrotation', 45);
    title('Train/Test Pre-Trial: All Trial Positions');
    
    for pos = 1:size(preTrialISC,3)
        subplot(4,1,pos+1);
        mlb{ani}.PlotMeanVarViolin(1, squeeze(preTrialISC{pos}), 1, 0.05, 'k', 'error', 'CI', 'binNum', binNum);
        mlb{ani}.PlotMeanVarViolin(2, squeeze(preTrialOSC{pos}), 1, 0.05, 'r', 'error', 'CI', 'binNum', binNum);
        mlb{ani}.PlotMeanVarViolin(3, squeeze(preTrialERR{pos}), 1, 0.05, 'b', 'error', 'CI', 'binNum', binNum);
        set(gca, 'xtick', 1:3, 'xticklabel', [{'InSeq-Corr'}, {'OutSeq-Corr'}, {'Errors'}], 'xticklabelrotation', 45);
        title(sprintf('Train/Test Pre-Trial: Pos%i', pos+1));
    end
    linkaxes
end
%% Evaluate Transitions During Trial & ITI Periods
% posVect = 1:mlb{ani}.seqLength; tit = 'All';
posVect = 1:mlb{ani}.seqLength-1; tit = 'All but Final';
% posVect = 2:3;    tit = 'Pos 2&3';
% posVect = 2:mlb{ani}.seqLength;   tit = 'SANSA';
% posVect = 3; tit = sprintf('Pos=%i',posVect);
    
tempTmat = trlD_tMat(:,:,al);
    
for al = 1:length(alignments)
    figure;
    blankTrlLogVect = false(size(obsvTimeVect{al}));
    if strcmp(alignments{al} , 'PokeIn')
        preTrlBound = find(obsvTimeVect{al}==0);%-(binSize/dsRate/2);
        preTrialLog = blankTrlLogVect;
        preTrialLog(1:preTrlBound) = true;
        
        trialLowBound = find(obsvTimeVect{al}==0);%+(binSize/dsRate/2);
        trialUpBound = find(obsvTimeVect{al}<=grpPiPoLat,1,'last');%-(binSize/dsRate/2);
        trialLog = blankTrlLogVect;
        trialLog(trialLowBound:trialUpBound) = true;
        
        pstTrlBound = find(obsvTimeVect{al}<=grpPiPoLat,1,'last');%+(binSize/dsRate/2);
        pstTrialLog = blankTrlLogVect;
        pstTrialLog(pstTrlBound:end) = true;
    elseif strcmp(alignments{al}, 'PokeOut')
        preTrlBound = find(obsvTimeVect{al}<=grpPoPiLat,1,'last');%-(binSize/dsRate/2);
        preTrialLog = blankTrlLogVect;
        preTrialLog(1:preTrlBound) = true;
        
        trialLowBound = find(obsvTimeVect{al}<=grpPoPiLat,1,'last');%+(binSize/dsRate/2);
        trialUpBound = find(obsvTimeVect{al}==0);%-(binSize/dsRate/2);
        trialLog = blankTrlLogVect;
        trialLog(trialLowBound:trialUpBound) = true;
        
        pstTrlBound = find(obsvTimeVect{al}==0);%+(binSize/dsRate/2);
        pstTrialLog = blankTrlLogVect;
        pstTrialLog(pstTrlBound:end) = true;
    end
    
    iscTmat = tempTmat(logical(eye(length(tempTmat))));
    iscTmatTrials = cell2mat(permute(iscTmat(posVect),[2,3,1]));
    nxtTmat = tempTmat(logical((triu(ones(length(tempTmat)),1) & tril(ones(length(tempTmat)),1))));
    prevTmat = tempTmat(logical((triu(ones(length(tempTmat)),-1) & tril(ones(length(tempTmat)),-1))));
    if strcmp(tit, 'SANSA')
        nxtTmat = nxtTmat(2:end);
    elseif strcmp(tit, 'Pos 2&3')
        nxtTmat = nxtTmat(2:end);
        prevTmat = prevTmat(1:end-1);
    elseif contains(tit, 'Pos=')
        if posVect~=4
            nxtTmat = nxtTmat(posVect);
        else
            nxtTmat = [];
        end
        if posVect~=1
            prevTmat = prevTmat(posVect-1);
        else
            prevTmat = [];
        end
    end
    nxtTmatTrials = cell2mat(permute(nxtTmat,[2,3,1]));
    prevTmatTrials = cell2mat(permute(prevTmat,[2,3,1]));
    
        % Trial Trained Models
    % Pre-trial period    
    prevTrPr = prevTmatTrials(preTrialLog,trialLog,:);
    trTr_TsPr = iscTmatTrials(preTrialLog,trialLog,:);
    nxtTrPr = nxtTmatTrials(preTrialLog,trialLog,:);
        
    prevTrPo = prevTmatTrials(pstTrialLog,trialLog,:);
    trTr_TsPo = iscTmatTrials(pstTrialLog,trialLog,:);
    nxtTrPo = nxtTmatTrials(pstTrialLog,trialLog,:);    
    
    trTr_TsPr_TrlMean = nan(1,size(trTr_TsPr,3));
    trTr_TsPo_TrlMean = nan(1,size(trTr_TsPo,3));
    for trl = 1:size(trTr_TsPr,3)
        trTr_TsPr_TrlMean(trl) = mean(mean(trTr_TsPr(:,:,trl),1));
        trTr_TsPo_TrlMean(trl) = mean(mean(trTr_TsPo(:,:,trl),1));
    end
        
    subplot(6,4,[1,5]);
    prevTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(prevTrPr,2),3,0.05,'b');
    curTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(trTr_TsPr,2),3,0.05,'r'); 
    nxtTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(nxtTrPr,2),3,0.05,'g');
    legend([prevTP, curTP,nxtTP], [{'Prev'}, {'Current'}, {'Next'}]);
%     legend([prevTP,nxtTP], [{'Prev'}, {'Next'}]);
    subplot(6,4,9);
    mlb{ani}.PlotMeanVarViolin(0, trTr_TsPr_TrlMean, 2, 0.05, 'k', 'error', 'CI', 'binNum', binNum);
    mlb{ani}.PlotMeanVarViolin(1, trTr_TsPo_TrlMean, 2, 0.05, 'r', 'error', 'CI', 'binNum', binNum);
    [~,p,~,stats] = ttest2(trTr_TsPr_TrlMean,trTr_TsPo_TrlMean);
    if p<0.01
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02i', p)}]);
    else
        text(min(get(gca,'xlim')), max(get(gca, 'ylim')), [{sprintf('t=%.02f', stats.tstat)};{sprintf('p=%.02f', p)}]);
    end

    % Trial Period
    subplot(6,4,[2,3,6,7]);
    pr = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trTr_TsPr,1),3,0.05,'k');
    hold on;
    pst = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trTr_TsPo,1),3,0.05,'r');
    legend([pr pst], [{'Pre-Trial'},{'Post-Trial'}], 'location', 'northwest');
    title('Decoding Outside the Trial Period (Train on Trial)');
    subplot(6,4,[10,11]);
    mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trTr_TsPo,1)-mean(trTr_TsPr,1),3,0.05,'k');
    plot(get(gca,'xlim'), [0 0], '-k', 'linewidth', 1);
    grid on;
    % Post-Trial Period
    subplot(6,4,[4,8]);
    prevTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(prevTrPo,2),3,0.05,'b');
    curTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(trTr_TsPo,2),3,0.05,'r'); 
    nxtTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(nxtTrPo,2),3,0.05,'g');
%     legend([prevTP, curTP,nxtTP], [{'Prev'}, {'Current'}, {'Next'}]);
        % ITD Trained Models
    %Pre-trial period
    prevPrTr = prevTmatTrials(trialLog,preTrialLog,:);
    trPr_TsTrl = iscTmatTrials(trialLog,preTrialLog,:);
    nxtPrTr = nxtTmatTrials(trialLog,preTrialLog,:);
    
    prevPoTr = prevTmatTrials(trialLog,pstTrialLog,:);
    trPo_TsTrl = iscTmatTrials(trialLog,pstTrialLog,:);
    nxtPoTr = nxtTmatTrials(trialLog,pstTrialLog,:);    
    
    subplot(6,4,[13,17]);
    prevTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(prevPrTr),3,0.05,'b');
    curTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(trPr_TsTrl),3,0.05,'r'); 
    nxtTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(nxtPrTr),3,0.05,'g');
    legend([prevTP, curTP,nxtTP], [{'Prev'}, {'Current'}, {'Next'}]);
%     legend([prevTP,nxtTP], [{'Prev'}, {'Next'}]);
    subplot(6,4,[14,15,18,19]);
    pr = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trPr_TsTrl,2),3,0.05,'k');
    hold on;
    pst = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trPo_TsTrl,2),3,0.05,'r');
    legend([pr pst], [{'Pre-Trial'},{'Post-Trial'}], 'location', 'northwest');
    title('Decoding the Trial Period (Train on ITD)');
    subplot(6,4,[22,23]);
    mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(trialLowBound):dsRate/1000:obsvTimeVect{al}(trialUpBound),mean(trPo_TsTrl,2)-mean(trPr_TsTrl,2),3,0.05,'k');
    plot(get(gca,'xlim'), [0 0], '-k', 'linewidth', 1);
    grid on;
    xlabel('Trial Time');
    % Post-Trial Period
    subplot(6,4,[16,20]);
    prevTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(prevPoTr),3,0.05,'b');
    curTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(trPo_TsTrl),3,0.05,'r'); 
    nxtTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(pstTrlBound):dsRate/1000:obsvTimeVect{al}(end),mean(nxtPoTr),3,0.05,'g');
%     legend([prevTP, curTP,nxtTP], [{'Prev'}, {'Current'}, {'Next'}]);
    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s: %s', alignments{al}, tit),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    linkaxes
end

%% Plot things
% Code taken from above when it was integrated into things... rework it so things look nice and work here too

    norm = figure;
    time = figure; 
    subplot(mlb{ani}.seqLength, mlb{ani}.seqLength, sub2ind([mlb{ani}.seqLength, mlb{ani}.seqLength],posD, posT));
            mlb{end}.PlotTrialPDM(trlD_tMat{posT,posD,al}, 'rotate', 'clim', [-1 1], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':w','linewidth', 1);

                plot(zeros(1,2),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':w','linewidth', 1);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--w','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 1);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--w','linewidth', 2);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':w','linewidth', 1);
            end
            title(sprintf('Pos=%i; Dec=%i',posT,posD));
            figure(time);
            subplot(mlb{ani}.seqLength, mlb{ani}.seqLength, sub2ind([mlb{ani}.seqLength, mlb{ani}.seqLength],posD, posT));
            mlb{end}.PlotTrialPDM(trlDtime_tMat{posT,posD,al}, 'rotate', 'clim', [-1 1], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':w','linewidth', 1);

                plot(zeros(1,2),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':w','linewidth', 1);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--w','linewidth', 1);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':w','linewidth', 1);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--w','linewidth', 1);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':w','linewidth', 1);
            end
            title(sprintf('Pos=%i; Dec=%i',posT,posD));
            
    annotation(norm,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment (Session FAR)', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(time,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment (Local FAR)', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');



                
    