% PFC_XTD_TrialWindows_Group_PosChance

%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'},...
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

% alignments = [{'PokeIn'}, {'PokeOut'}];
% trlWindows = [{[-1500 2000]}, {[-1500 1500]}];
% trlWindows = [{[-700 2000]}, {[-2000 700]}]; %% Non-Overlapping Pre/Post ITD

% alignments = {'PokeOut'};
% trlWindows = {[-1500 2000]};
% trlWindows = {[-2000 700]}; %% Non-Overlapping Pre/Post ITD
alignments = {'PokeIn'};
% trlWindows = {[-1500 2000]};
% trlWindows = {[-1200 2500]}; %% Overlapping Pre/Post ITD
trlWindows = {[-700 2000]}; %% Non-Overlapping Pre/Post ITD
    
lfpWindow = [16 32];
numChancePerms = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Chance Perm Nums

mlb = cell(size(fileDirs));
%%
realTic = tic;
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb{ani} = MLB_SM(fileDirs{ani});
    % Create Analysis Variables
    if ani == 1 
        % Behavior Variables
        piPokeOutLat = cell(length(fileDirs),1);
        piRwdLat = cell(length(fileDirs),1);
        poPokeInLat = cell(length(fileDirs),1);
        poRwdLat = cell(length(fileDirs),1);
        nxtTrlLat = cell(length(fileDirs),1);
        iscTAOlog = cell(mlb{ani}.seqLength,mlb{ani}.seqLength,length(fileDirs));
        iscSeqOSlog = cell(mlb{ani}.seqLength,mlb{ani}.seqLength,length(fileDirs));
        % Posteriors
        trlD = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        hitRate = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        faRate = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPr_TsTr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrTr_TsPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrTr_TsPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPo_TsTr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPo_TsPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPr_TsPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_PoPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrTr_TsTr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrTr_TsTr_Diag = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrTr_TsTr_OffDiag = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPr_TsPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrPo_TsPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        trlPersFit = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrITD_TsITD = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrITD_TsITD_Diag = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrITD_TsITD_OffDiag = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        window_TrlAltITD_TsAltITD = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        itiPersFit = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
        obsvTimeVect = cell(1,length(alignments));
        % LFP
        betaPower = cell(mlb{ani}.seqLength, length(alignments), length(fileDirs));
        thetaPower = cell(mlb{ani}.seqLength, length(alignments), length(fileDirs));
        trlMeanBeta = cell(mlb{ani}.seqLength, length(alignments), length(fileDirs));
        trlMeanTheta = cell(mlb{ani}.seqLength, length(alignments), length(fileDirs));
    end
    mlb{ani}.binSize = binSize;
    mlb{ani}.dsRate = dsRate;
    mlb{ani}.bayesType = bayesType;
   
    %% Extract Behavioral Variables
    piPokeOutLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex])'/1000;
    piRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex])'/1000;
    poPokeInLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeInIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex])'/1000;
    poRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials(~isnan(mlb{ani}.fiscTrials))).PokeOutIndex])'/1000;
%     nxtTrlLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(1:end-1,:)+1).PokeInIndex]-[mlb{ani}.trialInfo(mlb{ani}.fiscTrials(1:end-1,:)).PokeOutIndex])'/1000;
    
    %% Process Neural Data
    for al = 1:length(alignments)
        mlb{ani}.alignments = alignments(al);
        mlb{ani}.windows = trlWindows(al);
        mlb{ani}.SetLikes_ISC;
        mlb{ani}.Process_IterativeLikelyL1O;
        [~, tempBetaPower] = mlb{ani}.PP_TrialMatrix_LFP([16 32], trlWindows{al}, alignments{al});
        [~, tempThetaPower] = mlb{ani}.PP_TrialMatrix_LFP([4 12], trlWindows{al}, alignments{al});
        trlIDs = mlb{ani}.postTrlIDs(:);
        trlIDlog = ~isnan(trlIDs);
        posts = permute(mlb{ani}.post(:),[3,2,1]);
        posts = posts(:,:,trlIDlog);
        trlIDs = trlIDs(trlIDlog);
        trlPosVect = [mlb{ani}.trialInfo(trlIDs).Position];
        taoLog = false(size(trlIDs));
        seqOSlog = false(size(trlIDs));
        for trl = 2:length(trlIDs)
            if mlb{ani}.trialInfo(trlIDs(trl)-1).TranspositionDistance~=0 &&...
                    mlb{ani}.trialInfo(trlIDs(trl)-1).Performance == 1 &&...
                    mlb{ani}.trialInfo(trlIDs(trl)-1).Position ~=4
                taoLog(trl) = true;
            end
            if sum([mlb{ani}.trialInfo([mlb{ani}.trialInfo.SequenceNum]==mlb{ani}.trialInfo(trlIDs(trl)).SequenceNum & [mlb{ani}.trialInfo.TrialNum]<mlb{ani}.trialInfo(trlIDs(trl)).TrialNum).TranspositionDistance])~=0
                seqOSlog(trl) = true;
            end                
        end
        if ani == 1
            obsvTimeVect{al} = mlb{ani}.obsvTimeVect;
        end
        % Pull out replacement values to replace instances where mean posterior is P==0 or P=1
        unq = unique(cell2mat(posts));
        if unq(1)==0
            minProb = unq(2);
        else
            minProb = unq(1);
        end
        if unq(end)==1
            maxProb = unq(end-1);
        else
            maxProb = unq(end);
        end
        
        for pos = 1:mlb{ani}.seqLength
            % Extract out behavioral events for windowing purposes
            curTrlIDs = trlIDs(trlPosVect==pos);
            curTrlTAOlog = taoLog(trlPosVect==pos);
            curTrlSeqOSlog = seqOSlog(trlPosVect==pos);
            betaPower{pos,al,ani} = tempBetaPower(:,:,curTrlIDs);
            thetaPower{pos,al,ani} = tempThetaPower(:,:,curTrlIDs);
            curTrlPokeDurs = [mlb{ani}.trialInfo(curTrlIDs).PokeDuration];
            if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                nearestPIval = zeros(size(curTrlPokeDurs));
                nearestPOval = cell2mat(arrayfun(@(a){mlb{ani}.obsvTimeVect(find(a>mlb{ani}.obsvTimeVect,1,'last'))}, curTrlPokeDurs));
            elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                nearestPIval = cell2mat(arrayfun(@(a){mlb{ani}.obsvTimeVect(find(a>mlb{ani}.obsvTimeVect,1,'last'))}, curTrlPokeDurs))*-1;
                nearestPOval = zeros(size(curTrlPokeDurs));
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% +0.5 is a test!!!!
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Don't leave this here!!!
%                 nearestPOval = zeros(size(curTrlPokeDurs))+0.5;
            end
            % Now calculate decodability plots
            for p = 1:mlb{ani}.seqLength
                if al == 1
                    iscTAOlog{pos,p,ani} = curTrlTAOlog;
                    iscSeqOSlog{pos,p,ani} = curTrlSeqOSlog;
                end
                hits = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect==pos)));
                hitRate{pos,p,ani,al} = hits;
%                 if pos==p
                    fas = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect~=pos)));
%                 else
%                     fas = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect~=pos & trlPosVect~=p)));
%                 end
                faRate{pos,p,ani,al} = fas;
                meanFAR = mean(fas,3);
                meanFAR(meanFAR==0) = minProb;
                meanFAR(meanFAR==1) = maxProb;
                tempTrlD = nan(size(hits));                                     % Trial D matrix
                tempWindow_TrPr_TsTr = nan(size(hits,3),1);                     % Train Pre-Trial & Test Trial Window
                tempWindow_TrTr_TsPr = nan(size(hits,3),1);                     % Train Trial & Test Pre-Trial Window
                tempWindow_TrPr = nan(size(hits,3),1);                          % Pre-Trial vs Test Windows
                tempWindow_TrTr_TsPo = nan(size(hits,3),1);                     % Train Trial & Test Post-Trial Window
                tempWindow_TrPo_TsTr = nan(size(hits,3),1);                     % Train Post-Trial & Test Trial Window
                tempWindow_TrPo = nan(size(hits,3),1);                          % Trial vs Post-Trial Windows
                tempWindow_TrPo_TsPr = nan(size(hits,3),1);                     % Train Post-Trial & Test Pre-Trial Window
                tempWindow_TrPr_TsPo = nan(size(hits,3),1);                     % Train Pre-Trial & Test Post-Trial Window
                tempWindow_PoPr = nan(size(hits,3),1);                          % Post-Trial vs Pre-Trial Windows
                tempWindow_TrTr_TsTr = nan(size(hits,3),1);                     % Train & Test Trial Window
                tempWindow_TrTr_TsTr_Diag = nan(size(hits,3),1);                % Diagonal of the Train & Test Trial Window
                tempWindow_TrTr_TsTr_OffDiag = nan(size(hits,3),1);             % Off Diagonal of the Train & Test Trial Window
                tempWindow_TrPr_TsPr = nan(size(hits,3),1);                     % Train & Test Pre-Trial Window
                tempWindow_TrPo_TsPo = nan(size(hits,3),1);                     % Train & Test Post-Trial Window
                tempWindow_TrITD_TsITD = nan(size(hits,3),1);                   % Train & Test ITD Window
                tempWindow_TrITD_TsITD_Diag = nan(size(hits,3),1);              % Diagonal of the ITD Window
                tempWindow_TrITD_TsITD_OffDiag = nan(size(hits,3),1);           % Off Diagonal of the Train & Test ITD Window
                tempWindow_TrAltITD_TsAltITD = nan(size(hits,3),1);            % Train & Test on the Post-Trial Window
                tempTrialFit = nan(size(hits,3), length(mlb{ani}.obsvTimeVect));     % Persistence model fit Trial period
                tempITIfit = nan(size(hits,3), length(mlb{ani}.obsvTimeVect));       % Persistence model fit ITI period
                if p == pos
                    tempTrlBeta = nan(size(hits,3),1);
                    tempTrlTheta = nan(size(hits,3),1);
                end
                for trl = 1:size(hits,3)
                    %% Calculate D Matrix
                    tempHR = hits(:,:,trl);
                    tempHR(tempHR==0) = minProb;
                    tempHR(tempHR==1) = maxProb;
                    curTrlD = arrayfun(@(a,b)norminv(a)-norminv(b), tempHR, meanFAR);
                    tempTrlD(:,:,trl) = curTrlD;
                    
                    % Define logical vectors for trial periods
                    preTrialLog = mlb{ani}.obsvTimeVect<nearestPIval(trl);
                    trialLog = mlb{ani}.obsvTimeVect>=nearestPIval(trl) & mlb{ani}.obsvTimeVect<nearestPOval(trl);
                    pstTrialLog = mlb{ani}.obsvTimeVect>=nearestPOval(trl);
                    if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                        itiLog = preTrialLog;
                        altITIlog = pstTrialLog;
                    elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                        itiLog = pstTrialLog;
                        altITIlog = preTrialLog;
                    end
                                        
                    %% Quantify mean windows
                    % IMPORTANT NOTE: The organization of the posterior data is dim1/y-axis/rows = testing time, dim2/x-axis/columns = training time
                    % The organization of most presentations is reversed, i.e. train=rows; test=columns;
                    blankMask = false(size(curTrlD));
                    % Train Pre-Trial & Test Trial                    
                    prtrWindow = blankMask;
                    prtrWindow(trialLog,preTrialLog) = true;
                    tempWindow_TrPr_TsTr(trl) = mean(curTrlD(prtrWindow));
                    % Train Trial & Test Pre-Trial
                    trprWindow = blankMask;
                    trprWindow(preTrialLog,trialLog) = true;
                    tempWindow_TrTr_TsPr(trl) = mean(curTrlD(trprWindow));
                    % Pre-Trial vs Trial
                    trprprtrWindow = prtrWindow | trprWindow;
                    tempWindow_TrPr(trl) = mean(curTrlD(trprprtrWindow));
                    % Train Trial & Test Post-Trial
                    trpoWindow = blankMask;
                    trpoWindow(pstTrialLog,trialLog) = true;
                    tempWindow_TrTr_TsPo(trl) = mean(curTrlD(trpoWindow));
                    % Train Post-Trial & Test Trial
                    potrWindow = blankMask;
                    potrWindow(trialLog,pstTrialLog) = true;
                    tempWindow_TrPo_TsTr(trl) = mean(curTrlD(potrWindow));
                    % Trial vs Post-Trial
                    trpopotrWindow = trpoWindow | potrWindow;
                    tempWindow_TrPo(trl) = mean(curTrlD(trpopotrWindow));
                    % Train Post-Trial & Test Pre-Trial
                    poprWindow = blankMask;
                    poprWindow(preTrialLog, pstTrialLog) = true;
                    tempWindow_TrPo_TsPr(trl) = mean(curTrlD(poprWindow));
                    % Train Pre-Trial & Test Post-Trial
                    prpoWindow = blankMask;
                    prpoWindow(pstTrialLog, preTrialLog) = true;
                    tempWindow_TrPr_TsPo(trl) = mean(curTrlD(prpoWindow));
                    % Post-Trial vs Pre-Trial
                    poprprpoWindow = poprWindow | prpoWindow;
                    tempWindow_PoPr(trl) = mean(curTrlD(poprprpoWindow));      
                    
                    % Trial Window
                    trlWindow = blankMask;
                    trlWindow(trialLog, trialLog) = true;
                    tempWindow_TrTr_TsTr(trl) = mean(curTrlD(trlWindow));
                    
                    trlEpoch = curTrlD(trialLog, trialLog);
                    diagTrlEpoch = triu(ones(size(trlEpoch)), (binSize/2)/dsRate*-1) & tril(ones(size(trlEpoch)), (binSize/2)/dsRate);
                    tempWindow_TrTr_TsTr_Diag(trl) = mean(trlEpoch(diagTrlEpoch));
                    tempWindow_TrTr_TsTr_OffDiag(trl) = mean(trlEpoch(~diagTrlEpoch));
                    
                    % Pre-Trial Window
                    prTrlWindow = blankMask;
                    prTrlWindow(preTrialLog, preTrialLog) = true;
                    tempWindow_TrPr_TsPr(trl) = mean(curTrlD(prTrlWindow));
                    
                    % Post-Trial Window
                    poTrlWindow = blankMask;
                    poTrlWindow(pstTrialLog, pstTrialLog) = true;
                    tempWindow_TrPo_TsPo(trl) = mean(curTrlD(poTrlWindow));
                    
                    % ITD Window
                    itdWindow = blankMask;
                    itdWindow(itiLog, itiLog) = true;
                    tempWindow_TrITD_TsITD(trl) = mean(curTrlD(itdWindow));
                    
                    % Alternate ITD Window
                    altITDwindow = blankMask;
                    altITDwindow(altITIlog, altITIlog) = true;
                    tempWindow_TrAltITD_TsAltITD(trl) = mean(curTrlD(altITDwindow));
                    
                    itdEpoch = curTrlD(itiLog, itiLog);
                    diagITDepoch = triu(ones(size(itdEpoch)), (binSize/2)/dsRate*-1) & tril(ones(size(itdEpoch)), (binSize/2)/dsRate);
                    tempWindow_TrITD_TsITD_Diag(trl) = mean(itdEpoch(diagITDepoch));
                    tempWindow_TrITD_TsITD_OffDiag(trl) = mean(itdEpoch(~diagITDepoch));
                    
                    %% Fit the Trial & ITI periods to dynamic coding "models"
                    trialDynCap = sum(trialLog)-1;
                    trialPeriod = zscore(curTrlD(trialLog,trialLog), 0, 'all');
                    for d = 1:trialDynCap-1
                        tempDyn = zscore(triu(true(sum(trialLog)), d*-1) & tril(true(sum(trialLog)), d), 0, 'all');
                        tempTrialFit(trl,d) = pdist([trialPeriod(:)';tempDyn(:)'], 'cosine');
                    end
                    itiDynCap = sum(itiLog)-1;
                    itiPeriod = zscore(curTrlD(itiLog,itiLog), 0 , 'all');
                    for d = 1:itiDynCap-1
                        tempDyn = zscore(triu(true(sum(itiLog)), d*-1) & tril(true(sum(itiLog)), d), 0, 'all');
                        tempITIfit(trl,d) = pdist([itiPeriod(:)';tempDyn(:)'], 'cosine');
                    end
                    
                    %% Calculate average LFP power (if applicable)
                    if p == pos
                        tempTrlBeta(trl) = mean(betaPower{pos,al,ani}(trialLog,1,trl));
                        tempTrlTheta(trl) = mean(thetaPower{pos,al,ani}(trialLog,1,trl));
                    end
                end
                trlD{pos,p,ani,al} = tempTrlD;                
                window_TrPr_TsTr{pos,p,ani,al} = tempWindow_TrPr_TsTr;
                window_TrTr_TsPr{pos,p,ani,al} = tempWindow_TrTr_TsPr;
                window_TrPr{pos,p,ani,al} = tempWindow_TrPr;
                window_TrTr_TsPo{pos,p,ani,al} = tempWindow_TrTr_TsPo;
                window_TrPo_TsTr{pos,p,ani,al} = tempWindow_TrPo_TsTr;
                window_TrPo{pos,p,ani,al} = tempWindow_TrPo;
                window_TrPo_TsPr{pos,p,ani,al} = tempWindow_TrPo_TsPr;
                window_TrPr_TsPo{pos,p,ani,al} = tempWindow_TrPr_TsPo;
                window_PoPr{pos,p,ani,al} = tempWindow_PoPr;
                window_TrTr_TsTr{pos,p,ani,al} = tempWindow_TrTr_TsTr;
                window_TrTr_TsTr_Diag{pos,p,ani,al} = tempWindow_TrTr_TsTr_Diag;
                window_TrTr_TsTr_OffDiag{pos,p,ani,al} = tempWindow_TrTr_TsTr_OffDiag;
                window_TrPr_TsPr{pos,p,ani,al} = tempWindow_TrPr_TsPr;
                window_TrPo_TsPo{pos,p,ani,al} = tempWindow_TrPo_TsPo;
                window_TrITD_TsITD{pos,p,ani,al} = tempWindow_TrITD_TsITD;
                window_TrITD_TsITD_Diag{pos,p,ani,al} = tempWindow_TrITD_TsITD_Diag;
                window_TrITD_TsITD_OffDiag{pos,p,ani,al} = tempWindow_TrITD_TsITD_OffDiag;
                window_TrlAltITD_TsAltITD{pos,p,ani,al} = tempWindow_TrAltITD_TsAltITD;
                trlPersFit{pos,p,ani,al} = tempTrialFit;
                itiPersFit{pos,p,ani,al} = tempITIfit;
                if p==pos
                    trlMeanBeta{pos,al,ani} = tempTrlBeta;
                    trlMeanTheta{pos,al,ani} = tempTrlTheta;
                end
            end
        end
    end
end

%% Collapse trial data across animals
trialD = cell(size(trlD,1),size(trlD,2),size(trlD,4));
trialDz = cell(size(trlD,1),size(trlD,2),size(trlD,4));
trialHR = cell(size(hitRate,1),size(hitRate,2), size(hitRate,4));
trialFAR = cell(size(faRate,1),size(faRate,2),size(faRate,4));
trial_TAOlog = cell(size(iscTAOlog(:,:,1)));
trial_SeqOSlog = cell(size(iscSeqOSlog(:,:,1)));
trial_Window_TrPr_TsTr = cell(size(window_TrPr_TsTr,1), size(window_TrPr_TsTr,2), size(window_TrPr_TsTr,4));
trial_Window_TrTr_TsPr = cell(size(window_TrTr_TsPr,1), size(window_TrTr_TsPr,2), size(window_TrTr_TsPr,4));
trial_Window_TrPr = cell(size(window_TrPr,1), size(window_TrPr,2), size(window_TrPr,4));
trial_Window_TrTr_TsPo = cell(size(window_TrTr_TsPo,1), size(window_TrTr_TsPo,2), size(window_TrTr_TsPo,4));
trial_Window_TrPo_TsTr = cell(size(window_TrPo_TsTr,1), size(window_TrPo_TsTr,2), size(window_TrPo_TsTr,4));
trial_Window_TrPo = cell(size(window_TrPo,1), size(window_TrPo,2), size(window_TrPo,4));
trial_Window_TrPo_TsPr = cell(size(window_TrPo_TsPr,1), size(window_TrPo_TsPr,2), size(window_TrPo_TsPr,4));
trial_Window_TrPr_TsPo = cell(size(window_TrPr_TsPo,1), size(window_TrPr_TsPo,2), size(window_TrPr_TsPo,4));
trial_Window_PoPr = cell(size(window_PoPr,1), size(window_PoPr,2), size(window_PoPr,4));
trial_Window_TrTr_TsTr = cell(size(window_TrTr_TsTr,1), size(window_TrTr_TsTr,2), size(window_TrTr_TsTr,4));
trial_Window_TrTr_TsTr_Diag = cell(size(window_TrTr_TsTr_Diag,1), size(window_TrTr_TsTr_Diag,2), size(window_TrTr_TsTr_Diag,4));
trial_Window_TrTr_TsTr_OffDiag = cell(size(window_TrTr_TsTr_OffDiag,1), size(window_TrTr_TsTr_OffDiag,2), size(window_TrTr_TsTr_OffDiag,4));
trial_Window_TrPr_TsPr = cell(size(window_TrPr_TsPr,1), size(window_TrPr_TsPr,2), size(window_TrPr_TsPr,4));
trial_Window_TrPo_TsPo = cell(size(window_TrPo_TsPo,1), size(window_TrPo_TsPo,2), size(window_TrPo_TsPo,4));
trial_Window_TrITD_TsITD = cell(size(window_TrITD_TsITD,1), size(window_TrITD_TsITD,2), size(window_TrITD_TsITD,4));
trial_Window_TrITD_TsITD_Diag = cell(size(window_TrITD_TsITD_Diag,1), size(window_TrITD_TsITD_Diag,2), size(window_TrITD_TsITD_Diag,4));
trial_Window_TrITD_TsITD_OffDiag = cell(size(window_TrITD_TsITD_OffDiag,1), size(window_TrITD_TsITD_OffDiag,2), size(window_TrITD_TsITD_OffDiag,4));
trial_Window_TrAltITD_TsAltITD = cell(size(window_TrlAltITD_TsAltITD,1), size(window_TrlAltITD_TsAltITD,2), size(window_TrlAltITD_TsAltITD,4));
trial_TrialPersFit = cell(size(trlPersFit,1), size(trlPersFit,2), size(trlPersFit,4));
trial_ITIpersFit = cell(size(itiPersFit,1), size(itiPersFit,2), size(itiPersFit,4));
trial_BetaPowerTime = cell(size(betaPower,1), size(betaPower,2));
trial_ThetaPowerTime = cell(size(thetaPower,1), size(thetaPower,2));
trial_MeanBeta = cell(size(trlMeanBeta,1), size(trlMeanBeta,2));
trial_MeanTheta = cell(size(trlMeanTheta,1), size(trlMeanTheta,2));
for r = 1:mlb{1}.seqLength
    for al = 1:length(alignments)
        trial_BetaPowerTime{r,al} = cell2mat(betaPower(r,al,:));
        trial_ThetaPowerTime{r,al} = cell2mat(thetaPower(r,al,:));
        trial_MeanBeta{r,al} = cell2mat(permute(trlMeanBeta(r,al,:), [3,1,2]));
        trial_MeanTheta{r,al} = cell2mat(permute(trlMeanTheta(r,al,:), [3,1,2]));
    end
    for c = 1:mlb{1}.seqLength
        for al = 1:length(alignments)
            trialD{r,c,al} = cell2mat(trlD(r,c,:,al));  
            temp = nan(size(trialD{r,c,al}));
            for t = 1:size(trialD{r,c,al},3)
                temp(:,:,t) = zscore(trialD{r,c,al}(:,:,t),0,'all');
            end
            trialDz{r,c,al} = temp;
            trialHR{r,c,al} = cell2mat(hitRate(r,c,:,al));
            trialFAR{r,c,al} = cell2mat(faRate(r,c,:,al));
            trial_TAOlog{r,c} = cell2mat(permute(iscTAOlog(r,c,:), [3,1,2]));
            trial_SeqOSlog{r,c} = cell2mat(permute(iscSeqOSlog(r,c,:), [3,1,2]));
            trial_Window_TrPr_TsTr{r,c,al} = cell2mat(permute(window_TrPr_TsTr(r,c,:,al), [3,1,2]));            
            trial_Window_TrTr_TsPr{r,c,al} = cell2mat(permute(window_TrTr_TsPr(r,c,:,al), [3,1,2]));            
            trial_Window_TrTr_TsPo{r,c,al} = cell2mat(permute(window_TrTr_TsPo(r,c,:,al), [3,1,2]));
            trial_Window_TrPr{r,c,al} = cell2mat(permute(window_TrPr(r,c,:,al), [3,1,2]));
            trial_Window_TrPo_TsTr{r,c,al} = cell2mat(permute(window_TrPo_TsTr(r,c,:,al), [3,1,2]));
            trial_Window_TrPo{r,c,al} = cell2mat(permute(window_TrPo(r,c,:,al), [3,1,2]));
            trial_Window_TrPo_TsPr{r,c,al} = cell2mat(permute(window_TrPo_TsPr(r,c,:,al), [3,1,2]));
            trial_Window_TrPr_TsPo{r,c,al} = cell2mat(permute(window_TrPr_TsPo(r,c,:,al), [3,1,2]));
            trial_Window_PoPr{r,c,al} = cell2mat(permute(window_PoPr(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr{r,c,al} = cell2mat(permute(window_TrTr_TsTr(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr_Diag{r,c,al} = cell2mat(permute(window_TrTr_TsTr_Diag(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr_OffDiag{r,c,al} = cell2mat(permute(window_TrTr_TsTr_OffDiag(r,c,:,al), [3,1,2]));
            trial_Window_TrPr_TsPr{r,c,al} = cell2mat(permute(window_TrPr_TsPr(r,c,:,al), [3,1,2]));
            trial_Window_TrPo_TsPo{r,c,al} = cell2mat(permute(window_TrPo_TsPo(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD{r,c,al} = cell2mat(permute(window_TrITD_TsITD(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD_Diag{r,c,al} = cell2mat(permute(window_TrITD_TsITD_Diag(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD_OffDiag{r,c,al} = cell2mat(permute(window_TrITD_TsITD_OffDiag(r,c,:,al), [3,1,2]));
            trial_Window_TrAltITD_TsAltITD{r,c,al} = cell2mat(permute(window_TrlAltITD_TsAltITD(r,c,:,al), [3,1,2]));
            trial_TrialPersFit{r,c,al} = cell2mat(permute(trlPersFit(r,c,:,al), [3,1,2]));
            trial_ITIpersFit{r,c,al} = cell2mat(permute(itiPersFit(r,c,:,al), [3,1,2]));
        end
    end
end

trialDposZ = cell(size(trialD));
for p = 1:mlb{ani}.seqLength
    meanHRd = mean(trialD{p,p}(:));
    varHRd = std(trialD{p,p}(:));
    trialDposZ{p,p} = (trialD{p,p}-meanHRd)./varHRd;
    
    offDiag = find((1:mlb{ani}.seqLength)~=p);
    for od = 1:length(offDiag)
        trialDposZ{p,offDiag(od)} = (trialD{p,offDiag(od)}-meanHRd)./varHRd;
    end    
end

    
fprintf('Real data compiled after %.02fmin\n', toc(realTic)/60);

%% Z-Score Everything w/in windows
trialWiseRawEpochs = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
for al = 1:length(alignments)
    for p = 1:mlb{ani}.seqLength
        for p2 = 1:mlb{ani}.seqLength
            tempTrialWindowed = [trial_Window_TrPr_TsTr{p,p2,al}, trial_Window_TrTr_TsPr{p,p2,al},...
                trial_Window_TrTr_TsPo{p,p2,al}, trial_Window_TrPo_TsTr{p,p2,al},...
                trial_Window_TrPo_TsPr{p,p2,al}, trial_Window_TrPr_TsPo{p,p2,al},...
                trial_Window_TrTr_TsTr{p,p2,al}, trial_Window_TrITD_TsITD{p,p2,al},...
                trial_Window_TrAltITD_TsAltITD{p,p2,al}, trial_Window_TrPr_TsPr{p,p2,al},...
                trial_Window_TrPo_TsPo{p,p2,al}];
            trialWiseRawEpochs{p,p2,al} = tempTrialWindowed;
            [r,c] = find(isnan(tempTrialWindowed));
            colMean = mean(tempTrialWindowed,'omitnan');
            for val = 1:length(r)
                tempTrialWindowed(r(val),c(val)) = colMean(c(val));
            end
            tempTrialWindowed = zscore(tempTrialWindowed,0,'all');
            trial_Window_TrPr_TsTr{p,p2,al} = tempTrialWindowed(:,1);
            trial_Window_TrTr_TsPr{p,p2,al} = tempTrialWindowed(:,2);
            trial_Window_TrTr_TsPo{p,p2,al} = tempTrialWindowed(:,3);
            trial_Window_TrPo_TsTr{p,p2,al} = tempTrialWindowed(:,4);
            trial_Window_TrPo_TsPr{p,p2,al} = tempTrialWindowed(:,5);
            trial_Window_TrPr_TsPo{p,p2,al} = tempTrialWindowed(:,6);
            trial_Window_TrTr_TsTr{p,p2,al} = tempTrialWindowed(:,7);
            trial_Window_TrITD_TsITD{p,p2,al} = tempTrialWindowed(:,8);
            trial_Window_TrAltITD_TsAltITD{p,p2,al} = tempTrialWindowed(:,9);
            trial_Window_TrPr_TsPr{p,p2,al} = tempTrialWindowed(:,10);
            trial_Window_TrPo_TsPo{p,p2,al} = tempTrialWindowed(:,11);
        end
    end
end
%% Calculate Chance
% Chance Values
trialD_Chance = cell(size(trlD,1),size(trlD,2),numChancePerms,size(trlD,4));
trial_Window_TrPr_TsTr_Chance = cell(size(window_TrPr_TsTr,1), size(window_TrPr_TsTr,2), numChancePerms, size(window_TrPr_TsTr,4));
trial_Window_TrTr_TsPr_Chance = cell(size(window_TrTr_TsPr,1), size(window_TrTr_TsPr,2), numChancePerms, size(window_TrTr_TsPr,4));
trial_Window_TrPr_Chance = cell(size(window_TrPr,1), size(window_TrPr,2), numChancePerms, size(window_TrPr,4));
trial_Window_TrTr_TsPo_Chance = cell(size(window_TrTr_TsPo,1), size(window_TrTr_TsPo,2), numChancePerms, size(window_TrTr_TsPo,4));
trial_Window_TrPo_TsTr_Chance = cell(size(window_TrPo_TsTr,1), size(window_TrPo_TsTr,2), numChancePerms, size(window_TrPo_TsTr,4));
trial_Window_TrPo_Chance = cell(size(window_TrPo,1), size(window_TrPo,2), numChancePerms, size(window_TrPo,4));
trial_Window_TrPo_TsPr_Chance = cell(size(window_TrPo_TsPr,1), size(window_TrPo_TsPr,2), numChancePerms, size(window_TrPo_TsPr,4));
trial_Window_TrPr_TsPo_Chance = cell(size(window_TrPr_TsPo,1), size(window_TrPr_TsPo,2), numChancePerms, size(window_TrPr_TsPo,4));
trial_Window_PoPr_Chance = cell(size(window_PoPr,1), size(window_PoPr,2), numChancePerms, size(window_PoPr,4));
trial_TrialPersFit_Chance = cell(size(trlPersFit,1), size(trlPersFit,2), numChancePerms, size(trlPersFit,4));
trial_ITIpersFit_Chance = cell(size(itiPersFit,1), size(itiPersFit,2), numChancePerms, size(itiPersFit,4));
% Chance Beta Correlations
trial_Window_TrPr_TsTr_ChanceBetaCorr = cell(size(window_TrPr_TsTr,1), size(window_TrPr_TsTr,2), numChancePerms, size(window_TrPr_TsTr,4));
trial_Window_TrTr_TsPr_ChanceBetaCorr = cell(size(window_TrTr_TsPr,1), size(window_TrTr_TsPr,2), numChancePerms, size(window_TrTr_TsPr,4));
trial_Window_TrPr_ChanceBetaCorr = cell(size(window_TrPr,1), size(window_TrPr,2), numChancePerms, size(window_TrPr,4));
trial_Window_TrTr_TsPo_ChanceBetaCorr = cell(size(window_TrTr_TsPo,1), size(window_TrTr_TsPo,2), numChancePerms, size(window_TrTr_TsPo,4));
trial_Window_TrPo_TsTr_ChanceBetaCorr = cell(size(window_TrPo_TsTr,1), size(window_TrPo_TsTr,2), numChancePerms, size(window_TrPo_TsTr,4));
trial_Window_TrPo_ChanceBetaCorr= cell(size(window_TrPo,1), size(window_TrPo,2), numChancePerms, size(window_TrPo,4));
trial_Window_TrPo_TsPr_ChanceBetaCorr = cell(size(window_TrPo_TsPr,1), size(window_TrPo_TsPr,2), numChancePerms, size(window_TrPo_TsPr,4));
trial_Window_TrPr_TsPo_ChanceBetaCorr = cell(size(window_TrPr_TsPo,1), size(window_TrPr_TsPo,2), numChancePerms, size(window_TrPr_TsPo,4));
trial_Window_PoPr_ChanceBetaCorr = cell(size(window_PoPr,1), size(window_PoPr,2), numChancePerms, size(window_PoPr,4));
trial_TrialPersFit_ChanceBetaCorr = cell(size(trlPersFit,1), size(trlPersFit,2), numChancePerms, size(trlPersFit,4));
trial_ITIpersFit_ChanceBetaCorr = cell(size(itiPersFit,1), size(itiPersFit,2), numChancePerms, size(itiPersFit,4));
% Chance Theta Correlations
trial_Window_TrPr_TsTr_ChanceThetaCorr = cell(size(window_TrPr_TsTr,1), size(window_TrPr_TsTr,2), numChancePerms, size(window_TrPr_TsTr,4));
trial_Window_TrTr_TsPr_ChanceThetaCorr = cell(size(window_TrTr_TsPr,1), size(window_TrTr_TsPr,2), numChancePerms, size(window_TrTr_TsPr,4));
trial_Window_TrPr_ChanceThetaCorr = cell(size(window_TrPr,1), size(window_TrPr,2), numChancePerms, size(window_TrPr,4));
trial_Window_TrTr_TsPo_ChanceThetaCorr = cell(size(window_TrTr_TsPo,1), size(window_TrTr_TsPo,2), numChancePerms, size(window_TrTr_TsPo,4));
trial_Window_TrPo_TsTr_ChanceThetaCorr = cell(size(window_TrPo_TsTr,1), size(window_TrPo_TsTr,2), numChancePerms, size(window_TrPo_TsTr,4));
trial_Window_TrPo_ChanceThetaCorr= cell(size(window_TrPo,1), size(window_TrPo,2), numChancePerms, size(window_TrPo,4));
trial_Window_TrPo_TsPr_ChanceThetaCorr = cell(size(window_TrPo_TsPr,1), size(window_TrPo_TsPr,2), numChancePerms, size(window_TrPo_TsPr,4));
trial_Window_TrPr_TsPo_ChanceThetaCorr = cell(size(window_TrPr_TsPo,1), size(window_TrPr_TsPo,2), numChancePerms, size(window_TrPr_TsPo,4));
trial_Window_PoPr_ChanceThetaCorr = cell(size(window_PoPr,1), size(window_PoPr,2), numChancePerms, size(window_PoPr,4));
trial_TrialPersFit_ChanceThetaCorr = cell(size(trlPersFit,1), size(trlPersFit,2), numChancePerms, size(trlPersFit,4));
trial_ITIpersFit_ChanceThetaCorr = cell(size(itiPersFit,1), size(itiPersFit,2), numChancePerms, size(itiPersFit,4));

% Start tic
startTic = tic;
permDur = nan(1,numChancePerms);
for perm = 1:numChancePerms
    permTic = tic;
    perm_trlD = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPr_TsTr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrTr_TsPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrTr_TsPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPo_TsTr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPo_TsPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_TrPr_TsPo = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_window_PoPr = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_trlPersFit = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_itiPersFit = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), 2);
    perm_betaPower = cell(mlb{ani}.seqLength, 2, length(fileDirs));
    perm_thetaPower = cell(mlb{ani}.seqLength, 2, length(fileDirs));
    perm_trlMeanBeta = cell(mlb{ani}.seqLength, 2, length(fileDirs));
    perm_trlMeanTheta = cell(mlb{ani}.seqLength, 2, length(fileDirs));
    
    for ani = 1:length(fileDirs)
        %% Process Neural Data
        for al = 1:length(alignments)
            mlb{ani}.alignments = alignments(al);
            mlb{ani}.windows = trlWindows(al);
            % Set the likelihoods
            mlb{ani}.SetLikes_ISC;
            % Shuffle the trials around to destroy position information for the likelihoods
            tempTrlIDmtx = permute(cellfun(@(a)a(1,end),mlb{ani}.likeTrlIDs), [1,3,2]);
            chancePerm = sortrows([randperm(numel(tempTrlIDmtx(~isnan(tempTrlIDmtx))))',find(~isnan(tempTrlIDmtx))]);
            mlb{ani}.likeTrlSpikes(~isnan(tempTrlIDmtx)) = mlb{ani}.likeTrlSpikes(chancePerm(:,2));
            % Process the trials using L1O
            mlb{ani}.Process_IterativeLikelyL1O;
            % Pull LFP 
            [~, tempBetaPower] = mlb{ani}.PP_TrialMatrix_LFP([16 32], trlWindows{al}, alignments{al});
            [~, tempThetaPower] = mlb{ani}.PP_TrialMatrix_LFP([4 12], trlWindows{al}, alignments{al});
            trlIDs = mlb{ani}.postTrlIDs(:);
            trlIDlog = ~isnan(trlIDs);
            posts = permute(mlb{ani}.post(:),[3,2,1]);
            posts = posts(:,:,trlIDlog);
            trlIDs = trlIDs(trlIDlog);
            trlPosVect = [mlb{ani}.trialInfo(trlIDs).Position];
            % Pull out replacement values to replace instances where mean posterior is P==0 or P==1
            unq = unique(cell2mat(posts));
            if unq(1)==0
                minProb = unq(2);
            else
                minProb = unq(1);
            end
            if unq(end)==1
                maxProb = unq(end-1);
            else
                maxProb = unq(end);
            end
            
            for pos = 1:mlb{ani}.seqLength
                % Extract out behavioral events for windowing purposes
                curTrlNdx = mlb{ani}.postTrlIDs(chancePerm(trlPosVect==pos,2));
                perm_betaPower{pos,al,ani} = tempBetaPower(:,:,curTrlNdx);
                perm_thetaPower{pos,al,ani} = tempThetaPower(:,:,curTrlNdx);
                curTrlPokeDurs = [mlb{ani}.trialInfo(curTrlNdx).PokeDuration];
                if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                    nearestPIval = zeros(size(curTrlPokeDurs));
                    nearestPOval = cell2mat(arrayfun(@(a){mlb{ani}.obsvTimeVect(find(a>mlb{ani}.obsvTimeVect,1,'last'))}, curTrlPokeDurs));
                elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                    nearestPIval = cell2mat(arrayfun(@(a){mlb{ani}.obsvTimeVect(find(a>mlb{ani}.obsvTimeVect,1,'last'))}, curTrlPokeDurs))*-1;
                    nearestPOval = zeros(size(curTrlPokeDurs));
                end
                % Now calculate decodability plots
                for p = 1:mlb{ani}.seqLength
                    hits = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect==pos)));
                    fas = cell2mat(cellfun(@(a){a(:,:,p)}, posts(:,:,trlPosVect~=pos)));
                    meanFAR = mean(fas,3);
                    meanFAR(meanFAR==0) = minProb;
                    meanFAR(meanFAR==1) = maxProb;
                    tempTrlD = nan(size(hits));                                     % Trial D matrix
                    tempWindow_TrPr_TsTr = nan(size(hits,3),1);                     % Train Pre-Trial & Test Trial Window
                    tempWindow_TrTr_TsPr = nan(size(hits,3),1);                     % Train Trial & Test Pre-Trial Window
                    tempWindow_TrPr = nan(size(hits,3),1);                          % Pre-Trial vs Test Windows
                    tempWindow_TrTr_TsPo = nan(size(hits,3),1);                     % Train Trial & Test Post-Trial Window
                    tempWindow_TrPo_TsTr = nan(size(hits,3),1);                     % Train Post-Trial & Test Trial Window
                    tempWindow_TrPo = nan(size(hits,3),1);                          % Trial vs Post-Trial Windows
                    tempWindow_TrPo_TsPr = nan(size(hits,3),1);                     % Train Post-Trial & Test Pre-Trial Window
                    tempWindow_TrPr_TsPo = nan(size(hits,3),1);                     % Train Pre-Trial & Test Post-Trial Window
                    tempWindow_PoPr = nan(size(hits,3),1);                          % Post-Trial vs Pre-Trial Windows
                    tempTrialFit = nan(size(hits,3), length(mlb{ani}.obsvTimeVect));     % Persistence model fit Trial period
                    tempITIfit = nan(size(hits,3), length(mlb{ani}.obsvTimeVect));       % Persistence model fit ITI period
                    if p == pos
                        tempTrlBeta = nan(size(hits,3),1);
                        tempTrlTheta = nan(size(hits,3),1);
                    end
                    for trl = 1:size(hits,3)
                        %% Calculate D Matrix
                        tempHR = hits(:,:,trl);
                        tempHR(tempHR==0) = minProb;
                        tempHR(tempHR==1) = maxProb;
                        curTrlD = arrayfun(@(a,b)norminv(a)-norminv(b), tempHR, meanFAR);
                        tempTrlD(:,:,trl) = curTrlD;
                        
                        % Define logical vectors for trial periods
                        preTrialLog = mlb{ani}.obsvTimeVect<nearestPIval(trl);
                        trialLog = mlb{ani}.obsvTimeVect>=nearestPIval(trl) & mlb{ani}.obsvTimeVect<nearestPOval(trl);
                        pstTrialLog = mlb{ani}.obsvTimeVect>=nearestPOval(trl);
                        if strcmp(mlb{ani}.alignments{1}, 'PokeIn')
                            itiLog = preTrialLog;
                        elseif strcmp(mlb{ani}.alignments{1}, 'PokeOut')
                            itiLog = pstTrialLog;
                        end
                        
                        %% Quantify mean windows
                        % IMPORTANT NOTE: The organization of the posterior data is dim1/y-axis/rows = testing time, dim2/x-axis/columns = training time
                        % The organization of most presentations is reversed, i.e. train=rows; test=columns;
                        blankMask = false(size(curTrlD));
                        % Train Pre-Trial & Test Trial
                        prtrWindow = blankMask;
                        prtrWindow(trialLog,preTrialLog) = true;
                        tempWindow_TrPr_TsTr(trl) = mean(curTrlD(prtrWindow));
                        % Train Trial & Test Pre-Trial
                        trprWindow = blankMask;
                        trprWindow(preTrialLog,trialLog) = true;
                        tempWindow_TrTr_TsPr(trl) = mean(curTrlD(trprWindow));
                        % Pre-Trial vs Trial
                        trprprtrWindow = prtrWindow | trprWindow;
                        tempWindow_TrPr(trl) = mean(curTrlD(trprprtrWindow));
                        % Train Trial & Test Post-Trial
                        trpoWindow = blankMask;
                        trpoWindow(pstTrialLog,trialLog) = true;
                        tempWindow_TrTr_TsPo(trl) = mean(curTrlD(trpoWindow));
                        % Train Post-Trial & Test Trial
                        potrWindow = blankMask;
                        potrWindow(trialLog,pstTrialLog) = true;
                        tempWindow_TrPo_TsTr(trl) = mean(curTrlD(potrWindow));
                        % Trial vs Post-Trial
                        trpopotrWindow = trpoWindow | potrWindow;
                        tempWindow_TrPo(trl) = mean(curTrlD(trpopotrWindow));
                        % Train Post-Trial & Test Pre-Trial
                        poprWindow = blankMask;
                        poprWindow(preTrialLog, pstTrialLog) = true;
                        tempWindow_TrPo_TsPr(trl) = mean(curTrlD(poprWindow));
                        % Train Pre-Trial & Test Post-Trial
                        prpoWindow = blankMask;
                        prpoWindow(pstTrialLog, preTrialLog) = true;
                        tempWindow_TrPr_TsPo(trl) = mean(curTrlD(prpoWindow));
                        % Post-Trial vs Pre-Trial
                        poprprpoWindow = poprWindow | prpoWindow;
                        tempWindow_PoPr(trl) = mean(curTrlD(poprprpoWindow));
                        
                        %% Fit the Trial & ITI periods to dynamic coding "models"
                        trialDynCap = sum(trialLog)-1;
                        trialPeriod = zscore(curTrlD(trialLog,trialLog), 0, 'all');
                        for d = 1:trialDynCap-1
                            tempDyn = zscore(triu(true(sum(trialLog)), d*-1) & tril(true(sum(trialLog)), d), 0, 'all');
                            tempTrialFit(trl,d) = pdist([trialPeriod(:)';tempDyn(:)'], 'cosine');
                        end
                        itiDynCap = sum(itiLog)-1;
                        itiPeriod = zscore(curTrlD(itiLog,itiLog), 0 , 'all');
                        for d = 1:itiDynCap-1
                            tempDyn = zscore(triu(true(sum(itiLog)), d*-1) & tril(true(sum(itiLog)), d), 0, 'all');
                            tempITIfit(trl,d) = pdist([itiPeriod(:)';tempDyn(:)'], 'cosine');
                        end
                        
                        %% Calculate average LFP power (if applicable)
                        if p == pos
                            tempTrlBeta(trl) = mean(perm_betaPower{pos,al,ani}(trialLog,1,trl));
                            tempTrlTheta(trl) = mean(perm_thetaPower{pos,al,ani}(trialLog,1,trl));
                        end
                    end
                    perm_trlD{pos,p,ani,al} = tempTrlD;
                    perm_window_TrPr_TsTr{pos,p,ani,al} = tempWindow_TrPr_TsTr;
                    perm_window_TrTr_TsPr{pos,p,ani,al} = tempWindow_TrTr_TsPr;
                    perm_window_TrPr{pos,p,ani,al} = tempWindow_TrPr;
                    perm_window_TrTr_TsPo{pos,p,ani,al} = tempWindow_TrTr_TsPo;
                    perm_window_TrPo_TsTr{pos,p,ani,al} = tempWindow_TrPo_TsTr;
                    perm_window_TrPo{pos,p,ani,al} = tempWindow_TrPo;
                    perm_window_TrPo_TsPr{pos,p,ani,al} = tempWindow_TrPo_TsPr;
                    perm_window_TrPr_TsPo{pos,p,ani,al} = tempWindow_TrPr_TsPo;
                    perm_window_PoPr{pos,p,ani,al} = tempWindow_PoPr;
                    perm_trlPersFit{pos,p,ani,al} = tempTrialFit;
                    perm_itiPersFit{pos,p,ani,al} = tempITIfit;
                    if p==pos
                        perm_trlMeanBeta{pos,al,ani} = tempTrlBeta;
                        perm_trlMeanTheta{pos,al,ani} = tempTrlTheta;
                    end
                end
            end
        end
    end
    %% Collapse trial data across animals
    for r = 1:mlb{1}.seqLength
        tempMeanBeta_Chance = cell(1,length(alignments));
        tempMeanTheta_Chance = cell(1,length(alignments));
        for al = 1:length(alignments)
            tempMeanBeta_Chance{al} = cell2mat(permute(perm_trlMeanBeta(r,al,:), [3,1,2]));
            tempMeanTheta_Chance{al} = cell2mat(permute(perm_trlMeanTheta(r,al,:), [3,1,2]));
        end
        for c = 1:mlb{1}.seqLength
            for al = 1:length(alignments)
                trialD_Chance{r,c,perm,al} = mean(cell2mat(perm_trlD(r,c,:,al)),3);
                tempPerm_TrPr_TsTr = cell2mat(permute(perm_window_TrPr_TsTr(r,c,:,al), [3,1,2]));
                trial_Window_TrPr_TsTr_Chance{r,c,perm,al} = mean(tempPerm_TrPr_TsTr);
                trial_Window_TrPr_TsTr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPr_TsTr, tempMeanBeta_Chance{al});
                trial_Window_TrPr_TsTr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPr_TsTr, tempMeanTheta_Chance{al});
                
                tempPerm_TrTr_TsPr = cell2mat(permute(perm_window_TrTr_TsPr(r,c,:,al), [3,1,2]));
                trial_Window_TrTr_TsPr_Chance{r,c,perm,al} = mean(tempPerm_TrTr_TsPr);
                trial_Window_TrTr_TsPr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrTr_TsPr, tempMeanBeta_Chance{al});
                trial_Window_TrTr_TsPr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrTr_TsPr, tempMeanTheta_Chance{al});
                
                tempPerm_TrPr = cell2mat(permute(perm_window_TrPr(r,c,:,al), [3,1,2]));
                trial_Window_TrPr_Chance{r,c,perm,al} = mean(tempPerm_TrPr);
                trial_Window_TrPr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPr, tempMeanBeta_Chance{al});
                trial_Window_TrPr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPr, tempMeanTheta_Chance{al});
                
                tempPerm_TrTr_TsPo = cell2mat(permute(perm_window_TrTr_TsPo(r,c,:,al), [3,1,2]));
                trial_Window_TrTr_TsPo_Chance{r,c,perm,al} = mean(tempPerm_TrTr_TsPo);
                trial_Window_TrTr_TsPo_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrTr_TsPo, tempMeanBeta_Chance{al});
                trial_Window_TrTr_TsPo_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrTr_TsPo, tempMeanTheta_Chance{al});
                
                tempPerm_TrPo_TsTr = cell2mat(permute(perm_window_TrPo_TsTr(r,c,:,al), [3,1,2]));
                trial_Window_TrPo_TsTr_Chance{r,c,perm,al} = mean(tempPerm_TrPo_TsTr);
                trial_Window_TrPo_TsTr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPo_TsTr, tempMeanBeta_Chance{al});
                trial_Window_TrPo_TsTr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPo_TsTr, tempMeanTheta_Chance{al});
                
                tempPerm_TrPo = cell2mat(permute(perm_window_TrPo(r,c,:,al), [3,1,2]));
                trial_Window_TrPo_Chance{r,c,perm,al} = mean(tempPerm_TrPo);
                trial_Window_TrPo_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPo, tempMeanBeta_Chance{al});
                trial_Window_TrPo_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPo, tempMeanTheta_Chance{al});
                
                tempPerm_TrPo_TsPr = cell2mat(permute(perm_window_TrPo_TsPr(r,c,:,al), [3,1,2]));
                trial_Window_TrPo_TsPr_Chance{r,c,perm,al} = mean(tempPerm_TrPo_TsPr);
                trial_Window_TrPo_TsPr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPo_TsPr, tempMeanBeta_Chance{al});
                trial_Window_TrPo_TsPr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPo_TsPr, tempMeanTheta_Chance{al});
                                
                tempPerm_TrPr_TsPo = cell2mat(permute(perm_window_TrPr_TsPo(r,c,:,al), [3,1,2]));
                trial_Window_TrPr_TsPo_Chance{r,c,perm,al} = mean(tempPerm_TrPr_TsPo);
                trial_Window_TrPr_TsPo_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_TrPr_TsPo, tempMeanBeta_Chance{al});
                trial_Window_TrPr_TsPo_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_TrPr_TsPo, tempMeanTheta_Chance{al});
                
                tempPerm_PoPr = cell2mat(permute(perm_window_PoPr(r,c,:,al), [3,1,2]));
                trial_Window_PoPr_Chance{r,c,perm,al} = mean(tempPerm_PoPr);
                trial_Window_PoPr_ChanceBetaCorr{r,c,perm,al} = corr(tempPerm_PoPr, tempMeanBeta_Chance{al});
                trial_Window_PoPr_ChanceThetaCorr{r,c,perm,al} = corr(tempPerm_PoPr, tempMeanTheta_Chance{al});
                
                trial_TrialPersFit_Chance{r,c,perm,al} = cell2mat(permute(perm_trlPersFit(r,c,:,al), [3,1,2]));
                trial_ITIpersFit_Chance{r,c,perm,al} = cell2mat(permute(perm_itiPersFit(r,c,:,al), [3,1,2]));
            end
        end
    end
    clc;
    permDur(perm) = toc(permTic)/60;
    fprintf('Perm %i completed after %.02fmin, estimated %.02fhrs remaining\n', perm, permDur(perm), ((toc(startTic)/perm)*(numChancePerms-perm))/3600);
end
%% Colormap setup
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);

%% Plot things & do stats... you know, the sciency visual stuff!
% Create logicals for trial categories
isLog = logical(eye(mlb{1}.seqLength));
nxtPosLog = triu(ones(mlb{1}.seqLength),-1) & tril(ones(mlb{1}.seqLength), -1);
prvPosLog = triu(ones(mlb{1}.seqLength),1) & tril(ones(mlb{1}.seqLength), 1);

grpPiPoLat = median(cell2mat(piPokeOutLat));
grpPiRwdLat = median(cell2mat(piRwdLat));
grpPoPiLat = median(cell2mat(poPokeInLat));
grpPoRwdLat = median(cell2mat(poRwdLat));
% grpNxtTrlLat = median(cell2mat(nxtTrlLat));
%% Plot Trial D: All positions & decodings
for al = 1:length(alignments)
    figure;
    for r = 1:mlb{end}.seqLength
        for c = 1:mlb{end}.seqLength
            subplot(mlb{end}.seqLength, mlb{end}.seqLength, sub2ind([mlb{end}.seqLength, mlb{end}.seqLength], c,r))
            mlb{end}.PlotTrialPDM(trialD{r,c,al}, 'rotate', 'clim', [-1 1], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%             mlb{end}.PlotTrialPDM(trialDz{r,c,al}, 'rotate', 'clim', [-0.5 0.5], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%             mlb{end}.PlotTrialPDM(trialDposZ{r,c,al}, 'rotate', 'clim', [-1 1], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%             mlb{end}.PlotTrialPDM(trialHR{r,c,al}, 'rotate', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            title(sprintf('Obsv:%i, Decode:%i', r,c));
            hold on;
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);

                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            end
        end
    end    
    colormap(cMap)
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    linkaxes;
end
%% Plot Trial D: Select Prev/Next Trial Information
for al = 1:length(alignments)
    figure; 
    for pos = 1:mlb{end}.seqLength
        if pos>1
            subplot(mlb{end}.seqLength, 3, sub2ind([3, mlb{end}.seqLength], 1,pos));
            mlb{end}.PlotTrialPDM(trialD{pos-1,pos,al}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            title('Upcoming Position Decoding');
            hold on;
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);

                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            end
        end
        subplot(mlb{end}.seqLength, 3, sub2ind([3, mlb{end}.seqLength], 2,pos));
        mlb{end}.PlotTrialPDM(trialD{pos,pos,al}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
        title('Current Trial');
            hold on;
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);

                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            end
        if pos<mlb{end}.seqLength
            subplot(mlb{end}.seqLength, 3, sub2ind([3, mlb{end}.seqLength], 3,pos));
            mlb{end}.PlotTrialPDM(trialD{pos+1,pos,al}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            title('Previous Trial Decoding');
            hold on;
            if strcmp(alignments{al} , 'PokeIn')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);

                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            elseif strcmp(alignments{al}, 'PokeOut')
                plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
                plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
                
                plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
                plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
            end
        end
    end
    colormap(cMap)
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% Plot Trial D: Mean across positions
for al = 1:length(alignments)
    figure;
    tempTMat = trialD(:,:,al);
    subplot(1,3,1)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(mlb{1}.seqLength),1) & tril(ones(mlb{1}.seqLength),1)), [2,3,1])), 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    subplot(1,3,2)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(mlb{1}.seqLength))), [2,3,1])), 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    subplot(1,3,3)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(mlb{1}.seqLength),-1) & tril(ones(mlb{1}.seqLength),-1)), [2,3,1])), 'rotate', 'clim',  [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    colormap(cMap)
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% Plot Trial D: Alignment Specific Means
for al = 1:length(alignments)
    figure;
    tempTMat = trialD(:,:,al);
    if strcmp(alignments{al}, 'PokeIn')
        tempTMat = tempTMat(2:end,2:end);        
%         tempTMat = tempTMat(2:end-1,2:end-1);
    elseif strcmp(alignments{al}, 'PokeOut')
%                 tempTMat = tempTMat(2:end,2:end);        
%         tempTMat = tempTMat(1:end-1,1:end-1);
%         tempTMat = tempTMat(2:end-1,2:end-1);
    end
    subplot(1,3,1)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1)), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    subplot(1,3,2)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    set(gca,'clim', [-2 2]);
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    subplot(1,3,3)
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    set(gca,'clim', [-2 2]);
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    colormap(cMap)
    set(gca,'clim', [-2 2]);
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% Plot Specific Window Comparisons
% selType = 0; % Alignment includes all trial positions
selType = 1; % Alignment specific trial positions included
% selType = 2; % Alignment only common trial positions i.e. 2&3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Upcoming Trial Info
figure;
subplot(4,3,1)
tempTMat = trialD(:,:,1);
if selType == 1
    tempTMat = tempTMat(2:end,2:end);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [0, 0, grpPiPoLat, grpPiPoLat],...
    'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
    'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
    'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
title('Previous Position Info');

subplot(4,3,2)
tempTMat = trialD(:,:,2);
if selType == 1
    tempTMat = tempTMat(1:end-1,1:end-1);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [grpPoPiLat, grpPoPiLat, 0, 0],...
    'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);

    
subplot(4,3,3)
tempTrPr = trial_Window_TrPr(:,:,1);
tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
tempTrPo = trial_Window_TrPo(:,:,2);
tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
if selType == 1
    tempTrPr = tempTrPr(2:end,2:end);
    tempPrTr_Window = tempPrTr_Window(2:end,2:end);
    tempTrPr_Window = tempTrPr_Window(2:end,2:end);
    tempTrPo = tempTrPo(1:end-1,1:end-1);
    tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
    tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
elseif selType == 2
    tempTrPr = tempTrPr(2:end-1,2:end-1);
    tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
    tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
    tempTrPo = tempTrPo(2:end-1,2:end-1);
    tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
    tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
end
trPr = cell2mat(tempTrPr(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
win_PrTr = cell2mat(tempPrTr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
win_TrPr = cell2mat(tempTrPr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
trPo = cell2mat(tempTrPo(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
win_PoTr = cell2mat(tempPoTr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
win_TrPo = cell2mat(tempTrPo_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
hold on;
swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Current Trial
subplot(4,3,4)
tempTMat = trialD(:,:,1);
if selType == 1
    tempTMat = tempTMat(2:end,2:end);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [0, 0, grpPiPoLat, grpPiPoLat], 'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0],...
    'linestyle', '-', 'linewidth', 5);
patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
    'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
    'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
title('Current Position Info');

subplot(4,3,5)
tempTMat = trialD(:,:,2);
if selType == 1
    tempTMat = tempTMat(1:end-1,1:end-1);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [grpPoPiLat, grpPoPiLat, 0, 0], 'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255],...
    'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);

subplot(4,3,6)
tempTrPr = trial_Window_TrPr(:,:,1);
tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
tempTrPo = trial_Window_TrPo(:,:,2);
tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
if selType == 1
    tempTrPr = tempTrPr(2:end,2:end);
    tempPrTr_Window = tempPrTr_Window(2:end,2:end);
    tempTrPr_Window = tempTrPr_Window(2:end,2:end);
    tempTrPo = tempTrPo(1:end-1,1:end-1);
    tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
    tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
elseif selType == 2
    tempTrPr = tempTrPr(2:end-1,2:end-1);
    tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
    tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
    tempTrPo = tempTrPo(2:end-1,2:end-1);
    tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
    tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
end
trPr = cell2mat(tempTrPr(logical(eye(length(tempTrPr)))));
win_PrTr = cell2mat(tempPrTr_Window(logical(eye(length(tempTrPr)))));
win_TrPr = cell2mat(tempTrPr_Window(logical(eye(length(tempTrPr)))));
trPo = cell2mat(tempTrPo(logical(eye(length(tempTrPr)))));
win_TrPo = cell2mat(tempTrPo_Window(logical(eye(length(tempTrPr)))));
win_PoTr = cell2mat(tempPoTr_Window(logical(eye(length(tempTrPr)))));
bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
hold on;
swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_PrTr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);
% 
% [h,p,ci,stats] = ttest(trPr, trPo);
% [p,tbl,stats] = anova1([win_PrTr,win_TrPr,win_PoTr,win_TrPo]);
% 
% figure;
% multcompare(stats, 'CType', 'bonferroni');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Next Trial
subplot(4,3,7)
tempTMat = trialD(:,:,1);
if selType == 1
    tempTMat = tempTMat(2:end,2:end);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [0, 0, grpPiPoLat, grpPiPoLat],...
    'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
    'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
    'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
title('Next Position Info');

subplot(4,3,8)
tempTMat = trialD(:,:,2);
if selType == 1
    tempTMat = tempTMat(1:end-1,1:end-1);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [grpPoPiLat, grpPoPiLat, 0, 0],...
    'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '--', 'linewidth', 5);
patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
    'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
    'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);

subplot(4,3,9)
tempTrPr = trial_Window_TrPr(:,:,1);
tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
tempTrPo = trial_Window_TrPo(:,:,2);
tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
if selType == 1
    tempTrPr = tempTrPr(2:end,2:end);
    tempPrTr_Window = tempPrTr_Window(2:end,2:end);
    tempTrPr_Window = tempTrPr_Window(2:end,2:end);
    tempTrPo = tempTrPo(1:end-1,1:end-1);
    tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
    tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
elseif selType == 2
    tempTrPr = tempTrPr(2:end-1,2:end-1);
    tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
    tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
    tempTrPo = tempTrPo(2:end-1,2:end-1);
    tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
    tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
end

trPr = cell2mat(tempTrPr(triu(ones(length(tempTrPr)),1) & tril(ones(length(tempTrPr)),1)));
win_PrTr = cell2mat(tempPrTr_Window(triu(ones(length(tempPrTr_Window)),1) & tril(ones(length(tempPrTr_Window)),1)));
win_TrPr = cell2mat(tempTrPr_Window(triu(ones(length(tempTrPr_Window)),1) & tril(ones(length(tempTrPr_Window)),1)));
trPo = cell2mat(tempTrPo(triu(ones(length(tempTrPo)),1) & tril(ones(length(tempTrPo)),1)));
win_PoTr = cell2mat(tempPoTr_Window(triu(ones(length(tempPoTr_Window)),1) & tril(ones(length(tempPoTr_Window)),1)));
win_TrPo = cell2mat(tempTrPo_Window(triu(ones(length(tempTrPo_Window)),1) & tril(ones(length(tempTrPo_Window)),1)));
bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
hold on;
swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% ITI Generalization Indices
tempPoPr_Window = trial_Window_TrPo_TsPr(:,:,1);
tempPrPo_Window = trial_Window_TrPr_TsPo(:,:,1);
if selType == 1
    tempPoPr_Window = tempPoPr_Window(2:end,2:end);
    tempPrPo_Window = tempPrPo_Window(2:end,2:end);
elseif selType == 2
    tempPoPr_Window = tempPoPr_Window(2:end-1,2:end-1);
    tempPrPo_Window = tempPrPo_Window(2:end-1,2:end-1);
end
popr_Prev = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),-1) & tril(ones(length(tempPoPr_Window)),-1)));
popr_Curr = cell2mat(tempPoPr_Window(logical(eye(length(tempPoPr_Window)))));
popr_Next = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),1) & tril(ones(length(tempPoPr_Window)),1)));
prpo_Prev = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),-1) & tril(ones(length(tempPrPo_Window)),-1)));
prpo_Curr = cell2mat(tempPrPo_Window(logical(eye(length(tempPrPo_Window)))));
prpo_Next = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),1) & tril(ones(length(tempPrPo_Window)),1)));
subplot(4,3,10)
bar([1,2,3, 5,6,7], [mean(popr_Prev), mean(popr_Curr), mean(popr_Next), mean(prpo_Prev), mean(prpo_Curr), mean(prpo_Next)], 'k');
hold on;
swarmchart(ones(size(popr_Prev)), popr_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(popr_Curr))+1, popr_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(popr_Next))+2, popr_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Prev))+4, prpo_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Curr))+5, prpo_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Next))+6, prpo_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
set(gca, 'xtick', [1:3,5:7], 'xticklabel', [{'Po->Pr:Prev'}, {'Po->Pr:Curr'}, {'Po->Pr:Next'}, {'Pr->Po:Prev'}, {'Pr->Po:Curr'}, {'Pr->Po:Next'}], 'xticklabelrotation', 45);

tempPoPr_Window = trial_Window_TrPo_TsPr(:,:,2);
tempPrPo_Window = trial_Window_TrPr_TsPo(:,:,2);
if selType == 1
    tempPoPr_Window = tempPoPr_Window(1:end-1,1:end-1);
    tempPrPo_Window = tempPrPo_Window(1:end-1,1:end-1);
elseif selType == 2
    tempPoPr_Window = tempPoPr_Window(2:end-1,2:end-1);
    tempPrPo_Window = tempPrPo_Window(2:end-1,2:end-1);
end
popr_Prev = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),-1) & tril(ones(length(tempPoPr_Window)),-1)));
popr_Curr = cell2mat(tempPoPr_Window(logical(eye(length(tempPoPr_Window)))));
popr_Next = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),1) & tril(ones(length(tempPoPr_Window)),1)));
prpo_Prev = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),-1) & tril(ones(length(tempPrPo_Window)),-1)));
prpo_Curr = cell2mat(tempPrPo_Window(logical(eye(length(tempPrPo_Window)))));
prpo_Next = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),1) & tril(ones(length(tempPrPo_Window)),1)));
subplot(4,3,11)
bar([1,2,3, 5,6,7], [mean(popr_Prev), mean(popr_Curr), mean(popr_Next), mean(prpo_Prev), mean(prpo_Curr), mean(prpo_Next)], 'k');
hold on;
swarmchart(ones(size(popr_Prev)), popr_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(popr_Curr))+1, popr_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(popr_Next))+2, popr_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Prev))+4, prpo_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Curr))+5, prpo_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
swarmchart(ones(size(prpo_Next))+6, prpo_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
set(gca, 'xtick', [1:3,5:7], 'xticklabel', [{'Po->Pr:Prev'}, {'Po->Pr:Curr'}, {'Po->Pr:Next'}, {'Pr->Po:Prev'}, {'Pr->Po:Curr'}, {'Pr->Po:Next'}], 'xticklabelrotation', 45);
colormap(cMap)

if selType == 1
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'Select Positions Based on Windows',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
elseif selType == 2
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'Common Positions Across Alignments (pos 2&3)',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
else
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'All Trials',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

tempPoPr_PI = trial_Window_PoPr(:,:,1);
% tempPoPr_PI = tempPoPr_PI(2:end,2:end);
% tempPoPr_PI = tempPoPr_PI(2:end-1,2:end-1);
tempPoPr_PO = trial_Window_PoPr(:,:,2);
% tempPoPr_PO = tempPoPr_PO(1:end-1,1:end-1);
% tempPoPr_PO = tempPoPr_PO(2:end-1,2:end-1);




% trial_Window_TrPo_TsPr
% trial_Window_TrPr_TsPo
% trial_Window_PoPr
%% Plot Fits
% Plot poke in aligned D data
figure;
subplot(2,3,1)
tempTMat = trialD(:,:,1);
if selType == 1
    tempTMat = tempTMat(2:end,2:end);
elseif selType == 2
    tempTMat = tempTMat(2:end-1,2:end-1);
end
mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
hold on;
plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
title('Poke In Aligned');
% Plot poke in aligned fits
tempTrialFits = trial_TrialPersFit(:,:,1);
tempITIfits = trial_ITIpersFit(:,:,1);
if selType == 1
    tempTrialFits = tempTrialFits(2:end,2:end);
    tempITIfits = tempITIfits(2:end,2:end);
elseif selType == 2
    tempTrialFits = tempTrialFits(2:end-1,2:end-1);
    tempITIfits = tempITIfits(2:end-1,2:end-1);
end
tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
[minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
minTrialLat = minTrialLat./(1/dsRate);
[minITIfit,minITIlat] = min(tempITIfits,[],2);
minITIlat = minITIlat./(1/dsRate);
subplot(2,3,2)
swarmchart(ones(size(minTrialLat)), minTrialLat, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
hold on;
swarmchart(ones(size(minITIlat))+1, minITIlat, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
[~,p,ci,stats] = ttest(minTrialLat, minITIlat);
if p<0.05
    title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
else
    title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
end
set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
ylabel('Latency to Best Fit (ms)');
subplot(2,3,3)
swarmchart(ones(size(minTrialFit)), minTrialFit, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
hold on;
swarmchart(ones(size(minITIfit))+1, minITIfit, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
[~,p,~,stats] = ttest(minTrialFit, minITIfit);
if p<0.05
    title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
else
    title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
end
set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
ylabel('Best Fit (Cos Sim)');

% Plot poke out aligned D data
if length(alignments)>1
    subplot(2,3,4)
    tempTMat = trialD(:,:,2);
    if selType == 1
        tempTMat = tempTMat(1:end-1,1:end-1);
    elseif selType == 2
        tempTMat = tempTMat(2:end-1,2:end-1);
    end
    mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
    plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
    plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
    plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    title('Poke Out Aligned');
    % Plot poke out aligned fits
    tempTrialFits = trial_TrialPersFit(:,:,2);
    tempITIfits = trial_ITIpersFit(:,:,2);
    if selType == 1
        tempTrialFits = tempTrialFits(1:end-1,1:end-1);
        tempITIfits = tempITIfits(1:end-1,1:end-1);
    elseif selType == 2
        tempTrialFits = tempTrialFits(2:end-1,2:end-1);
        tempITIfits = tempITIfits(2:end-1,2:end-1);
    end
    tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
    tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
    [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
    minTrialLat = minTrialLat./(1/dsRate);
    [minITIfit,minITIlat] = min(tempITIfits,[],2);
    minITIlat = minITIlat./(1/dsRate);
    subplot(2,3,5)
    swarmchart(ones(size(minTrialLat)), minTrialLat, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
    hold on;
    swarmchart(ones(size(minITIlat))+1, minITIlat, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
    [~,p,~,stats] = ttest(minTrialLat, minITIlat);
    if p<0.05
        title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
    else
        title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
    end
    set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
    ylabel('Latency to Best Fit (ms)');
    subplot(2,3,6)
    swarmchart(ones(size(minTrialFit)), minTrialFit, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
    hold on;
    swarmchart(ones(size(minITIfit))+1, minITIfit, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
    [~,p,~,stats] = ttest(minTrialFit, minITIfit);
    if p<0.05
        title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
    else
        title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
    end
    set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
    ylabel('Best Fit (Cos Sim)');
end 
    colormap(cMap)


if selType == 1
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'Select Positions Based on Windows',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
elseif selType == 2
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'Common Positions Across Alignments (pos 2&3)',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
else
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'All Trials',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% Holy shit we have a beta effect!
for al = 1:length(alignments)
    figure;    
    betaLog = true(1,4);
    tempLog = logical(eye(mlb{ani}.seqLength));
    decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
    if selType == 1
%         if al == 1
%             tempLog(4,4) = false;
%             betaLog(4) = false;
%         elseif al == 2
            tempLog(1,1) = false;
            betaLog(1) = false;
%         end
    elseif selType == 2        
        tempLog(1,1) = false;
        tempLog(4,4) = false;
        betaLog(1) = false;
        betaLog(4) = false;
    elseif selType == 3
        tempLog((1:mlb{ani}.seqLength)~=pos,(1:mlb{ani}.seqLength)~=pos) = false;
        betaLog((1:mlb{ani}.seqLength)~=pos) = false;
    end
        
    decLog(:,:,al) = tempLog;
    subplot(4,3,1);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr_TsTr(decLog)), 'Beta', 'PrTr', []);
    subplot(4,3,2);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrTr_TsPr(decLog)), 'Beta', 'TrPr', []);
    subplot(4,3,3);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'PrTrTrPr', []);
    subplot(4,3,4);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrTr_TsPo(decLog)), 'Beta', 'TrPo', []);
    subplot(4,3,5);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo_TsTr(decLog)), 'Beta', 'PoTr', []);
    subplot(4,3,6);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'TrPoPoTr', []);
    subplot(4,3,7);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo_TsPr(decLog)), 'Beta', 'PoPr', []);
    subplot(4,3,8);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr_TsPo(decLog)), 'Beta', 'PrPo', []);
    subplot(4,3,9);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_PoPr(decLog)), 'Beta', 'PoPrPrPo', []);
    
    tempTrialFits = cell2mat(trial_TrialPersFit(decLog));
    [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
    minTrialLat = minTrialLat./(1/dsRate);
    tempITIfits = cell2mat(trial_ITIpersFit(decLog));
    [minITIfit,minITIlat] = min(tempITIfits,[],2);
    minITIlat = minITIlat./(1/dsRate);
    subplot(4,3,10);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), minTrialLat, 'Beta', 'Trial Pers', []);
    subplot(4,3,11);
    corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), minITIlat, 'Beta', 'ITI Pers', []);
    
    if selType == 1
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf('%s: Select Positions Based on Windows', alignments{al}),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    elseif selType == 2
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf('%s: Common Positions Across Alignments (pos 2&3)', alignments{al}),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    else
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf('%s: All Trials', alignments{al}),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    end
end
% trial_TrialPersFit = cell(size(trlPersFit,1), size(trlPersFit,2), size(trlPersFit,4));
% trial_ITIpersFit = cell(size(itiPersFit,1), size(itiPersFit,2), size(itiPersFit,4));

%% Pre-Trial vs Trial Port Entry Aligned
figure; 
tempLog = logical(eye(mlb{ani}.seqLength));
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(:,:,1) = tempLog;
subplot(4,4,[1:3,5:7,9:11,13:15])
corrScatPlot(cell2mat(trial_MeanBeta(:,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('All Positions');

subplot(4,4,4);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(1,1,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 1');

subplot(4,4,8);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(2,2,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 2');

subplot(4,4,12);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(3,3,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 3');

subplot(4,4,16);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(4,4,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 4');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Port Entry Aligned: Pre-Trial vs. Trial; Current Trial Info',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Post-Trial vs Trial Port Entry Aligned
figure; 
tempLog = logical(eye(mlb{ani}.seqLength));
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(:,:,2) = tempLog;
subplot(4,4,[1:3,5:7,9:11,13:15])
corrScatPlot(cell2mat(trial_MeanBeta(:,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('All Positions');

subplot(4,4,4);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(1,1,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 1');

subplot(4,4,8);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(2,2,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 2');

subplot(4,4,12);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(3,3,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 3');

subplot(4,4,16);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(4,4,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Position 4');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Port Exit Aligned: Post-Trial vs. Trial; Current Trial Info',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');


%% Pre-Trial vs Trial Port Entry Aligned Previous Trial Info
figure; 
tempLog = triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(:,:,1) = tempLog;
subplot(3,3,[1:2,4:5,7:8])
corrScatPlot(cell2mat(trial_MeanBeta(2:end,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('All Positions');

subplot(3,3,3);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(2,1,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 1 in Pos 2');

subplot(3,3,6);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(3,2,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 2 in Pos 3');

subplot(3,3,9);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(4,3,1) = true;
corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 3 in Pos 4');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Port Entry Aligned: Pre-Trial vs. Trial; Prev Trial Info',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Pre-Trial vs Trial Port Entry Aligned Previous Trial Info
figure; 
tempLog = triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(:,:,2) = tempLog;
subplot(3,3,[1:2,4:5,7:8])
corrScatPlot(cell2mat(trial_MeanBeta(1:end-1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('All Positions');

subplot(3,3,3);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(1,2,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 1 in Pos 2');

subplot(3,3,6);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(2,3,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 2 in Pos 3');

subplot(3,3,9);
decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
decLog(3,4,2) = true;
corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
title('Pos 3 in Pos 4');

annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
    'String', 'Port Entry Aligned: Post-Trial vs. Trial; Next Trial Info',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% NEW FIGURES 08/06/22
% %%%%%% Variables %%%%%% % 
% NOTE: Visuals depict plots with training time as rows and testing time as columns. In the actual data, everything's set up to depict
% testing time as the rows and training time as the columns.
% %%%%%%%%%%%%%  Cross-Epoch Variables %%%%%%%%%%%%% %
%       --------------
%       |    |    |  |
%       |____|____|__|
%       |%%%%|    |  |
%       |%%%%|    |  |
%       |%%%%|____|__|
%       |    |%%%%|  |
%       |    |%%%%|  |
%       |    |%%%%|  |
%       --------------
% trial_Window_TrPr_TsTr
% trial_Window_TrTr_TsPr
% trial_Window_TrPr
%       --------------
%       |  |%%%%|    |
%       |  |%%%%|    |
%       |__|%%%%|____|
%       |  |    |%%%%|
%       |  |    |%%%%|
%       |__|____|%%%%|
%       |  |    |    |
%       |  |    |    |
%       --------------
% trial_Window_TrTr_TsPo
% trial_Window_TrPo_TsTr
% trial_Window_TrPo 

% %%%%%%%%%%%%%  Within-Epoch Variables %%%%%%%%%%%%% %
%       --------------
%       |    |    |  |
%       |____|____|__|
%       |    |%%%%|  |
%       |    |%%%%|  |
%       |    |%%%%|  |
%       |____|%%%%|__|
%       |%%%%|    |  |
%       |%%%%|    |  |
%       |%%%%|    |  |
%       |%%%%|    |  |
%       --------------
% trial_Window_TrTr_TsTr 
% trial_Window_TrITD_TsITD 
%       --------------
%       |    |    |  |
%       |____|____|__|
%       |    |   %|  |
%       |    |  % |  |
%       |    | %  |  |
%       |____|%___|__|
%       |   %|    |  |
%       |  % |    |  |
%       | %  |    |  |
%       |%   |    |  |
%       --------------
% trial_Window_TrTr_TsTr_Diag
% trial_Window_TrITD_TsITD_Diag 
%       --------------
%       |    |    |  |
%       |____|____|__|
%       |    |%%% |  |
%       |    |%% %|  |
%       |    |% %%|  |
%       |____|_%%%|__|
%       |%%% |    |  |
%       |%% %|    |  |
%       |% %%|    |  |
%       | %%%|    |  |
%       --------------
% trial_Window_TrTr_TsTr_OffDiag 
% trial_Window_TrITD_TsITD_OffDiag 

% %%%%%%%%%%%%%  Persistence Fit Values %%%%%%%%%%%%% %
% trial_TrialPersFit 
% trial_ITIpersFit 

% %%%%%%%%%%%%%  LFP Variables %%%%%%%%%%%%% %
% trial_BetaPowerTime 
% trial_ThetaPowerTime 
% trial_MeanBeta 
% trial_MeanTheta 

%% New Fig 1: Within vs across epochs per alignment
% Port Entry Aligned
for al = 1:length(alignments)
    piWIn_Trl = trial_Window_TrTr_TsTr(:,:,al);
    piWIn_ITD = trial_Window_TrITD_TsITD(:,:,al);
    if strcmp(alignments{al}, 'PokeIn')
        piX_PrTr = trial_Window_TrPr_TsTr(:,:,al);
        piX_TrPr = trial_Window_TrTr_TsPr(:,:,al);
        piX_TrPo = trial_Window_TrTr_TsPo(:,:,al);
        piX_PoTr = trial_Window_TrPo_TsTr(:,:,al);
    elseif strcmp(alignments{al}, 'PokeOut')
        piX_PrTr = trial_Window_TrPo_TsTr(:,:,al);
        piX_TrPr = trial_Window_TrTr_TsPo(:,:,al);
        piX_TrPo = trial_Window_TrTr_TsPr(:,:,al);
        piX_PoTr = trial_Window_TrPr_TsTr(:,:,al);
    end
    piWIn_Trl = piWIn_Trl(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    piWIn_ITD = piWIn_ITD(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    piX_PrTr = piX_PrTr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    piX_TrPr = piX_TrPr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    piX_TrPo = piX_TrPo(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    piX_PoTr = piX_PoTr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    
    spIDs = reshape(1:mlb{ani}.seqLength*(mlb{ani}.seqLength+1), [mlb{ani}.seqLength+1, mlb{ani}.seqLength])';
    figure;
    subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,reshape(spIDs(:,1:end-1), [1,mlb{ani}.seqLength^2]));
    mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(piWIn_Trl),1,0.05,'r');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(piWIn_ITD),1,0.05,'b');
    mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(piX_PrTr),1,0.05,[191/255, 191/255, 0]);
    mlb{ani}.PlotMeanVarSwarmBar(4,cell2mat(piX_PoTr),1,0.05,[191/255, 191/255, 191/255]);
    mlb{ani}.PlotMeanVarSwarmBar(5,cell2mat(piX_TrPr),1,0.05,[191/255, 191/255, 0]);
    mlb{ani}.PlotMeanVarSwarmBar(6,cell2mat(piX_TrPo),1,0.05,[191/255, 191/255, 191/255]);
    if strcmp(alignments{al}, 'PokeIn')
        set(gca, 'xtick', 1:6, 'xticklabel', [{'Trial'}, {'ITD'}, {'Pre->Trial'}, {'Post->Trial'}, {'Trial->Pre'}, {'Trial->Post'}], 'xticklabelrotation', 45);
    elseif strcmp(alignments{al}, 'PokeOut')
        set(gca, 'xtick', 1:6, 'xticklabel', [{'Trial'}, {'ITD'}, {'Post->Trial'}, {'Pre->Trial'}, {'Trial->Post'}, {'Trial->Pre'}], 'xticklabelrotation', 45);
    end
    for sp = 1:mlb{ani}.seqLength
        subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,spIDs(sp,end));
        mlb{ani}.PlotMeanVarSwarmBar(1,piWIn_Trl{sp},1,0.05,'r');
        hold on;
        mlb{ani}.PlotMeanVarSwarmBar(2,piWIn_ITD{sp},1,0.05,'b');
        mlb{ani}.PlotMeanVarSwarmBar(3,piX_PrTr{sp},1,0.05,[191/255, 191/255, 0]);
        mlb{ani}.PlotMeanVarSwarmBar(4,piX_PoTr{sp},1,0.05,[191/255, 191/255, 191/255]);
        mlb{ani}.PlotMeanVarSwarmBar(5,piX_TrPr{sp},1,0.05,[191/255, 191/255, 0]);
        mlb{ani}.PlotMeanVarSwarmBar(6,piX_TrPo{sp},1,0.05,[191/255, 191/255, 191/255]);
        title(sp);
    end
    linkaxes
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('%s Aligned', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% New Fig 2: Diagonal vs Off-Diagonal w/in epochs
for al = 1:length(alignments)
    tr_Diag = trial_Window_TrTr_TsTr_Diag(:,:,al);
    tr_OffDiag = trial_Window_TrTr_TsTr_OffDiag(:,:,al);
    itd_Diag = trial_Window_TrITD_TsITD_Diag(:,:,al);
    itd_OffDiag = trial_Window_TrITD_TsITD_OffDiag(:,:,al);
    tr_Diag = tr_Diag(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    tr_OffDiag = tr_OffDiag(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    itd_Diag = itd_Diag(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    itd_OffDiag = itd_OffDiag(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    
    figure;
    subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,reshape(spIDs(:,1:end-1), [1,mlb{ani}.seqLength^2])); 
    mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(tr_Diag),1,0.05, 'r');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(tr_OffDiag),1,0.05, 'r');
    mlb{ani}.PlotMeanVarSwarmBar(4,cell2mat(itd_Diag),1,0.05, 'b');
    mlb{ani}.PlotMeanVarSwarmBar(5,cell2mat(itd_OffDiag),1,0.05, 'b');
    set(gca, 'xtick', [1,2,4,5], 'xticklabel', [{'Diag'}, {'OffDiag'}, {'Diag'}, {'OffDiag'}]);
    title('All Positions')
    for sp = 1:mlb{ani}.seqLength
        subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,spIDs(sp,end));
        mlb{ani}.PlotMeanVarSwarmBar(1,tr_Diag{sp},1,0.05, 'r');
        hold on;
        mlb{ani}.PlotMeanVarSwarmBar(2,tr_OffDiag{sp},1,0.05, 'r');
        mlb{ani}.PlotMeanVarSwarmBar(4,itd_Diag{sp},1,0.05, 'b');
        mlb{ani}.PlotMeanVarSwarmBar(5,itd_OffDiag{sp},1,0.05, 'b');
        title(sp);
    end
    linkaxes    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', alignments{al},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% New Fig 3: Pre/Post w/in vs cross decode
for al = 1:length(alignments)
    win_PrPr = trial_Window_TrITD_TsITD(:,:,al);
    win_PoPo = trial_Window_TrAltITD_TsAltITD(:,:,al);
    x_PrPo = trial_Window_TrTr_TsPo(:,:,al);
    x_PoPr = trial_Window_TrPo_TsPr(:,:,al);    
    win_PrPr = win_PrPr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    win_PoPo = win_PoPo(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    x_PrPo = x_PrPo(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    x_PoPr = x_PoPr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    
    spIDs = reshape(1:mlb{ani}.seqLength*(mlb{ani}.seqLength+1), [mlb{ani}.seqLength+1, mlb{ani}.seqLength])';
    figure;
    subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,reshape(spIDs(:,1:end-1), [1,mlb{ani}.seqLength^2]));
    mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(win_PrPr),1,0.05,'r');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(win_PoPo),1,0.05,'b');
    mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(x_PrPo),1,0.05,[191/255, 191/255, 0]);
    mlb{ani}.PlotMeanVarSwarmBar(4,cell2mat(x_PoPr),1,0.05,[191/255, 191/255, 0]);
    if strcmp(alignments{al}, 'PokeIn')
        set(gca, 'xtick', 1:4, 'xticklabel', [{'Pre-Trial'}, {'Post-Trial'}, {'Pre->Post'}, {'Post->Pre'}], 'xticklabelrotation', 45);
    else
        set(gca, 'xtick', 1:4, 'xticklabel', [{'Post-Trial'}, {'Pre-Trial'}, {'Pre->Post'}, {'Post->Pre'}], 'xticklabelrotation', 45);
    end
    for sp = 1:mlb{ani}.seqLength
        subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,spIDs(sp,end));
        mlb{ani}.PlotMeanVarSwarmBar(1,win_PrPr{sp},1,0.05,'r');
        hold on;
        mlb{ani}.PlotMeanVarSwarmBar(2,win_PoPo{sp},1,0.05,'b');
        mlb{ani}.PlotMeanVarSwarmBar(3,x_PrPo{sp},1,0.05,[191/255, 191/255, 0]);
        mlb{ani}.PlotMeanVarSwarmBar(4,x_PoPr{sp},1,0.05,[191/255, 191/255, 0]);
        title(sp);
    end
    linkaxes
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', alignments{al},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% New Fig 4: ITD XTD By Lag
% NOTE: Set to only take inter-mediate positions that have both previous & upcoming positions
for al = 1:length(alignments)
    tempPrPo = trial_Window_TrPr_TsPo(:,:,al);
    next_PrPo = tempPrPo(triu(true(size(tempPrPo)),1) & tril(true(size(tempPrPo)),1));
    curr_PrPo = tempPrPo(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    prev_PrPo = tempPrPo(tril(true(size(tempPrPo)),-1) & triu(true(size(tempPrPo)),-1));
    tempPoPr = trial_Window_TrPo_TsPr(:,:,al);
    next_PoPr = tempPoPr(triu(true(size(tempPoPr)),1) & tril(true(size(tempPoPr)),1));
    curr_PoPr = tempPoPr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
    prev_PoPr = tempPoPr(tril(true(size(tempPoPr)),-1) & triu(true(size(tempPoPr)),-1));
    
    figure;
    subplot(2,1,1)
    mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(prev_PrPo(1:end-1)),1,0.05,'k');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(curr_PrPo(2:3)),1,0.05,'k');
    mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(next_PrPo(2:3)),1,0.05,'k');
    set(gca, 'xtick', 1:3, 'xticklabel', [{'Previous'}, {'Current'}, {'Next'}]);
    title('Pre->Post');
    
    subplot(2,1,2)
    mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(prev_PoPr(1:end-1)),1,0.05,'k');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(curr_PoPr(2:3)),1,0.05,'k');
    mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(next_PoPr(2:3)),1,0.05,'k');
    set(gca, 'xtick', 1:3, 'xticklabel', [{'Previous'}, {'Current'}, {'Next'}]);
    title('Post->Pre');
end

%% New Fig 5: Persistance Model Fit
% Plot poke in aligned fits
for al = 1:length(alignments)
    tempTrialFits = trial_TrialPersFit(:,:,al);
%     tempTrialFits = tempTrialFits(2:end,2:end);
    tempTrialFits = tempTrialFits(1:end-1,1:end-1);
%     tempTrialFits = tempTrialFits(2:end-1,2:end-1);
    tempITIfits = trial_ITIpersFit(:,:,al);
%     tempITIfits = tempITIfits(2:end,2:end);
    tempITIfits = tempITIfits(1:end-1,1:end-1);
%     tempITIfits = tempITIfits(2:end-1,2:end-1);
    tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
    tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
    itiDur = find(sum(~isnan(tempITIfits),1)./size(tempITIfits,1)==1,1,'last');
    
    [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
    minTrialLat = minTrialLat./(1/dsRate);
    [minITIfit,minITIlat] = min(tempITIfits,[],2);
    minITIlat = minITIlat./(1/dsRate);
    
    figure;
    spIDs = reshape(1:4^2, [4, 4])';
    sp(1) = subplot(4,4,reshape(spIDs(1:end-1,:), [numel(spIDs(:,2:end)),1]));
    trlPlot = mlb{ani}.PlotMeanVarLine((1:itiDur)./(1/dsRate), tempTrialFits(:,1:itiDur),1,0.05,'r');
    hold on;
    itiPlot = mlb{ani}.PlotMeanVarLine((1:itiDur)./(1/dsRate), tempITIfits(:,1:itiDur),1,0.05,'k');
%     scatter(minTrialLat, minTrialFit, 'or', 'filled', 'markerfacealpha', 0.2); 
%     scatter(minITIlat, minITIfit, 'ok', 'filled', 'markerfacealpha', 0.2);
    set(gca, 'xlim', [0 itiDur*dsRate]);
    legend([trlPlot, itiPlot], [{'Trial'}, {'ITD'}])
%     
%     sp(2) = subplot(4,4,spIDs(1:end-1,1));
%     histogram(minTrialFit, 0.3:0.05:1.2, 'FaceColor', 'r', 'orientation', 'horizontal', 'FaceAlpha', 0.5);
%     hold on;
%     histogram(minITIfit, 0.3:0.05:1.2, 'FaceColor', 'k', 'orientation', 'horizontal', 'FaceAlpha', 0.5);
    
    sp(3) = subplot(4,4,spIDs(end,:));
    histogram(minTrialLat, 0:dsRate:itiDur*dsRate,'FaceColor', 'r', 'FaceAlpha', 0.5);
    hold on;
    histogram(minITIlat, 0:dsRate:itiDur*dsRate,'FaceColor', 'k', 'FaceAlpha', 0.5);
    set(gca, 'xlim', [0 itiDur*dsRate]);
    
%     linkaxes(sp(1:2), 'y');
    linkaxes(sp([1,3]),'x');
end
%% Plot Lag Values
for al = 1:length(alignments)
%     colLim = [-0.5 0.5];
    colLim = [-0.75 0.75];
%     colLim = [-1 1];
    figure;
    tempTMat = trialD(:,:,al);
%     tempTMat = trialDz(:,:,al);
%     tempTMat = trialDposZ(:,:,al);
    
    subplot(2,3,1)
    prevPosD = tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1));
    mlb{end}.PlotTrialPDM(cell2mat(permute(prevPosD(1:end-1), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%     mlb{end}.PlotTrialPDM(cell2mat(permute(prevPosD, [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    set(gca, 'clim', colLim );
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    title('Previous Position');
    subplot(2,3,2)
    currPosD = tempTMat(logical(eye(length(tempTMat))));
    mlb{end}.PlotTrialPDM(cell2mat(permute(currPosD(2:end-1), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%     mlb{end}.PlotTrialPDM(cell2mat(permute(currPosD(2:end), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%     mlb{end}.PlotTrialPDM(cell2mat(permute(currPosD, [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    set(gca, 'clim', colLim );
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    title('Current Position');
    subplot(2,3,3)
    nextPosD = tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1));
    mlb{end}.PlotTrialPDM(cell2mat(permute(nextPosD(2:end), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%     mlb{end}.PlotTrialPDM(cell2mat(permute(nextPosD, [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    set(gca, 'clim', colLim);
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
%     set(gca,'clim', [-0.5 0.5]);
    title('Next Position');
    colormap(cMap)
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    subplot(2,2,3);    
    prevPosPoPr = trial_Window_TrPo_TsPr(:,:,al);
    prevPosPoPr = prevPosPoPr(triu(ones(length(trial_Window_TrPo_TsPr)),-1) & tril(ones(length(trial_Window_TrPo_TsPr)),-1));
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(prevPosPoPr(1:end-1)),1,0.05,'k');
    curPosPoPr = trial_Window_TrPo_TsPr(:,:,al);
    curPosPoPr = curPosPoPr(logical(eye(length(trial_Window_TrPo_TsPr))));
    mlb{ani}.PlotMeanVarSwarmBar(2, cell2mat(curPosPoPr(2:end-1)),1,0.05,'k');
    nextPosPoPr = trial_Window_TrPo_TsPr(:,:,al);
    nextPosPoPr = nextPosPoPr(triu(ones(length(trial_Window_TrPo_TsPr)),1) & tril(ones(length(trial_Window_TrPo_TsPr)),1));
    mlb{ani}.PlotMeanVarSwarmBar(3, cell2mat(nextPosPoPr(2:end)),1,0.05,'k');
    set(gca,'xtick', 1:3, 'xticklabels', [{'Prev'},{'Current'},{'Next'}]);
    title([{'Train:Post'};{'Test:Pre'}]);
    
    subplot(2,2,4);
    prevPosPrPo = trial_Window_TrPr_TsPo(:,:,al);
    prevPosPrPo = prevPosPrPo(triu(ones(length(trial_Window_TrPr_TsPo)),-1) & tril(ones(length(trial_Window_TrPr_TsPo)),-1));
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(prevPosPrPo(1:end-1)),1,0.05,'k');
    curPosPrPo = trial_Window_TrPr_TsPo(:,:,al);
    curPosPrPo = curPosPrPo(logical(eye(length(trial_Window_TrPr_TsPo))));
    mlb{ani}.PlotMeanVarSwarmBar(2, cell2mat(curPosPrPo(2:end-1)),1,0.05,'k');
    nextPosPrPo = trial_Window_TrPr_TsPo(:,:,al);
    nextPosPrPo = nextPosPrPo(triu(ones(length(trial_Window_TrPr_TsPo)),1) & tril(ones(length(trial_Window_TrPr_TsPo)),1));
    mlb{ani}.PlotMeanVarSwarmBar(3, cell2mat(nextPosPrPo(2:end)),1,0.05,'k');
    set(gca,'xtick', 1:3, 'xticklabels', [{'Prev'},{'Current'},{'Next'}]);
    title([{'Train:Pre'};{'Test:Post'}]);
end
%% Future-Past Decodability
% for al = 1:length(alignments)
%     figure;
%     spIDs = reshape(1:(mlb{ani}.seqLength-2)*3,[3,mlb{ani}.seqLength-2])';
%     tempTMat = trialD(:,:,al);
%     diff = cell(1,1,mlb{ani}.seqLength-2);
%     for p = 2:mlb{ani}.seqLength-1
%         diff{1,1,p-1} = tempTMat{p,p-1} - tempTMat{p,p+1};
%         subplot(2,3,spIDs(p-1));
%         mlb{end}.PlotTrialPDM(diff{1,1,p-1}, 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%         title(sprintf('Pos %i',p));
%         hold on;
%         if strcmp(alignments{al} , 'PokeIn')
%             plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
%             plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
%             plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
%             
%             plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
%             plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%             plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%         elseif strcmp(alignments{al}, 'PokeOut')
%             plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
%             plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
%             plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
%             
%             plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
%             plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%             plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%         end
%         colormap(cMap);
%     end
%     subplot(2,3,reshape(spIDs(:,2:end), [1,numel(spIDs(:,2:end))]));
%     mlb{end}.PlotTrialPDM(cell2mat(diff), 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
%     hold on;
%     if strcmp(alignments{al} , 'PokeIn')
%         plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
%         
%         plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
%         plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%         plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     elseif strcmp(alignments{al}, 'PokeOut')
%         plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
%         plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
%         
%         plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
%         plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
%         plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
%     end
%     colormap(cMap);
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', sprintf('Decodability of Prev - Next : %s Alignment', alignments{al}),...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% end
%% Plot XTD by Compartment
for al = 1:length(alignments)
%     posVect = 1:mlb{ani}.seqLength;
%     posVect = 1:mlb{ani}.seqLength-1;
%     posVect = 2:3;
    posVect = 2:mlb{ani}.seqLength;
%     posVect = 1;
    figure;
    % PoPr
    subplot(3,3,1);
    curPosPoPr = trial_Window_TrPo_TsPr(:,:,al);
    curPosPoPr = curPosPoPr(logical(eye(length(curPosPoPr))));
    curPosPoPr = curPosPoPr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPoPr),1,0.05,'k');
    title([{'Train:Post'};{'Test:Pre'}]);
    % PoTr
    subplot(3,3,2);
    curPosPoTr = trial_Window_TrPo_TsTr(:,:,al);
    curPosPoTr = curPosPoTr(logical(eye(length(curPosPoTr))));
    curPosPoTr = curPosPoTr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPoTr),1,0.05,'k');
    title([{'Train:Post'};{'Test:Trial'}]);
    % PoPo
    subplot(3,3,3);
    curPosPoPo = trial_Window_TrPo_TsPo(:,:,al);
    curPosPoPo = curPosPoPo(logical(eye(length(curPosPoPo))));
    curPosPoPo = curPosPoPo(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPoPo),1,0.05,'k');
    title([{'Train:Post'};{'Test:Post'}]);
    % TrPr
    subplot(3,3,4);
    curPosTrPr = trial_Window_TrTr_TsPr(:,:,al);
    curPosTrPr = curPosTrPr(logical(eye(length(curPosTrPr))));
    curPosTrPr = curPosTrPr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosTrPr),1,0.05,'k');
    title([{'Train:Trial'};{'Test:Pre'}]);
    % TrTr
    subplot(3,3,5);
    curPosTrTr = trial_Window_TrTr_TsTr(:,:,al);
    curPosTrTr = curPosTrTr(logical(eye(length(curPosTrTr))));
    curPosTrTr = curPosTrTr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosTrTr),1,0.05,'k');
    title([{'Train:Trial'};{'Test:Trial'}]);
    % TrPo
    subplot(3,3,6);
    curPosTrPo = trial_Window_TrTr_TsPo(:,:,al);
    curPosTrPo = curPosTrPo(logical(eye(length(curPosTrPo))));
    curPosTrPo = curPosTrPo(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosTrPo),1,0.05,'k');
    title([{'Train:Trial'};{'Test:Post'}]);
    % PrPr
    subplot(3,3,7);
    curPosPrPr = trial_Window_TrPr_TsPr(:,:,al);
    curPosPrPr = curPosPrPr(logical(eye(length(curPosPrPr))));
    curPosPrPr = curPosPrPr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPrPr),1,0.05,'k');
    title([{'Train:Pre'};{'Test:Pre'}]);
    % PrTr
    subplot(3,3,8);
    curPosPrTr = trial_Window_TrPr_TsTr(:,:,al);
    curPosPrTr = curPosPrTr(logical(eye(length(curPosPrTr))));
    curPosPrTr = curPosPrTr(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPrTr),1,0.05,'k');
    title([{'Train:Pre'};{'Test:Trial'}]);
    % PrPo
    subplot(3,3,9);
    curPosPrPo = trial_Window_TrPr_TsPo(:,:,al);
    curPosPrPo = curPosPrPo(logical(eye(length(curPosPrPo))));
    curPosPrPo = curPosPrPo(posVect);
    mlb{ani}.PlotMeanVarSwarmBar(1, cell2mat(curPosPrPo),1,0.05,'k');
    title([{'Train:Pre'};{'Test:Post'}]);
    linkaxes;
    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
end
%% Plot Session Trial Windows
% % trlWindows{1} = [-700 2000];
% for ani = 1:length(mlb)
%     curAniBehav = mlb{ani}.trialInfo;
%     trlPokeIns = [curAniBehav.PokeInIndex];
%     trlPokeOuts = [curAniBehav.PokeOutIndex];
%     seqTransIndex = logical(diff([curAniBehav.SequenceNum]));
%     isLog = [curAniBehav.TranspositionDistance]==0;
%     itiDur = trlPokeIns(2:end)-trlPokeOuts(1:end-1);
%     figure; 
%     subplot(2,2,1:2);
%     for trl = 1:length(trlPokeIns)
%         patch('XData', [trlPokeIns(trl)+trlWindows{1}(1), trlPokeIns(trl)+trlWindows{1}(2),trlPokeIns(trl)+trlWindows{1}(2),trlPokeIns(trl)+trlWindows{1}(1)],...
%             'YData', [-0.5 -0.5 0.5 0.5],...
%             'FaceColor', 'k', 'FaceAlpha', 0.2);
%     end
%     hold on;
%     scatter(trlPokeIns(isLog),zeros(1,length(trlPokeIns(isLog))), '|k');
%     scatter(trlPokeIns(~isLog),zeros(1,length(trlPokeIns(~isLog))), 'xr');
%     set(gca, 'ylim', [-2 2]);
%     title(mlb{ani}.pathDir, 'interpreter', 'none');
%     subplot(2,2,3);
%     histogram(itiDur(seqTransIndex)/1000, 0:1:20);
%     subplot(2,2,4);
%     histogram(itiDur(~seqTransIndex)/1000, 0:0.2:5)
%     title(sprintf('Mean = %.02f', mean(itiDur(~seqTransIndex))));
%     
% end 
%% Plot False Alarm Rate
for al = 1:length(alignments)
    figure;
    curFAR = faRate(:,:,:,al);
    for ani = 1:length(fileDirs)
        aniFAR = curFAR(:,:,ani);
        posFAR = aniFAR(logical(eye(size(aniFAR))));
        tempAniFAR = [];
        for pos = 1:mlb{ani}.seqLength
            subplot(mlb{ani}.seqLength+1,length(fileDirs), sub2ind([mlb{ani}.seqLength+1,length(fileDirs)],ani,pos))
            mlb{end}.PlotTrialPDM(posFAR{pos}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            tempAniFAR(:,:,pos) = mean(posFAR{pos},3,'omitnan'); %#ok<SAGROW>
%             tempAniFAR(:,:,pos) = max(posFAR{pos},[],3); %#ok<SAGROW>
        end
        meanDecode = nan(size(tempAniFAR(:,:,1)));
        for p = 1:size(tempAniFAR,1)
            for pp = 1:size(tempAniFAR,2)
                meanDecode(p,pp) = find(tempAniFAR(p,pp,:)==max(tempAniFAR(p,pp,:)),1,'first');
            end
        end
        subplot(mlb{ani}.seqLength+1,length(fileDirs), sub2ind([mlb{ani}.seqLength+1,length(fileDirs)],ani,mlb{ani}.seqLength+1))
        imagesc(obsvTimeVect{al}, obsvTimeVect{al}, meanDecode');
        xlabel('Test Time');
        ylabel('Train Time');
        colorbar
        set(gca,'ydir', 'normal')
    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('False Alarm Rate: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    end
end

%% Plot Hit Rate
for al = 1:length(alignments)
    figure;
    curHR = hitRate(:,:,:,al);
    for ani = 1:length(fileDirs)
        aniHR = curHR(:,:,ani);
        posHR = aniHR(logical(eye(size(aniHR))));
        tempAniFAR = [];
        for pos = 1:mlb{ani}.seqLength
            subplot(mlb{ani}.seqLength+1,length(fileDirs), sub2ind([mlb{ani}.seqLength+1,length(fileDirs)],ani,pos))
            mlb{end}.PlotTrialPDM(posHR{pos}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
            tempAniHR(:,:,pos) = mean(posHR{pos},3,'omitnan'); %#ok<SAGROW>
%             tempAniHR(:,:,pos) = max(posHR{pos},[],3); %#ok<SAGROW>
        end
        meanDecode = nan(size(tempAniHR(:,:,1)));
        for p = 1:size(tempAniHR,1)
            for pp = 1:size(tempAniHR,2)
                meanDecode(p,pp) = find(tempAniHR(p,pp,:)==max(tempAniHR(p,pp,:)),1,'first');
            end
        end
        subplot(mlb{ani}.seqLength+1,length(fileDirs), sub2ind([mlb{ani}.seqLength+1,length(fileDirs)],ani,mlb{ani}.seqLength+1))
        imagesc(obsvTimeVect{al}, obsvTimeVect{al}, meanDecode');
        xlabel('Test Time');
        ylabel('Train Time');
        colorbar
        set(gca,'ydir', 'normal')
    
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Hit Rate: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    end
end

%% Evaluate similarity in XTD by position
% for al = 1:length(alignments)
%     tempTmat = trialD(:,:,al);
%     iscXTD = tempTmat (logical(eye(size(tempTmat))));
%     iscXTDvect = cell(1,size(iscXTD,1));
%     trlIDvect = cellfun(@(a){ones(1,size(a,3))},iscXTD)';
%     for pos = 1:length(iscXTD)
%         trlIDvect{pos} = trlIDvect{pos}.*pos;
%         tempXTDvect = nan(numel(iscXTD{pos}(:,:,1)),size(iscXTD{pos},3));
%         for trl = 1:size(iscXTD{pos},3)
%             tempXTDvect(:,trl) = reshape(iscXTD{pos}(:,:,trl),[numel(iscXTD{pos}(:,:,trl)),1]);
%         end
%         iscXTDvect{pos} = tempXTDvect;
%     end
%     iscXTDcorrMTX = corr(cell2mat(iscXTDvect));
%     trlIDvect = cell2mat(trlIDvect);
%     
%     figure;
%     imagesc(iscXTDcorrMTX);
%     hold on; plot(repmat(find(trlIDvect==2,1,'first'),[1,2]), get(gca,'ylim'),'-k')
%     hold on; plot(repmat(find(trlIDvect==3,1,'first'),[1,2]), get(gca,'ylim'),'-k')
%     hold on; plot(repmat(find(trlIDvect==4,1,'first'),[1,2]), get(gca,'ylim'),'-k')
%     hold on; plot(get(gca,'xlim'), repmat(find(trlIDvect==4,1,'first'),[1,2]),'-k')
%     hold on; plot(get(gca,'xlim'), repmat(find(trlIDvect==3,1,'first'),[1,2]),'-k')
%     hold on; plot(get(gca,'xlim'), repmat(find(trlIDvect==2,1,'first'),[1,2]),'-k')
%     hold on; plot(get(gca,'xlim'), repmat(find(trlIDvect==2,1,'first'),[1,2]),'-k')
%     
%     figure;
%     for p1 = 1:4
%         for p2 = 1:4
%             blankMask = false(size(iscXTDcorrMTX));
%             blankMask(trlIDvect==p1, trlIDvect==p2) = true;
%             blankMask(logical(eye(length(blankMask)))) = false;
%             subplot(4,4,sub2ind([4,4],p2,p1));
%             mlb{ani}.PlotMeanVarSwarmBar(1, reshape(iscXTDcorrMTX(blankMask),[sum(blankMask(:)),1]),1,0.05,'k');
%         end
%     end
%     linkaxes
% end
%% Evaluate Trial Train Pre vs Post
for al = 1:length(alignments)
%     posVect = 1:mlb{ani}.seqLength; tit = 'All';
%     posVect = 1:mlb{ani}.seqLength-1; tit = 'All but Final';
%     posVect = 2:3;    tit = 'Pos 2&3';
    posVect = 2:mlb{ani}.seqLength;   tit = 'SANSA';
%     posVect = 4; tit = sprintf('Pos=%i',posVect);
    figure;
    % TrPr
    curPosTrPr = trial_Window_TrTr_TsPr(:,:,al);
    curPosTrPr = curPosTrPr(logical(eye(length(curPosTrPr))));
    curPosTrPr = cell2mat(curPosTrPr(posVect));
    % TrPo
    curPosTrPo = trial_Window_TrTr_TsPo(:,:,al);
    curPosTrPo = curPosTrPo(logical(eye(length(curPosTrPo))));
    curPosTrPo = cell2mat(curPosTrPo(posVect));
    
    trTRdiffPoPr = curPosTrPo-curPosTrPr;    
    mlb{ani}.PlotMeanVarSwarmBar(1,trTRdiffPoPr,1,0.05,'k', 'error', 'CI');
    [~,p(1),~,stats(1)] = ttest(trTRdiffPoPr);
    
    % PoTr
    curPosPoTr = trial_Window_TrPo_TsTr(:,:,al);
    curPosPoTr = curPosPoTr(logical(eye(length(curPosPoTr))));
    curPosPoTr = cell2mat(curPosPoTr(posVect));
    
    % PrTr
    curPosPrTr = trial_Window_TrPr_TsTr(:,:,al);
    curPosPrTr = curPosPrTr(logical(eye(length(curPosPrTr))));
    curPosPrTr = cell2mat(curPosPrTr(posVect));
    
    trITDdiffTr = curPosPoTr-curPosPrTr;
    mlb{ani}.PlotMeanVarSwarmBar(2,trITDdiffTr,1,0.05,'k', 'error', 'CI');    
    [~,p(2),~,stats(2)] = ttest(trITDdiffTr);
    
    % PoPr
    curPosPoPr = trial_Window_TrPo_TsPr(:,:,al);
    curPosPoPr = curPosPoPr(logical(eye(length(curPosPoPr))));
    curPosPoPr = cell2mat(curPosPoPr(posVect));
    % PrPo
    curPosPrPo = trial_Window_TrPr_TsPo(:,:,al);
    curPosPrPo = curPosPrPo(logical(eye(length(curPosPrPo))));
    curPosPrPo = cell2mat(curPosPrPo(posVect));

    trlITDdiffITD = curPosPoPr-curPosPrPo;
    mlb{ani}.PlotMeanVarSwarmBar(3, trlITDdiffITD,1,0.05,'k', 'error', 'CI'); 
    [~,p(3),~,stats(3)] = ttest(trlITDdiffITD);
    
    set(gca,'xtick', 1:3, 'xticklabels',...
        [{'Train Trial: Decode ITD (Post-Pre)'},...
        {'Decode Trial: Train ITD (Post-Pre)'},...
        {'X Decode ITD (Train Post-Pre)'}],...
        'xticklabelrotation', 45);
    title(tit);
    
    set(gca,'ylim', [min(get(gca,'ylim')), max(get(gca,'ylim'))+2]);
    for tst = 1:3
        if p(tst)<0.05
            text(tst,max(get(gca,'ylim'))-1, [{sprintf('t_{(%i)}=%.02f',stats(tst).df, stats(tst).tstat)};{sprintf('p=%.02i',p(tst))}],...
                'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color','r');
        else
            text(tst,max(get(gca,'ylim'))-1, [{sprintf('t_{(%i)}=%.02f',stats(tst).df, stats(tst).tstat)};{sprintf('p=%.02i',p(tst))}],...
                'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color','k');
        end
    end
end
%% Future-Past Decodability
for al = 1:length(alignments)
    figure;
    spIDs = reshape(1:(mlb{ani}.seqLength-2)*3,[3,mlb{ani}.seqLength-2])';
    tempTMat = trialD(:,:,al);
%     tempTMat = trialDz(:,:,al);
    diff = cell(1,1,mlb{ani}.seqLength-2);
    for p = 2:mlb{ani}.seqLength-1
        diff{1,1,p-1} = tempTMat{p,p-1} - tempTMat{p,p+1};
        subplot(2,3,spIDs(p-1));
        mlb{end}.PlotTrialPDM(diff{1,1,p-1}, 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
        title(sprintf('Pos %i',p));
        hold on;
        if strcmp(alignments{al} , 'PokeIn')
            plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
            
            plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
        elseif strcmp(alignments{al}, 'PokeOut')
            plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
            plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
            
            plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
        end
        colormap(cMap);
    end
    subplot(2,3,reshape(spIDs(:,2:end), [1,numel(spIDs(:,2:end))]));
    mlb{end}.PlotTrialPDM(cell2mat(diff), 'rotate', 'clim', [-2 2], 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
    hold on;
    if strcmp(alignments{al} , 'PokeIn')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    elseif strcmp(alignments{al}, 'PokeOut')
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
        
        plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
    end
    colormap(cMap);
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability of Prev - Next : %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% Evaluate Transitions During Trial & ITI Periods
%     posVect = 1:mlb{ani}.seqLength; tit = 'All';
%     posVect = 1:mlb{ani}.seqLength-1; tit = 'All but Final';
%     posVect = 2:3;    tit = 'Pos 2&3';
    posVect = 2:mlb{ani}.seqLength;   tit = 'SANSA';
%     posVect = 3; tit = sprintf('Pos=%i',posVect);
    
    tempTmat = trialD(:,:,al);
%     tempTmat = trialDz(:,:,al);
%     tempTmat = trialDposZ(:,:,al);
%     tempTmat = trialHR(:,:,al);
%     tempTmat = trialFAR(:,:,al);

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
        preTrlBound = find(obsvTimeVect{al}<=grpPoPiLat,1,'last')-(binSize/dsRate/2);
        preTrialLog = blankTrlLogVect;
        preTrialLog(1:preTrlBound) = true;
        
        trialLowBound = find(obsvTimeVect{al}<=grpPoPiLat,1,'last')+(binSize/dsRate/2);
        trialUpBound = find(obsvTimeVect{al}==0)-(binSize/dsRate/2);
        trialLog = blankTrlLogVect;
        trialLog(trialLowBound:trialUpBound) = true;
        
        pstTrlBound = find(obsvTimeVect{al}==0)+(binSize/dsRate/2);
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
    
    subplot(6,4,[1,5]);
    prevTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(prevTrPr,2),3,0.05,'b');
    curTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(trTr_TsPr,2),3,0.05,'r'); 
    nxtTP = mlb{ani}.PlotMeanVarLine(obsvTimeVect{al}(1):dsRate/1000:obsvTimeVect{al}(preTrlBound),mean(nxtTrPr,2),3,0.05,'g');
    legend([prevTP, curTP,nxtTP], [{'Prev'}, {'Current'}, {'Next'}]);
%     legend([prevTP,nxtTP], [{'Prev'}, {'Next'}]);
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
%% Trial/ITI Fits by Position
trialPers = trial_TrialPersFit(logical(eye(length(trial_TrialPersFit))));
itiPers = trial_ITIpersFit(logical(eye(length(trial_ITIpersFit))));

trlLat = cell(size(trialPers));
itiLat = cell(size(itiPers));
posLog = cell(size(trialPers));
figure;
for pos = 1:size(trialPers,1)
    [~,itiLat{pos}] = min(itiPers{pos},[],2);
    tempITIplt = mlb{ani}.PlotMeanVarSwarmBar(((pos-1)*2)+1, itiLat{pos},1,0.05,mlb{ani}.PositionColors(pos,:), 'error', 'CI');
    set(tempITIplt, 'linewidth', 1.5, 'linestyle','--');
    hold on;
    [~,trlLat{pos}] = min(trialPers{pos},[],2);
    tempTrlPlt = mlb{ani}.PlotMeanVarSwarmBar(((pos-1)*2)+2, trlLat{pos},1,0.05,mlb{ani}.PositionColors(pos,:), 'error', 'CI');
    set(tempTrlPlt , 'linewidth', 1.5, 'linestyle','-');
    posLog{pos} = ones(size(trialPers{pos},1),1).*pos;
end
set(gca,'xtick', 1:2, 'xticklabels', [{'ITD'}, {'Trial'}]);
anovaDta = [cell2mat(itiLat);cell2mat(trlLat)];
groupLog = [repmat(cell2mat(posLog),[2,1]), [ones(size(cell2mat(posLog))); ones(size(cell2mat(posLog))).*2]];
[p,tbl,stats] = anovan(anovaDta,groupLog, 'model', 'interaction', 'varnames', {'Position', 'Epoch'}, 'display','off');
figure;
multcompare(stats,'Dimension', [1,2]);

% trialLog = trial_TAOlog(logical(eye(length(trial_TAOlog))));
trialLog = trial_SeqOSlog(logical(eye(length(trial_SeqOSlog))));
trialLog = cell2mat(trialLog(4));
tLat = cell2mat(trlLat(4));
iLat = cell2mat(itiLat(4));
figure;
mlb{ani}.PlotMeanVarSwarmBar(1,tLat(trialLog),1,0.05,'k','error', 'CI');
hold on;
mlb{ani}.PlotMeanVarSwarmBar(2,tLat(~trialLog),1,0.05,'r','error','CI');
mlb{ani}.PlotMeanVarSwarmBar(4,iLat(trialLog),1,0.05,'k','error','CI');
mlb{ani}.PlotMeanVarSwarmBar(5,iLat(~trialLog),1,0.05,'r','error','CI');



%% Save output
% save('PFC_XTD_Windows_PosChance.mat', '-v7.3');