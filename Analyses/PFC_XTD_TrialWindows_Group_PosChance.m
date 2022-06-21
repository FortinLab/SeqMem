% PFC_XTD_TrialWindows_Group_PosChance

%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

% % CA1 Data
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% tets = [1,22,17,18,17]; % Lateral/Distal
% % tets = [7,3,1,5,5]; % Medial/Proximal

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
binSize = 200;
dsRate = 50;
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

% alignments = [{'PokeIn'}, {'PokeOut'}];
% trlWindows = [{[-1500 2000]}, {[-2000 2000]}];
alignments = {'PokeIn'};
trlWindows = {[-1200 2500]};
    
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
        % Posteriors
        trlD = cell(mlb{ani}.seqLength, mlb{ani}.seqLength, length(fileDirs), length(alignments));
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
    piPokeOutLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeOutIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeInIndex])'/1000;
    piRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeInIndex])'/1000;
    poPokeInLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeInIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeOutIndex])'/1000;
    poRwdLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials).RewardIndex] - [mlb{ani}.trialInfo(mlb{ani}.fiscTrials).PokeOutIndex])'/1000;
    nxtTrlLat{ani} = ([mlb{ani}.trialInfo(mlb{ani}.fiscTrials(1:end-1,:)+1).PokeInIndex]-[mlb{ani}.trialInfo(mlb{ani}.fiscTrials(1:end-1,:)).PokeOutIndex])'/1000;
    
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
                tempWindow_TrTr_TsTr = nan(size(hits,3),1);                     % Train & Test Trial Window
                tempWindow_TrTr_TsTr_Diag = nan(size(hits,3),1);                % Diagonal of the Train & Test Trial Window
                tempWindow_TrTr_TsTr_OffDiag = nan(size(hits,3),1);             % Off Diagonal of the Train & Test Trial Window
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
            trial_Window_TrPr_TsTr{r,c,al} = cell2mat(permute(window_TrPr_TsTr(r,c,:,al), [3,1,2]));            
            trial_Window_TrTr_TsPr{r,c,al} = cell2mat(permute(window_TrTr_TsPr(r,c,:,al), [3,1,2]));            
            trial_Window_TrPr{r,c,al} = cell2mat(permute(window_TrPr(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsPo{r,c,al} = cell2mat(permute(window_TrTr_TsPo(r,c,:,al), [3,1,2]));
            trial_Window_TrPo_TsTr{r,c,al} = cell2mat(permute(window_TrPo_TsTr(r,c,:,al), [3,1,2]));
            trial_Window_TrPo{r,c,al} = cell2mat(permute(window_TrPo(r,c,:,al), [3,1,2]));
            trial_Window_TrPo_TsPr{r,c,al} = cell2mat(permute(window_TrPo_TsPr(r,c,:,al), [3,1,2]));
            trial_Window_TrPr_TsPo{r,c,al} = cell2mat(permute(window_TrPr_TsPo(r,c,:,al), [3,1,2]));
            trial_Window_PoPr{r,c,al} = cell2mat(permute(window_PoPr(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr{r,c,al} = cell2mat(permute(window_TrTr_TsTr(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr_Diag{r,c,al} = cell2mat(permute(window_TrTr_TsTr_Diag(r,c,:,al), [3,1,2]));
            trial_Window_TrTr_TsTr_OffDiag{r,c,al} = cell2mat(permute(window_TrTr_TsTr_OffDiag(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD{r,c,al} = cell2mat(permute(window_TrITD_TsITD(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD_Diag{r,c,al} = cell2mat(permute(window_TrITD_TsITD_Diag(r,c,:,al), [3,1,2]));
            trial_Window_TrITD_TsITD_OffDiag{r,c,al} = cell2mat(permute(window_TrITD_TsITD_OffDiag(r,c,:,al), [3,1,2]));
            trial_Window_TrAltITD_TsAltITD{r,c,al} = cell2mat(permute(window_TrlAltITD_TsAltITD(r,c,:,al), [3,1,2]));
            trial_TrialPersFit{r,c,al} = cell2mat(permute(trlPersFit(r,c,:,al), [3,1,2]));
            trial_ITIpersFit{r,c,al} = cell2mat(permute(itiPersFit(r,c,:,al), [3,1,2]));
        end
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
                trial_Window_TrAltITD_TsAltITD{p,p2,al}];
            trialWiseRawEpochs{p,p2,al} = tempTrialWindowed;
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
grpNxtTrlLat = median(cell2mat(nxtTrlLat));
%% Plot Trial D: All positions & decodings
for al = 1:length(alignments)
    figure;
    for r = 1:mlb{end}.seqLength
        for c = 1:mlb{end}.seqLength
            subplot(mlb{end}.seqLength, mlb{end}.seqLength, sub2ind([mlb{end}.seqLength, mlb{end}.seqLength], c,r))
            mlb{end}.PlotTrialPDM(trialD{r,c,al}, 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
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
%         tempTMat = tempTMat(2:end,2:end);        
%         tempTMat = tempTMat(2:end-1,2:end-1);
    elseif strcmp(alignments{al}, 'PokeOut')
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
% %% Plot Specific Window Comparisons
% selType = 0; % Alignment includes all trial positions
% % selType = 1; % Alignment specific trial positions included
% % selType = 2; % Alignment only common trial positions i.e. 2&3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Upcoming Trial Info
% figure;
% subplot(4,3,1)
% tempTMat = trialD(:,:,1);
% if selType == 1
%     tempTMat = tempTMat(2:end,2:end);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [0, 0, grpPiPoLat, grpPiPoLat],...
%     'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
%     'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
%     'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% title('Previous Position Info');
% 
% subplot(4,3,2)
% tempTMat = trialD(:,:,2);
% if selType == 1
%     tempTMat = tempTMat(1:end-1,1:end-1);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [grpPoPiLat, grpPoPiLat, 0, 0],...
%     'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% 
%     
% subplot(4,3,3)
% tempTrPr = trial_Window_TrPr(:,:,1);
% tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
% tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
% tempTrPo = trial_Window_TrPo(:,:,2);
% tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
% tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
% if selType == 1
%     tempTrPr = tempTrPr(2:end,2:end);
%     tempPrTr_Window = tempPrTr_Window(2:end,2:end);
%     tempTrPr_Window = tempTrPr_Window(2:end,2:end);
%     tempTrPo = tempTrPo(1:end-1,1:end-1);
%     tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
%     tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
% elseif selType == 2
%     tempTrPr = tempTrPr(2:end-1,2:end-1);
%     tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
%     tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
%     tempTrPo = tempTrPo(2:end-1,2:end-1);
%     tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
%     tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
% end
% trPr = cell2mat(tempTrPr(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% win_PrTr = cell2mat(tempPrTr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% win_TrPr = cell2mat(tempTrPr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% trPo = cell2mat(tempTrPo(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% win_PoTr = cell2mat(tempPoTr_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% win_TrPo = cell2mat(tempTrPo_Window(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1)));
% bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
% hold on;
% swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Current Trial
% subplot(4,3,4)
% tempTMat = trialD(:,:,1);
% if selType == 1
%     tempTMat = tempTMat(2:end,2:end);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [0, 0, grpPiPoLat, grpPiPoLat], 'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0],...
%     'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
%     'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
%     'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% title('Current Position Info');
% 
% subplot(4,3,5)
% tempTMat = trialD(:,:,2);
% if selType == 1
%     tempTMat = tempTMat(1:end-1,1:end-1);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [grpPoPiLat, grpPoPiLat, 0, 0], 'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255],...
%     'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% 
% subplot(4,3,6)
% tempTrPr = trial_Window_TrPr(:,:,1);
% tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
% tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
% tempTrPo = trial_Window_TrPo(:,:,2);
% tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
% tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
% if selType == 1
%     tempTrPr = tempTrPr(2:end,2:end);
%     tempPrTr_Window = tempPrTr_Window(2:end,2:end);
%     tempTrPr_Window = tempTrPr_Window(2:end,2:end);
%     tempTrPo = tempTrPo(1:end-1,1:end-1);
%     tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
%     tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
% elseif selType == 2
%     tempTrPr = tempTrPr(2:end-1,2:end-1);
%     tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
%     tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
%     tempTrPo = tempTrPo(2:end-1,2:end-1);
%     tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
%     tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
% end
% trPr = cell2mat(tempTrPr(logical(eye(length(tempTrPr)))));
% win_PrTr = cell2mat(tempPrTr_Window(logical(eye(length(tempTrPr)))));
% win_TrPr = cell2mat(tempTrPr_Window(logical(eye(length(tempTrPr)))));
% trPo = cell2mat(tempTrPo(logical(eye(length(tempTrPr)))));
% win_TrPo = cell2mat(tempTrPo_Window(logical(eye(length(tempTrPr)))));
% win_PoTr = cell2mat(tempPoTr_Window(logical(eye(length(tempTrPr)))));
% bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
% hold on;
% swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_PrTr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);
% % 
% % [h,p,ci,stats] = ttest(trPr, trPo);
% % [p,tbl,stats] = anova1([win_PrTr,win_TrPr,win_PoTr,win_TrPo]);
% % 
% % figure;
% % multcompare(stats, 'CType', 'bonferroni');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Next Trial
% subplot(4,3,7)
% tempTMat = trialD(:,:,1);
% if selType == 1
%     tempTMat = tempTMat(2:end,2:end);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [0, 0, grpPiPoLat, grpPiPoLat],...
%     'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, grpPiPoLat, grpPiPoLat, 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),0, 0],...
%     'facecolor', 'none', 'edgecolor', [191/255, 191/255, 0], 'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeIn')}), 0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'YData', [grpPiPoLat, grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [grpPiPoLat, max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), max(obsvTimeVect{strcmp(alignments, 'PokeIn')}), grpPiPoLat],...
%     'YData', [0, 0, min(obsvTimeVect{strcmp(alignments, 'PokeIn')}),  min(obsvTimeVect{strcmp(alignments, 'PokeIn')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% title('Next Position Info');
% 
% subplot(4,3,8)
% tempTMat = trialD(:,:,2);
% if selType == 1
%     tempTMat = tempTMat(1:end-1,1:end-1);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1)), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% patch('XData', [grpPoPiLat, 0, 0, grpPoPiLat],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [grpPoPiLat, grpPoPiLat, 0, 0],...
%     'facecolor', 'none', 'edgecolor', [196/255, 33/255, 255/255], 'linestyle', '--', 'linewidth', 5);
% patch('XData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat, min(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'YData', [0, 0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')})],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '-', 'linewidth', 5);
% patch('XData', [0, max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), max(obsvTimeVect{strcmp(alignments, 'PokeOut')}), 0],...
%     'YData', [min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), min(obsvTimeVect{strcmp(alignments, 'PokeOut')}), grpPoPiLat, grpPoPiLat],...
%     'facecolor', 'none', 'edgecolor', [191/255, 96/255, 0], 'linestyle', '--', 'linewidth', 5);
% 
% subplot(4,3,9)
% tempTrPr = trial_Window_TrPr(:,:,1);
% tempPrTr_Window = trial_Window_TrPr_TsTr(:,:,1);
% tempTrPr_Window = trial_Window_TrTr_TsPr(:,:,1);
% tempTrPo = trial_Window_TrPo(:,:,2);
% tempPoTr_Window = trial_Window_TrPo_TsTr(:,:,2);
% tempTrPo_Window = trial_Window_TrTr_TsPo(:,:,2);
% if selType == 1
%     tempTrPr = tempTrPr(2:end,2:end);
%     tempPrTr_Window = tempPrTr_Window(2:end,2:end);
%     tempTrPr_Window = tempTrPr_Window(2:end,2:end);
%     tempTrPo = tempTrPo(1:end-1,1:end-1);
%     tempPoTr_Window = tempPoTr_Window(1:end-1,1:end-1);
%     tempTrPo_Window = tempTrPo_Window(1:end-1,1:end-1);
% elseif selType == 2
%     tempTrPr = tempTrPr(2:end-1,2:end-1);
%     tempPrTr_Window = tempPrTr_Window(2:end-1,2:end-1);
%     tempTrPr_Window = tempTrPr_Window(2:end-1,2:end-1);
%     tempTrPo = tempTrPo(2:end-1,2:end-1);
%     tempPoTr_Window = tempPoTr_Window(2:end-1,2:end-1);
%     tempTrPo_Window = tempTrPo_Window(2:end-1,2:end-1);
% end
% 
% trPr = cell2mat(tempTrPr(triu(ones(length(tempTrPr)),1) & tril(ones(length(tempTrPr)),1)));
% win_PrTr = cell2mat(tempPrTr_Window(triu(ones(length(tempPrTr_Window)),1) & tril(ones(length(tempPrTr_Window)),1)));
% win_TrPr = cell2mat(tempTrPr_Window(triu(ones(length(tempTrPr_Window)),1) & tril(ones(length(tempTrPr_Window)),1)));
% trPo = cell2mat(tempTrPo(triu(ones(length(tempTrPo)),1) & tril(ones(length(tempTrPo)),1)));
% win_PoTr = cell2mat(tempPoTr_Window(triu(ones(length(tempPoTr_Window)),1) & tril(ones(length(tempPoTr_Window)),1)));
% win_TrPo = cell2mat(tempTrPo_Window(triu(ones(length(tempTrPo_Window)),1) & tril(ones(length(tempTrPo_Window)),1)));
% bar([1,2,3,4,5,6], [mean(trPr), mean(win_PrTr), mean(win_TrPr), mean(trPo), mean(win_PoTr), mean(win_TrPo)], 'k');
% hold on;
% swarmchart(ones(size(trPr)), trPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPr))+1, win_PrTr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPr))+2, win_TrPr, 5, 'markerfacecolor', [191/255, 191/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(trPo))+3, trPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+4, win_PoTr, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(win_TrPo))+5, win_TrPo, 5, 'markerfacecolor', [196/255, 33/255, 255/255], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% set(gca, 'xticklabel', [{'TrPr'}, {'Pr->Tr'}, {'Tr->Pr'}, {'TrPo'}, {'Po->Tr'}, {'Tr->Po'}], 'xticklabelrotation', 45);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% ITI Generalization Indices
% tempPoPr_Window = trial_Window_TrPo_TsPr(:,:,1);
% tempPrPo_Window = trial_Window_TrPr_TsPo(:,:,1);
% if selType == 1
%     tempPoPr_Window = tempPoPr_Window(2:end,2:end);
%     tempPrPo_Window = tempPrPo_Window(2:end,2:end);
% elseif selType == 2
%     tempPoPr_Window = tempPoPr_Window(2:end-1,2:end-1);
%     tempPrPo_Window = tempPrPo_Window(2:end-1,2:end-1);
% end
% popr_Prev = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),-1) & tril(ones(length(tempPoPr_Window)),-1)));
% popr_Curr = cell2mat(tempPoPr_Window(logical(eye(length(tempPoPr_Window)))));
% popr_Next = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),1) & tril(ones(length(tempPoPr_Window)),1)));
% prpo_Prev = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),-1) & tril(ones(length(tempPrPo_Window)),-1)));
% prpo_Curr = cell2mat(tempPrPo_Window(logical(eye(length(tempPrPo_Window)))));
% prpo_Next = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),1) & tril(ones(length(tempPrPo_Window)),1)));
% subplot(4,3,10)
% bar([1,2,3, 5,6,7], [mean(popr_Prev), mean(popr_Curr), mean(popr_Next), mean(prpo_Prev), mean(prpo_Curr), mean(prpo_Next)], 'k');
% hold on;
% swarmchart(ones(size(popr_Prev)), popr_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(popr_Curr))+1, popr_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(popr_Next))+2, popr_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Prev))+4, prpo_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Curr))+5, prpo_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Next))+6, prpo_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% set(gca, 'xtick', [1:3,5:7], 'xticklabel', [{'Po->Pr:Prev'}, {'Po->Pr:Curr'}, {'Po->Pr:Next'}, {'Pr->Po:Prev'}, {'Pr->Po:Curr'}, {'Pr->Po:Next'}], 'xticklabelrotation', 45);
% 
% tempPoPr_Window = trial_Window_TrPo_TsPr(:,:,2);
% tempPrPo_Window = trial_Window_TrPr_TsPo(:,:,2);
% if selType == 1
%     tempPoPr_Window = tempPoPr_Window(1:end-1,1:end-1);
%     tempPrPo_Window = tempPrPo_Window(1:end-1,1:end-1);
% elseif selType == 2
%     tempPoPr_Window = tempPoPr_Window(2:end-1,2:end-1);
%     tempPrPo_Window = tempPrPo_Window(2:end-1,2:end-1);
% end
% popr_Prev = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),-1) & tril(ones(length(tempPoPr_Window)),-1)));
% popr_Curr = cell2mat(tempPoPr_Window(logical(eye(length(tempPoPr_Window)))));
% popr_Next = cell2mat(tempPoPr_Window(triu(ones(length(tempPoPr_Window)),1) & tril(ones(length(tempPoPr_Window)),1)));
% prpo_Prev = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),-1) & tril(ones(length(tempPrPo_Window)),-1)));
% prpo_Curr = cell2mat(tempPrPo_Window(logical(eye(length(tempPrPo_Window)))));
% prpo_Next = cell2mat(tempPrPo_Window(triu(ones(length(tempPrPo_Window)),1) & tril(ones(length(tempPrPo_Window)),1)));
% subplot(4,3,11)
% bar([1,2,3, 5,6,7], [mean(popr_Prev), mean(popr_Curr), mean(popr_Next), mean(prpo_Prev), mean(prpo_Curr), mean(prpo_Next)], 'k');
% hold on;
% swarmchart(ones(size(popr_Prev)), popr_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(popr_Curr))+1, popr_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(popr_Next))+2, popr_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Prev))+4, prpo_Prev, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Curr))+5, prpo_Curr, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% swarmchart(ones(size(prpo_Next))+6, prpo_Next, 5, 'markerfacecolor', [191/255, 96/255, 0], 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% set(gca, 'xtick', [1:3,5:7], 'xticklabel', [{'Po->Pr:Prev'}, {'Po->Pr:Curr'}, {'Po->Pr:Next'}, {'Pr->Po:Prev'}, {'Pr->Po:Curr'}, {'Pr->Po:Next'}], 'xticklabelrotation', 45);
% colormap(cMap)
% 
% if selType == 1
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'Select Positions Based on Windows',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% elseif selType == 2
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'Common Positions Across Alignments (pos 2&3)',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% else
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'All Trials',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% end
% 
% tempPoPr_PI = trial_Window_PoPr(:,:,1);
% % tempPoPr_PI = tempPoPr_PI(2:end,2:end);
% % tempPoPr_PI = tempPoPr_PI(2:end-1,2:end-1);
% tempPoPr_PO = trial_Window_PoPr(:,:,2);
% % tempPoPr_PO = tempPoPr_PO(1:end-1,1:end-1);
% % tempPoPr_PO = tempPoPr_PO(2:end-1,2:end-1);
% 
% 
% 
% 
% % trial_Window_TrPo_TsPr
% % trial_Window_TrPr_TsPo
% % trial_Window_PoPr
% %% Plot Fits
% % Plot poke in aligned D data
% figure;
% subplot(2,3,1)
% tempTMat = trialD(:,:,1);
% if selType == 1
%     tempTMat = tempTMat(2:end,2:end);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeIn')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiPoLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPiRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiPoLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPiRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% title('Poke In Aligned');
% % Plot poke in aligned fits
% tempTrialFits = trial_TrialPersFit(:,:,1);
% tempITIfits = trial_ITIpersFit(:,:,1);
% if selType == 1
%     tempTrialFits = tempTrialFits(2:end,2:end);
%     tempITIfits = tempITIfits(2:end,2:end);
% elseif selType == 2
%     tempTrialFits = tempTrialFits(2:end-1,2:end-1);
%     tempITIfits = tempITIfits(2:end-1,2:end-1);
% end
% tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
% tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
% [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
% minTrialLat = minTrialLat./(1/dsRate);
% [minITIfit,minITIlat] = min(tempITIfits,[],2);
% minITIlat = minITIlat./(1/dsRate);
% subplot(2,3,2)
% swarmchart(ones(size(minTrialLat)), minTrialLat, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% hold on;
% swarmchart(ones(size(minITIlat))+1, minITIlat, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% [~,p,~,stats] = ttest(minTrialLat, minITIlat);
% if p<0.05
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
% else
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
% end
% set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
% ylabel('Latency to Best Fit (ms)');
% subplot(2,3,3)
% swarmchart(ones(size(minTrialFit)), minTrialFit, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% hold on;
% swarmchart(ones(size(minITIfit))+1, minITIfit, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% [~,p,~,stats] = ttest(minTrialFit, minITIfit);
% if p<0.05
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
% else
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
% end
% set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
% ylabel('Best Fit (Cos Sim)');
% 
% % Plot poke out aligned D data
% subplot(2,3,4)
% tempTMat = trialD(:,:,2);
% if selType == 1
%     tempTMat = tempTMat(1:end-1,1:end-1);
% elseif selType == 2
%     tempTMat = tempTMat(2:end-1,2:end-1);
% end
% mlb{end}.PlotTrialPDM(cell2mat(permute(tempTMat(logical(eye(length(tempTMat)))), [2,3,1])), 'rotate', 'clim', [-1.5 1.5], 'x', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'y', obsvTimeVect{strcmp(alignments, 'PokeOut')}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
% hold on;
% plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoPiLat, [1,2]), '--k','linewidth', 2);
% plot(get(gca, 'xlim'),repmat(grpPoRwdLat, [1,2]), ':k','linewidth', 2);
% plot(zeros(1,2),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoPiLat, [1,2]),get(gca, 'ylim'), '--k','linewidth', 2);
% plot(repmat(grpPoRwdLat, [1,2]),get(gca, 'ylim'), ':k','linewidth', 2);
% title('Poke Out Aligned');
% % Plot poke out aligned fits
% tempTrialFits = trial_TrialPersFit(:,:,2);
% tempITIfits = trial_ITIpersFit(:,:,2);
% if selType == 1
%     tempTrialFits = tempTrialFits(1:end-1,1:end-1);
%     tempITIfits = tempITIfits(1:end-1,1:end-1);
% elseif selType == 2
%     tempTrialFits = tempTrialFits(2:end-1,2:end-1);
%     tempITIfits = tempITIfits(2:end-1,2:end-1);
% end
% tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
% tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
% [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
% minTrialLat = minTrialLat./(1/dsRate);
% [minITIfit,minITIlat] = min(tempITIfits,[],2);
% minITIlat = minITIlat./(1/dsRate);
% subplot(2,3,5)
% swarmchart(ones(size(minTrialLat)), minTrialLat, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% hold on;
% swarmchart(ones(size(minITIlat))+1, minITIlat, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% [~,p,~,stats] = ttest(minTrialLat, minITIlat);
% if p<0.05
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
% else
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
% end
% set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
% ylabel('Latency to Best Fit (ms)');
% subplot(2,3,6)
% swarmchart(ones(size(minTrialFit)), minTrialFit, 5, 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% hold on;
% swarmchart(ones(size(minITIfit))+1, minITIfit, 5, 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markerfacealpha', 0.5);
% [~,p,~,stats] = ttest(minTrialFit, minITIfit);
% if p<0.05
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p), 'color','r');
% else
%     title(sprintf('t(%.00f) = %.02f; p = %.02i', stats.df, stats.tstat, p));
% end
% set(gca, 'xtick', 1:2, 'xticklabel', [{'Trial'}, {'ITI'}]);
% ylabel('Best Fit (Cos Sim)');
% 
% colormap(cMap)
% 
% 
% if selType == 1
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'Select Positions Based on Windows',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% elseif selType == 2
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'Common Positions Across Alignments (pos 2&3)',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% else
%     annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'All Trials',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% end
% 
% %% Holy shit we have a beta effect!
% for al = 1:length(alignments)
%     figure;    
%     betaLog = true(1,4);
%     tempLog = logical(eye(mlb{ani}.seqLength));
%     decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
%     if selType == 1
% %         if al == 1
% %             tempLog(4,4) = false;
% %             betaLog(4) = false;
% %         elseif al == 2
%             tempLog(1,1) = false;
%             betaLog(1) = false;
% %         end
%     elseif selType == 2        
%         tempLog(1,1) = false;
%         tempLog(4,4) = false;
%         betaLog(1) = false;
%         betaLog(4) = false;
%     elseif selType == 3
%         tempLog((1:mlb{ani}.seqLength)~=pos,(1:mlb{ani}.seqLength)~=pos) = false;
%         betaLog((1:mlb{ani}.seqLength)~=pos) = false;
%     end
%         
%     decLog(:,:,al) = tempLog;
%     subplot(4,3,1);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr_TsTr(decLog)), 'Beta', 'PrTr', []);
%     subplot(4,3,2);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrTr_TsPr(decLog)), 'Beta', 'TrPr', []);
%     subplot(4,3,3);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'PrTrTrPr', []);
%     subplot(4,3,4);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrTr_TsPo(decLog)), 'Beta', 'TrPo', []);
%     subplot(4,3,5);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo_TsTr(decLog)), 'Beta', 'PoTr', []);
%     subplot(4,3,6);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'TrPoPoTr', []);
%     subplot(4,3,7);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPo_TsPr(decLog)), 'Beta', 'PoPr', []);
%     subplot(4,3,8);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_TrPr_TsPo(decLog)), 'Beta', 'PrPo', []);
%     subplot(4,3,9);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), cell2mat(trial_Window_PoPr(decLog)), 'Beta', 'PoPrPrPo', []);
%     
%     tempTrialFits = cell2mat(trial_TrialPersFit(decLog));
%     [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
%     minTrialLat = minTrialLat./(1/dsRate);
%     tempITIfits = cell2mat(trial_ITIpersFit(decLog));
%     [minITIfit,minITIlat] = min(tempITIfits,[],2);
%     minITIlat = minITIlat./(1/dsRate);
%     subplot(4,3,10);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), minTrialLat, 'Beta', 'Trial Pers', []);
%     subplot(4,3,11);
%     corrScatPlot(cell2mat(trial_MeanBeta(betaLog,al)), minITIlat, 'Beta', 'ITI Pers', []);
%     
%     if selType == 1
%         annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%             'String', sprintf('%s: Select Positions Based on Windows', alignments{al}),...
%             'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     elseif selType == 2
%         annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%             'String', sprintf('%s: Common Positions Across Alignments (pos 2&3)', alignments{al}),...
%             'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     else
%         annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%             'String', sprintf('%s: All Trials', alignments{al}),...
%             'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     end
% end
% % trial_TrialPersFit = cell(size(trlPersFit,1), size(trlPersFit,2), size(trlPersFit,4));
% % trial_ITIpersFit = cell(size(itiPersFit,1), size(itiPersFit,2), size(itiPersFit,4));
% 
% %% Pre-Trial vs Trial Port Entry Aligned
% figure; 
% tempLog = logical(eye(mlb{ani}.seqLength));
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(:,:,1) = tempLog;
% subplot(4,4,[1:3,5:7,9:11,13:15])
% corrScatPlot(cell2mat(trial_MeanBeta(:,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('All Positions');
% 
% subplot(4,4,4);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(1,1,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 1');
% 
% subplot(4,4,8);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(2,2,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 2');
% 
% subplot(4,4,12);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(3,3,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 3');
% 
% subplot(4,4,16);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(4,4,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 4');
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', 'Port Entry Aligned: Pre-Trial vs. Trial; Current Trial Info',...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% 
% %% Post-Trial vs Trial Port Entry Aligned
% figure; 
% tempLog = logical(eye(mlb{ani}.seqLength));
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(:,:,2) = tempLog;
% subplot(4,4,[1:3,5:7,9:11,13:15])
% corrScatPlot(cell2mat(trial_MeanBeta(:,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('All Positions');
% 
% subplot(4,4,4);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(1,1,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 1');
% 
% subplot(4,4,8);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(2,2,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 2');
% 
% subplot(4,4,12);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(3,3,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 3');
% 
% subplot(4,4,16);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(4,4,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Position 4');
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', 'Port Exit Aligned: Post-Trial vs. Trial; Current Trial Info',...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% 
% 
% %% Pre-Trial vs Trial Port Entry Aligned Previous Trial Info
% figure; 
% tempLog = triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(:,:,1) = tempLog;
% subplot(3,3,[1:2,4:5,7:8])
% corrScatPlot(cell2mat(trial_MeanBeta(2:end,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('All Positions');
% 
% subplot(3,3,3);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(2,1,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 1 in Pos 2');
% 
% subplot(3,3,6);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(3,2,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 2 in Pos 3');
% 
% subplot(3,3,9);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(4,3,1) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(4,al)), cell2mat(trial_Window_TrPr(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 3 in Pos 4');
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', 'Port Entry Aligned: Pre-Trial vs. Trial; Prev Trial Info',...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% 
% %% Pre-Trial vs Trial Port Entry Aligned Previous Trial Info
% figure; 
% tempLog = triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(:,:,2) = tempLog;
% subplot(3,3,[1:2,4:5,7:8])
% corrScatPlot(cell2mat(trial_MeanBeta(1:end-1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('All Positions');
% 
% subplot(3,3,3);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(1,2,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(1,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 1 in Pos 2');
% 
% subplot(3,3,6);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(2,3,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(2,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 2 in Pos 3');
% 
% subplot(3,3,9);
% decLog = false(mlb{ani}.seqLength, mlb{ani}.seqLength, length(alignments));
% decLog(3,4,2) = true;
% corrScatPlot(cell2mat(trial_MeanBeta(3,al)), cell2mat(trial_Window_TrPo(decLog)), 'Beta', 'Cross Epoch Generalization', []);
% title('Pos 3 in Pos 4');
% 
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%     'String', 'Port Entry Aligned: Post-Trial vs. Trial; Next Trial Info',...
%     'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

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
piWIn_Trl = trial_Window_TrTr_TsTr(:,:,1);
piWIn_ITD = trial_Window_TrITD_TsITD(:,:,1);
piX_PrTr = trial_Window_TrPr_TsTr(:,:,1);
piX_TrPr = trial_Window_TrTr_TsPr(:,:,1);
piWIn_Trl = piWIn_Trl(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
piWIn_ITD = piWIn_ITD(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
piX_PrTr = piX_PrTr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
piX_TrPr = piX_TrPr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));

spIDs = reshape(1:mlb{ani}.seqLength*(mlb{ani}.seqLength+1), [mlb{ani}.seqLength+1, mlb{ani}.seqLength])';
figure;
subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,reshape(spIDs(:,1:end-1), [1,mlb{ani}.seqLength^2])); 
mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(piWIn_Trl),1,0.05,'r');
hold on;
mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(piWIn_ITD),1,0.05,'b');
mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(piX_PrTr),1,0.05,[191/255, 191/255, 0]);
mlb{ani}.PlotMeanVarSwarmBar(4,cell2mat(piX_TrPr),1,0.05,[191/255, 191/255, 0]);
set(gca, 'xtick', 1:4, 'xticklabel', [{'Trial'}, {'ITD'}, {'Pre->Trial'}, {'Trial->Pre'}], 'xticklabelrotation', 45);
for sp = 1:mlb{ani}.seqLength
    subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,spIDs(sp,end));
    mlb{ani}.PlotMeanVarSwarmBar(1,piWIn_Trl{sp},1,0.05,'r');
    hold on;
    mlb{ani}.PlotMeanVarSwarmBar(2,piWIn_ITD{sp},1,0.05,'b');
    mlb{ani}.PlotMeanVarSwarmBar(3,piX_PrTr{sp},1,0.05,[191/255, 191/255, 0]);
    mlb{ani}.PlotMeanVarSwarmBar(4,piX_TrPr{sp},1,0.05,[191/255, 191/255, 0]);
    title(sp);
end
linkaxes
annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', 'Port Entry Aligned',...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% Port Exit Aligned
% poWIn_Trl = trial_Window_TrTr_TsTr(:,:,2);
% poWIn_ITD = trial_Window_TrITD_TsITD(:,:,2);
% poX_TrPo = trial_Window_TrTr_TsPo(:,:,2);
% poX_PoTr = trial_Window_TrPo_TsTr(:,:,2);
% poWIn_Trl = poWIn_Trl(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
% poWIn_ITD = poWIn_ITD(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
% poX_TrPo = poX_TrPo(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
% poX_PoTr = poX_PoTr(logical(eye(mlb{ani}.seqLength,mlb{ani}.seqLength)));
% 
% figure;
% subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,reshape(spIDs(:,1:end-1), [1,mlb{ani}.seqLength^2])); 
% mlb{ani}.PlotMeanVarSwarmBar(1,cell2mat(poWIn_Trl),1,0.05,'r');
% hold on;
% mlb{ani}.PlotMeanVarSwarmBar(2,cell2mat(poWIn_ITD),1,0.05,'b');
% mlb{ani}.PlotMeanVarSwarmBar(3,cell2mat(poX_PoTr),1,0.05,[191/255, 191/255, 0]);
% mlb{ani}.PlotMeanVarSwarmBar(4,cell2mat(poX_TrPo),1,0.05,[191/255, 191/255, 0]);
% set(gca, 'xtick', 1:4, 'xticklabel', [{'Trial'}, {'ITD'}, {'Post->Trial'}, {'Trial->Post'}], 'xticklabelrotation', 45);
% for sp = 1:mlb{ani}.seqLength
%     subplot(mlb{ani}.seqLength, mlb{ani}.seqLength+1,spIDs(sp,end));
%     mlb{ani}.PlotMeanVarSwarmBar(1,poWIn_Trl{sp},1,0.05,'r');
%     hold on;
%     mlb{ani}.PlotMeanVarSwarmBar(2,poWIn_ITD{sp},1,0.05,'b');
%     mlb{ani}.PlotMeanVarSwarmBar(3,poX_PoTr{sp},1,0.05,[191/255, 191/255, 0]);
%     mlb{ani}.PlotMeanVarSwarmBar(4,poX_TrPo{sp},1,0.05,[191/255, 191/255, 0]);
%     title(sp);
% end
% linkaxes
% annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
%         'String', 'Port Exit Aligned',...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
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
    tempITIfits = trial_ITIpersFit(:,:,al);
    tempTrialFits = cell2mat(tempTrialFits(logical(eye(length(tempTrialFits)))));
    tempITIfits = cell2mat(tempITIfits(logical(eye(length(tempITIfits)))));
    itiDur = find(sum(~isnan(tempITIfits),1)./size(tempITIfits,1)==1,1,'last');
    
    [minTrialFit,minTrialLat] = min(tempTrialFits,[],2);
    minTrialLat = minTrialLat./(1/dsRate);
    [minITIfit,minITIlat] = min(tempITIfits,[],2);
    minITIlat = minITIlat./(1/dsRate);
    
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
    histogram(minTrialLat, 0:dsRate*2:itiDur*dsRate,'FaceColor', 'r', 'FaceAlpha', 0.5);
    hold on;
    histogram(minITIlat, 0:dsRate*2:itiDur*dsRate,'FaceColor', 'k', 'FaceAlpha', 0.5);
    set(gca, 'xlim', [0 itiDur*dsRate]);
    
%     linkaxes(sp(1:2), 'y');
    linkaxes(sp([1,3]),'x');
end
%% Plot Lag Values
for al = 1:length(alignments)
    figure;
    tempTMat = trialD(:,:,al);
    
    subplot(1,3,1)
    nextPosD = tempTMat(triu(ones(length(tempTMat)),1) & tril(ones(length(tempTMat)),1));
    mlb{end}.PlotTrialPDM(cell2mat(permute(nextPosD(2:end), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
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
    title('Next Position');
    subplot(1,3,2)
    currPosD = tempTMat(logical(eye(length(tempTMat))));
    mlb{end}.PlotTrialPDM(cell2mat(permute(currPosD(2:end-1), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
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
    subplot(1,3,3)
    prevPosD = tempTMat(triu(ones(length(tempTMat)),-1) & tril(ones(length(tempTMat)),-1));
    mlb{end}.PlotTrialPDM(cell2mat(permute(prevPosD(1:end-1), [2,3,1])), 'rotate', 'clim', 'Max', 'x', obsvTimeVect{al}, 'y', obsvTimeVect{al}, 'xlabel', 'Test Time', 'ylabel', 'Train Time');
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
    colormap(cMap)
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', sprintf('Decodability: %s Alignment', alignments{al}),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
%% Save output
save('PFC_XTD_Windows_PosChance.mat', '-v7.3');