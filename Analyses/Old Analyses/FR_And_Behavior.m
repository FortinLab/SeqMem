fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
setupSeqLength = 4; 
% 
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% tets = [1,22,17,18,17];
% setupSeqLength = 5;
% 
% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4; 

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4; 

binSize = 200;
dsRate = 1;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
binType = 'gauss';
% binType = 'box';
sfpYN = 0;
% behavLatType = 'rel';
behavLatType = 'abs';

% frThresh = 0.2/(binSize/1000);
frThresh = 0.00001;
%% Data Vectors
% Neural
popVects = cell(length(fileDirs),1);
popVectsSortVect = cell(length(fileDirs),1);
popVectsThreshVect = cell(length(fileDirs),1);
% Behavior
fiscPokeOutLat = cell(length(fileDirs),setupSeqLength);
fiscRwdSigLat = cell(length(fileDirs),setupSeqLength);
fiscRwdDelivLat = cell(length(fileDirs),setupSeqLength);
itiDur = cell(length(fileDirs),setupSeqLength);
pokeLatRaw = cell(length(fileDirs),2,2);
pokeLatOPrawC = cell(length(fileDirs),setupSeqLength,2);
pokeLatOPrawI = cell(length(fileDirs),setupSeqLength,2);
smi = nan(length(fileDirs),1);
dPrm = nan(length(fileDirs),1);
ri = nan(length(fileDirs),1);
smiByOP = nan(length(fileDirs),setupSeqLength,2);
dPrmByOP = nan(length(fileDirs),setupSeqLength,2);
riByOP = nan(length(fileDirs),setupSeqLength,2);
tMatAcc = nan(setupSeqLength,setupSeqLength,length(fileDirs));
lagAcc = nan(length(fileDirs),setupSeqLength*2-1);
tMatLatC = cell(setupSeqLength,setupSeqLength,length(fileDirs));
tMatLatI = cell(setupSeqLength,setupSeqLength,length(fileDirs));
lagLat = cell(length(fileDirs),setupSeqLength*2-1,2);
roc = nan(length(fileDirs),2);
rocOPpos = nan(length(fileDirs),setupSeqLength,2);
rocOPodr = nan(length(fileDirs),setupSeqLength,2);
taoAcc = nan(length(fileDirs),1);
taoLat = cell(length(fileDirs),2);
betaModCells = nan(length(fileDirs),1);
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
% %     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.binType = binType;
    %% Pull out spiking data
    [trlSpikes, timeVect] = mlb.PP_TrialMatrix_Spiking(trlWindow{1}, alignment{1});
    nsmblSpks = nan(size(trlSpikes,2),size(trlSpikes,1), mlb.seqLength);
    peakNdx = nan(size(trlSpikes,2),mlb.seqLength);
    peakRt = nan(size(trlSpikes,2),mlb.seqLength);
    for op = 1:mlb.seqLength
        tempNsmbl = mean(trlSpikes(:,:,mlb.fiscTrials(op,:)),3);
        for uni = 1:size(tempNsmbl,2)
            peakNdx(uni,op) = find(tempNsmbl(:,uni)==max(tempNsmbl(:,uni)),1,'first');
            peakRt(uni,op) = max(tempNsmbl(:,uni));
            tempNsmbl(:,uni) = tempNsmbl(:,uni)./max(tempNsmbl(:,uni));
%             tempNsmbl(:,uni) = tempNsmbl(:,uni);
        end
        nsmblSpks(:,:,op) = tempNsmbl';
    end
    popVects{ani} = nsmblSpks;
    popVectsSortVect{ani} = peakNdx;
    popVectsThreshVect{ani} = peakRt;
    %% Extract behavioral variables    
    if strcmp(behavLatType, 'rel')
        relPokeDur = arrayfun(@(a)a.PokeDuration - a.TargetDuration, mlb.trialInfo);
        for trl = 1:length(mlb.trialInfo)
            mlb.trialInfo(trl).PokeDuration = relPokeDur(trl);
        end
        mlb.SummarizeSessionBehavior;
    end
    tmLatC = mlb.transMatLatRaw(:,:,1);
    tmLatIC = mlb.transMatLatRaw(:,:,2);
    for op = 1:mlb.seqLength
        fiscPokeOutLat{ani,op} = ([mlb.trialInfo(mlb.fiscTrials(op,:)).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials(op,:)).PokeInIndex])'/1000;
        fiscRwdSigLat{ani,op} = ([mlb.trialInfo(mlb.fiscTrials(op,:)).RewardSignalIndex] - [mlb.trialInfo(mlb.fiscTrials(op,:)).PokeInIndex])'/1000;
        fiscRwdDelivLat{ani,op} = ([mlb.trialInfo(mlb.fiscTrials(op,:)).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials(op,:)).PokeInIndex])'/1000;
        if op ~= 1
            itiDur{ani,op} = ([mlb.trialInfo(mlb.fiscTrials(op,:)).PokeInIndex] - [mlb.trialInfo(mlb.fiscTrials(op,:)-1).PokeOutIndex])'/1000;
        end
        smiByOP(ani,:,1) = mlb.smiByPos;
        smiByOP(ani,:,2) = mlb.smiByOdr;
        dPrmByOP(ani,:,1) = mlb.dPrimeByPos;
        dPrmByOP(ani,:,2) = mlb.dPrimeByOdr;
        riByOP(ani,:,1) = mlb.riByPos;
        riByOP(ani,:,2) = mlb.riByOdr;
        rocOPpos(ani,op,:) = [mlb.responseMatrixByPos{op}(1,1)/sum(mlb.responseMatrixByPos{op}(1,:)), mlb.responseMatrixByPos{op}(2,1)/sum(mlb.responseMatrixByPos{op}(2,:))];
        rocOPodr(ani,op,:) = [mlb.responseMatrixByOdr{op}(1,1)/sum(mlb.responseMatrixByOdr{op}(1,:)), mlb.responseMatrixByOdr{op}(2,1)/sum(mlb.responseMatrixByOdr{op}(2,:))];
    end
    tMatAcc(:,:,ani) = mlb.transMatAcc;
    tMatLatC(:,:,ani) = mlb.transMatLatRaw(:,:,1);
    tMatLatI(:,:,ani) = mlb.transMatLatRaw(:,:,2);
    taoAcc(ani) = mlb.taoAcc;
    taoLat(ani,:) = mlb.taoLatRaw;
    
    if sfpYN == 0
        pokeLatRaw{ani,1,1} = cell2mat(tmLatC(logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,2,1} = cell2mat(tmLatIC(logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,1,2} = cell2mat(tmLatC(~logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,2,2} = cell2mat(tmLatIC(~logical(eye(mlb.seqLength))));
        smi(ani) = mlb.smi;
        dPrm(ani) = mlb.dPrime;
        ri(ani) = mlb.ri;
        lagAcc(ani,:) = mlb.lagAccVect(1,:,1)./sum(mlb.lagAccVect,3);
        lagLat(ani,:,1) = mlb.lagLatVectRaw(:,:,1);
        lagLat(ani,:,2) = mlb.lagLatVectRaw(:,:,2);
        roc(ani,:) = [mlb.responseMatrix(1,1)/sum(mlb.responseMatrix(1,:)), mlb.responseMatrix(2,1)/sum(mlb.responseMatrix(2,:))];
    else
        tmLatC(1,:) = cell(1,mlb.seqLength);
        tmLatIC(1,:) = cell(1,mlb.seqLength);
        pokeLatRaw{ani,1,1} = cell2mat(tmLatC(logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,2,1} = cell2mat(tmLatIC(logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,1,2} = cell2mat(tmLatC(~logical(eye(mlb.seqLength))));
        pokeLatRaw{ani,2,2} = cell2mat(tmLatIC(~logical(eye(mlb.seqLength))));
        smi(ani) = mlb.smiSFP;
        dPrm(ani) = mlb.dPrimeSFP;
        ri(ani) = mlb.riSFP;
        lagAcc(ani,:) = mlb.lagAccVectSFP(1,:,1)./sum(mlb.lagAccVectSFP,3);
        lagLat(ani,:,1) = mlb.lagLatVectSFPraw(:,:,1);
        lagLat(ani,:,2) = mlb.lagLatVectSFPraw(:,:,2);
        roc(ani,:) = [mlb.responseMatrixSFP(1,1)/sum(mlb.responseMatrixSFP(1,:)), mlb.responseMatrixSFP(2,1)/sum(mlb.responseMatrixSFP(2,:))];
    end
end

%% Plot Behavior
%%%%% Bandaid to mask an error in the poke duration code for one trial... probably more affected... but until I get the curation code up this is how we'll process things
tMatLatC = cellfun(@(a){a(a<10)}, tMatLatC);
tMatLatI = cellfun(@(a){a(a<10)}, tMatLatI);
lagLat = cellfun(@(a){a(a<10)}, lagLat);
pokeLatRaw = cellfun(@(a){a(a<10)}, pokeLatRaw);
% Behavior Summary 1
figure;
% Session ROC 
rocMean_Plot = subplot(6,3,[1,4]);
scatter(roc(:,2), roc(:,1), 20, 'ko');
hold on;
errorbar(mean(roc(:,2)), mean(roc(:,1)),...
    std(roc(:,1)), std(roc(:,1)),...
    std(roc(:,2)), std(roc(:,2)), 'o', 'linewidth', 2, 'color', 'k', 'MarkerFaceColor', 'k');
plot([0 1], [0 1], '--k');
set(gca, 'xlim', [0 1], 'ylim', [0 1]);
xlabel('False Alarm Rate');
ylabel('Hit Rate');
text(mean(roc(:,2))+0.05, mean(roc(:,1))+0.05, sprintf('d''=%.02f +/- %.02fStD', mean(dPrm), std(dPrm)));
% Session SMI
smiMean_Bar = subplot(6,3,7);
bar(mean(smi), 'facecolor', 'none');
hold on;
scatter(normrnd(1,0.05, size(smi)), smi, 20, 'ko');
errorbar(mean(smi), std(smi), 'color', 'k', 'capsize', 0);
title(sprintf('SMI = %.02f +/- %.02fStD', mean(smi), std(smi)));
set(gca, 'ylim', [0 1]);
ylabel('SMI');
% Session D-Prime
dPrmMean_Bar = subplot(6,3,10);
bar(mean(dPrm), 'facecolor', 'none');
hold on;
scatter(normrnd(1,0.05, size(dPrm)), dPrm, 20, 'ko');
errorbar(mean(dPrm), std(dPrm), 'color', 'k', 'capsize', 0);
title(sprintf('d''=%.02f +/- %.02fStD', mean(dPrm), std(dPrm)));
set(gca, 'ylim', [0 5]);
ylabel('d''');
% Session RI
riMean_Bar = subplot(6,3,13);
bar(mean(ri), 'facecolor', 'none');
hold on;
scatter(normrnd(1,0.05, size(ri)), ri, 20, 'ko');
errorbar(mean(ri), std(ri), 'color', 'k', 'capsize', 0);
title(sprintf('d''=%.02f +/- %.02fStD', mean(ri), std(ri)));
set(gca, 'ylim', [-1 1]);
ylabel('RI');
% Session Poke Latencies
lat_Bar = subplot(6,3,16);
aniPokeLat = cellfun(@(a)mean(a, 'omitnan'), pokeLatRaw);
aniPokeLat = [aniPokeLat(:,:,1), aniPokeLat(:,:,2)];
pokesISC = cell2mat(pokeLatRaw(:,1,1));
pokesISI = cell2mat(pokeLatRaw(:,2,1));
pokesOSC = cell2mat(pokeLatRaw(:,1,2));
pokesOSI = cell2mat(pokeLatRaw(:,2,2));
bar([mean(pokesISC), mean(pokesISI), mean(pokesOSC), mean(pokesOSI)], 'facecolor', 'none');
hold on;
scatter(lat_Bar,normrnd(1,0.05, size(pokesISC)), pokesISC, 20, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2);
scatter(lat_Bar,normrnd(2,0.05, size(pokesISI)), pokesISI, 20, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2);
scatter(lat_Bar,normrnd(3,0.05, size(pokesOSC)), pokesOSC, 20, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2);
scatter(lat_Bar,normrnd(4,0.05, size(pokesOSI)), pokesOSI, 20, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2);
for ani = 1:size(aniPokeLat,1)
    plot(aniPokeLat(ani,:), 'marker', 'o', 'MarkerEdgeColor', 'k', 'color', 'k', 'MarkerFaceColor', 'k');
end
if strcmp(behavLatType,'rel')
    set(gca, 'ylim', [-1.5 1.5], 'xtick', 1:4, 'xticklabel', [{'ISC'}, {'ISI'}, {'OSC'}, {'OSI'}], 'xticklabelrotation', 45);
else
    set(gca, 'ylim', [0 2], 'xtick', 1:4, 'xticklabel', [{'ISC'}, {'ISI'}, {'OSC'}, {'OSI'}], 'xticklabelrotation', 45);
end

% OP ROC
rocOdr_Plot = subplot(6,3,[2,5]);
plot(rocOdr_Plot, [0 1], [0 1], '--k');
hold(rocOdr_Plot, 'on');
rocPos_Plot = subplot(6,3,[3,6]);
plot(rocPos_Plot, [0 1], [0 1], '--k');
hold(rocPos_Plot, 'on');
% OP SMI
smiOdr_Bar = subplot(6,3,8);
hold(smiOdr_Bar, 'on');
smiPos_Bar = subplot(6,3,9);
hold(smiPos_Bar, 'on');
% OP D-Prime
dPrmOdr_Bar = subplot(6,3,11);
hold(dPrmOdr_Bar, 'on');
dPrmPos_Bar = subplot(6,3,12);
hold(dPrmPos_Bar, 'on');
% OP RI
riOdr_Bar = subplot(6,3,14);
hold(riOdr_Bar, 'on');
riPos_Bar = subplot(6,3,15);
hold(riPos_Bar, 'on');
% OP Latencies
latOdr_Bar = subplot(6,3,17);
hold(latOdr_Bar, 'on');
latPos_Bar = subplot(6,3,18);
hold(latPos_Bar, 'on');
for op = 1:4
    % OP ROC
    scatter(rocOdr_Plot, rocOPodr(:,op,2), rocOPodr(:,op,1), 20,...
        'markerfacecolor', 'none', 'markeredgecolor', mlb.PositionColors(op,:), 'marker', 'o');    
    errorbar(rocOdr_Plot, mean(rocOPodr(:,op,2)), mean(rocOPodr(:,op,1)),...
        std(rocOPodr(:,op,1)), std(rocOPodr(:,op,1)),...
        std(rocOPodr(:,op,2)), std(rocOPodr(:,op,2)), 'o', 'linewidth', 2, 'color', mlb.PositionColors(op,:), 'MarkerFaceColor', mlb.PositionColors(op,:));
    if op~=1
        scatter(rocPos_Plot, rocOPpos(:,op,2), rocOPpos(:,op,1), 20,...
            'markerfacecolor', 'none', 'markeredgecolor', mlb.PositionColors(op,:), 'marker', 'o');
        errorbar(rocPos_Plot, mean(rocOPpos(:,op,2)), mean(rocOPpos(:,op,1)),...
            std(rocOPpos(:,op,1)), std(rocOPpos(:,op,1)),...
            std(rocOPpos(:,op,2)), std(rocOPpos(:,op,2)), 'o', 'linewidth', 2, 'color', mlb.PositionColors(op,:), 'MarkerFaceColor', mlb.PositionColors(op,:));
    end
    % OP SMI
    bar(smiOdr_Bar, op, mean(smiByOP(:,op,2)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(smiOdr_Bar, ones([size(smiByOP,1),1]).*op, smiByOP(:,op,2), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
    errorbar(smiOdr_Bar, op, mean(smiByOP(:,op,2)), std(smiByOP(:,op,2)), 'color', 'k', 'capsize', 0);
    if op~=1        
        bar(smiPos_Bar, op, mean(smiByOP(:,op,1)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
        scatter(smiPos_Bar, ones([size(smiByOP,1),1]).*op, smiByOP(:,op,1), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
        errorbar(smiPos_Bar, op, mean(smiByOP(:,op,1)), std(smiByOP(:,op,1)), 'color', 'k', 'capsize', 0);        
    end
    % OP D-Prime
    bar(dPrmOdr_Bar, op, mean(dPrmByOP(:,op,2)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(dPrmOdr_Bar, ones([size(dPrmByOP,1),1]).*op, dPrmByOP(:,op,2), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
    errorbar(dPrmOdr_Bar, op, mean(dPrmByOP(:,op,2)), std(dPrmByOP(:,op,2)), 'color', 'k', 'capsize', 0);
    if op~=1
        bar(dPrmPos_Bar, op, mean(dPrmByOP(:,op,1)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
        scatter(dPrmPos_Bar, ones([size(dPrmByOP,1),1]).*op, dPrmByOP(:,op,1), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
        errorbar(dPrmPos_Bar, op, mean(dPrmByOP(:,op,1)), std(dPrmByOP(:,op,1)), 'color', 'k', 'capsize', 0);
    end
    % OP RI
    bar(riOdr_Bar, op, mean(riByOP(:,op,2)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(riOdr_Bar, ones([size(riByOP,1),1]).*op, riByOP(:,op,2), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
    errorbar(riOdr_Bar, op, mean(riByOP(:,op,2)), std(riByOP(:,op,2)), 'color', 'k', 'capsize', 0);
    if op~=1
        bar(riPos_Bar, op, mean(riByOP(:,op,1)), 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
        scatter(riPos_Bar, ones([size(riByOP,1),1]).*op, riByOP(:,op,1), 20, 'markerfacecolor', mlb.PositionColors(op,:), 'markeredgecolor', 'k', 'marker', 'o');
        errorbar(riPos_Bar, op, mean(riByOP(:,op,1)), std(riByOP(:,op,1)), 'color', 'k', 'capsize', 0);
    end
    % OP Latencies
    iscAni = permute(cellfun(@(a)mean(a), tMatLatC(op,op,:)), [3,1,2]);
    isiAni = permute(cellfun(@(a)mean(a), tMatLatI(op,op,:)), [3,1,2]);
    osLog = true(1,mlb.seqLength);
    osLog(op) = false;
    tempOSCani = permute(tMatLatC(op,osLog,:), [2,1,3]);
    tempOSIani = permute(tMatLatI(op,osLog,:), [2,1,3]);
    oscAni = nan(size(tempOSCani,1),1);
    osiAni = nan(size(tempOSIani,1),1);
    for a = 1:size(tempOSCani,3)
        oscAni(a) = mean(cell2mat(tempOSCani(:,:,a)));
        osiAni(a) = mean(cell2mat(tempOSIani(:,:,a)));
    end

    tempISClat = cell2mat(reshape(tMatLatC(op,op,:), [size(tMatLatC,3),1]));
    bar(latOdr_Bar, op-0.3, mean(tempISClat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latOdr_Bar, normrnd(op-0.3,0.01, size(tempISClat)), tempISClat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    
    scatter(latOdr_Bar, normrnd(op-0.3,0, size(iscAni)), iscAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    tempISIlat = cell2mat(reshape(tMatLatI(op,op,:), [size(tMatLatI,3),1]));
    bar(latOdr_Bar, op-0.1, mean(tempISIlat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latOdr_Bar, normrnd(op-0.1,0.01, size(tempISIlat)), tempISIlat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    scatter(latOdr_Bar, normrnd(op-0.1,0, size(isiAni)), isiAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    tempOSClat = cell2mat(reshape(tMatLatC(op,osLog,:), [size(tMatLatC,3)*sum(osLog),1]));
    bar(latOdr_Bar, op+0.1, mean(tempOSClat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latOdr_Bar, normrnd(op+0.1,0.01, size(tempOSClat)), tempOSClat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    scatter(latOdr_Bar, normrnd(op+0.1,0, size(oscAni)), oscAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    tempOSIlat = cell2mat(reshape(tMatLatI(op,osLog,:), [size(tMatLatI,3)*sum(osLog),1]));
    bar(latOdr_Bar, op+0.3, mean(tempOSIlat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latOdr_Bar, normrnd(op+0.3,0.01, size(tempOSIlat)), tempOSIlat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    scatter(latOdr_Bar, normrnd(op+0.3,0, size(osiAni)), osiAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    plot(latOdr_Bar, repmat(op-0.3:0.2:op+0.3,[length(fileDirs),1])', [iscAni, isiAni, oscAni, osiAni]', 'k');
    
    bar(latPos_Bar, op-0.3, mean(tempISClat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latPos_Bar, normrnd(op-0.3,0.01, size(tempISClat)), tempISClat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    scatter(latPos_Bar, normrnd(op-0.3,0, size(iscAni)), iscAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    bar(latPos_Bar, op-0.1, mean(tempISIlat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
    scatter(latPos_Bar, normrnd(op-0.1,0.01, size(tempISIlat)), tempISIlat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2);
    scatter(latPos_Bar, normrnd(op-0.1,0, size(isiAni)), isiAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    if op~=1
        osLog = true(mlb.seqLength,1);
        osLog(op) = false;
        tempOSCani = tMatLatC(osLog,op,:);
        tempOSIani = tMatLatI(osLog,op,:);
        oscAni = nan(size(tempOSCani,1),1);
        osiAni = nan(size(tempOSIani,1),1);
        for a = 1:size(tempOSCani,3)
            oscAni(a) = mean(cell2mat(tempOSCani(:,:,a)));
            osiAni(a) = mean(cell2mat(tempOSIani(:,:,a)));
        end
        tempOSClat = cell2mat(reshape(tMatLatC(osLog,op,:), [size(tMatLatC,3)*sum(osLog),1]));
        bar(latPos_Bar, op+0.1, mean(tempOSClat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
        scatter(latPos_Bar, normrnd(op+0.1,0.01, size(tempOSClat)), tempOSClat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2)
        scatter(latPos_Bar, normrnd(op+0.1,0, size(oscAni)), oscAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        tempOSIlat = cell2mat(reshape(tMatLatI(osLog,op,:), [size(tMatLatI,3)*sum(osLog),1]));
        bar(latPos_Bar, op+0.3, mean(tempOSIlat), 0.2, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(op,:));
        scatter(latPos_Bar, normrnd(op+0.3,0.01, size(tempOSIlat)), tempOSIlat, 20, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceColor', mlb.PositionColors(op,:), 'MarkerFaceAlpha', 0.2)
        scatter(latPos_Bar, normrnd(op+0.3,0, size(osiAni)), osiAni, 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        plot(latPos_Bar, repmat(op-0.3:0.2:op+0.3,[length(fileDirs),1])', [iscAni, isiAni, oscAni, osiAni]', 'k');
    end
end
plot(smiOdr_Bar, smiByOP(:,:,2)', 'k');
plot(smiPos_Bar, smiByOP(:,:,1)', 'k');
plot(dPrmOdr_Bar, dPrmByOP(:,:,2)', 'k');
plot(dPrmPos_Bar, dPrmByOP(:,:,1)', 'k');
plot(riOdr_Bar, riByOP(:,:,2)', 'k');
plot(riPos_Bar, riByOP(:,:,1)', 'k');
title(rocOdr_Plot, 'ROC by Odor');
title(rocPos_Plot, 'ROC by Position');
xlabel([rocOdr_Plot,rocPos_Plot], 'False Alarm Rate');
ylabel([rocOdr_Plot,rocPos_Plot], 'Hit Rate');
title([smiOdr_Bar,smiPos_Bar], 'SMI');
ylabel([smiOdr_Bar,smiPos_Bar], 'SMI');
title([dPrmOdr_Bar,dPrmPos_Bar], 'd''');
ylabel([dPrmOdr_Bar,dPrmPos_Bar], 'd''');
title([riOdr_Bar,riPos_Bar], 'RI');
ylabel([riOdr_Bar,riPos_Bar], 'RI');
title([lat_Bar,latOdr_Bar, latPos_Bar], 'Poke Duration');
ylabel([lat_Bar,latOdr_Bar, latPos_Bar], 'Duration (ms)');
set([rocOdr_Plot,rocPos_Plot], 'xlim', [0 1], 'ylim', [0 1]);
set([smiOdr_Bar, smiPos_Bar], 'xlim', [0 mlb.seqLength+1], 'xtick', 1:mlb.seqLength, 'ylim', [0 1]);
set([dPrmOdr_Bar, dPrmPos_Bar], 'xlim', [0 mlb.seqLength+1], 'xtick', 1:mlb.seqLength, 'ylim', [0 5]);
set([riOdr_Bar, riPos_Bar], 'xlim', [0 mlb.seqLength+1], 'xtick', 1:mlb.seqLength, 'ylim', [-1 1]);
if strcmp(behavLatType,'rel')
    set([latOdr_Bar, latPos_Bar], 'xlim', [0 mlb.seqLength+1], 'xtick', 1:mlb.seqLength, 'ylim', [-1.5 1.5]);
else
    set([latOdr_Bar, latPos_Bar], 'xlim', [0 mlb.seqLength+1], 'xtick', 1:mlb.seqLength, 'ylim', [0 2.5]);
end
set([smiOdr_Bar, dPrmOdr_Bar, riOdr_Bar, latOdr_Bar], 'xticklabel', mlb.Rosetta(1:mlb.seqLength))
xlabel([smiOdr_Bar, dPrmOdr_Bar, riOdr_Bar, latOdr_Bar], 'Odor')
xlabel([smiPos_Bar, dPrmPos_Bar, riPos_Bar, latPos_Bar], 'Position')

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 1: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Behavior Summary 2
figure; 
subplot(1,2,1)
bar(mean(taoAcc), 'facecolor', 'none', 'edgecolor', 'k');
hold on;
scatter(normrnd(1,0.1,size(taoAcc)), taoAcc, 20, 'markeredgecolor', 'none', 'markerfacecolor','k');
errorbar(mean(taoAcc), std(taoAcc), 'capsize', 0, 'color', 'k');
title('TAO Accuracy (Mean +/- StD)');
ylabel('Accuracy (% correct trials)');

subplot(1,2,2)
bar(1,mean(cell2mat(taoLat(:,1))), 'facecolor', 'none', 'edgecolor', 'k');
hold on;
scatter(normrnd(1,0.1, size(cell2mat(taoLat(:,1)))), cell2mat(taoLat(:,1)), 20, 'markeredgecolor','none', 'markerfacecolor','k', 'markerfacealpha', 0.2);
bar(2, mean(cell2mat(taoLat(:,2))), 'facecolor', 'none', 'edgecolor', 'k');
scatter(normrnd(2,0.1, size(cell2mat(taoLat(:,2)))), cell2mat(taoLat(:,2)), 20, 'markeredgecolor','none', 'markerfacecolor','k', 'markerfacealpha', 0.2);
errorbar(1:2, [mean(cell2mat(taoLat(:,1))), mean(cell2mat(taoLat(:,2)))], [std(cell2mat(taoLat(:,1))), std(cell2mat(taoLat(:,2)))], 'linestyle', 'none', 'color', 'k', 'capsize', 0);
title('TAO Response Latency');
ylabel('Poke Duration (s)');
set(gca, 'xtick', 1:2, 'xticklabel', [{'Correct'}, {'Incorrect'}], 'xticklabelrotation', 45);

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 2: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%  Behavior Summary 3
figure;
sp1 = subplot(2,1,1);
bar(mlb.lagVect, mean(lagAcc), 'edgecolor', 'k', 'facecolor', 'none');
hold on;
errorbar(mlb.lagVect, mean(lagAcc), std(lagAcc), 'linestyle', 'none', 'color', 'k', 'capsize', 0);
for lag = 1:length(mlb.lagVect)
    scatter(normrnd(mlb.lagVect(lag), 0.01, size(lagAcc(:,lag))), lagAcc(:,lag), 'markeredgecolor', 'none', 'markerfacecolor', 'k');
end
    
plot(mlb.lagVect, lagAcc', 'color', 'k');
ylabel('Accuracy (% correct trials)');
title('Accuracy By Lag (Mean +/- StD)');
box off;

sp2 = subplot(2,1,2);
hold on;
for lag = 1:length(mlb.lagVect)
    curCorrLag = cell2mat(lagLat(:,lag,1));
    curInCorrLag = cell2mat(lagLat(:,lag,2));
    bar(mlb.lagVect(lag)-0.1, mean(curCorrLag), 0.2, 'facecolor', 'none', 'edgecolor', 'k');
    scatter(normrnd(mlb.lagVect(lag)-0.1, 0.01, size(curCorrLag)), curCorrLag, 20, 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.2);
    bar(mlb.lagVect(lag)+0.1, mean(curInCorrLag), 0.2, 'facecolor', 'none', 'edgecolor', 'k');
    scatter(normrnd(mlb.lagVect(lag)+0.1, 0.01, size(curInCorrLag)), curInCorrLag, 20, 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.2);
end
title('Poke Duration By Lag');
linkaxes([sp1 sp2], 'x');
set([sp1 sp2], 'xlim', [min(mlb.lagVect)-1, max(mlb.lagVect)+1], 'xtick', mlb.lagVect);
xlabel('Lag');
ylabel('Poke Duration (s)');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 3: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%  Behavior Summary 4
figure;
sp1 = subplot(2,1,1);
bar(mlb.lagVect, mean(lagAcc), 'edgecolor', 'k', 'facecolor', 'none');
hold on;
errorbar(mlb.lagVect, mean(lagAcc), std(lagAcc), 'linestyle', 'none', 'color', 'k', 'capsize', 0);
for lag = 1:length(mlb.lagVect)
    scatter(normrnd(mlb.lagVect(lag), 0, size(lagAcc(:,lag))), lagAcc(:,lag), 'markeredgecolor', 'none', 'markerfacecolor', 'k');
end
    
plot(mlb.lagVect, lagAcc', 'color', 'k');
ylabel('Accuracy (% correct trials)');
title('Accuracy By Lag (Mean +/- StD)');
box off;

sp2 = subplot(2,1,2);
hold on;
meanAniLag = cellfun(@(a)mean(a), lagLat);
for lag = 1:length(mlb.lagVect)
    curCorrLag = meanAniLag(:,lag,1);
    curInCorrLag = meanAniLag(:,lag,2);
    bar(mlb.lagVect(lag)-0.1, mean(curCorrLag, 'omitnan'), 0.2, 'facecolor', 'none', 'edgecolor', 'k');
    scatter(normrnd(mlb.lagVect(lag)-0.1, 0, size(curCorrLag)), curCorrLag, 20, 'markeredgecolor', 'none', 'markerfacecolor', 'k');
    bar(mlb.lagVect(lag)+0.1, mean(curInCorrLag, 'omitnan'), 0.2, 'facecolor', 'none', 'edgecolor', 'k');
    scatter(normrnd(mlb.lagVect(lag)+0.1, 0, size(curInCorrLag)), curInCorrLag, 20, 'markeredgecolor', 'none', 'markerfacecolor', 'k');
end
plot(mlb.lagVect-0.1, meanAniLag(:,:,1)', 'color', 'k');
plot(mlb.lagVect+0.1, meanAniLag(:,:,2)', '--k');
title('Poke Duration By Lag');
linkaxes([sp1 sp2], 'x');
set([sp1 sp2], 'xlim', [min(mlb.lagVect)-1, max(mlb.lagVect)+1], 'xtick', mlb.lagVect);
xlabel('Lag');
ylabel('Poke Duration (s)');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 4: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Behavior Summary 5
nonISlagAcc = lagAcc;
nonISlagAcc(:,mlb.lagVect==0) = [];
nonISClagLat = meanAniLag(:,:,1);
nonISClagLat(:,mlb.lagVect==0) = [];
figure;
corrScatPlot(nonISlagAcc(:), nonISClagLat(:), 'Accuracy', 'OSC Latency', []);

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 5: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Behavior Summary 6
figure; 
if strcmp(behavLatType,'rel')
    histBins = -1.5:0.05:1;
else
    histBins = 0:0.05:2;
end
histogram(cell2mat(pokeLatRaw(:,1,1)),histBins); 
hold on; 
histogram(cell2mat(pokeLatRaw(:,2,1)),histBins); 
histogram(cell2mat(pokeLatRaw(:,1,2)),histBins); 
histogram(cell2mat(pokeLatRaw(:,2,2)),histBins);
legend([{'InSeq Correct'}, {'InSeq Incorrect'}, {'OutSeq Correct'}, {'OutSeq Incorrect'}]);
xlabel('Time (s)');
ylabel('Trial Counts');
title('Poke Duration by Trial Type');

annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Behavior Summary 6: SFP = %i", sfpYN),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
%% Plot Neural
grpPV = cell2mat(popVects);
% grpPVsort = timeVect(cell2mat(popVectsSortVect));
grpPVsort = cell2mat(popVectsSortVect);
% grpPVmaxFR(sum(grpPVmaxFR>=frThresh,2)>=1,:) = 100;
% grpPVmaxFR(grpPVmaxFR~=100)=0;
grpPVmaxFR = cell2mat(popVectsThreshVect);

% Plot Ensemble Vectors
pokeOutLats = cell2mat(fiscPokeOutLat(:));
rwdSigLat = cell2mat(fiscRwdSigLat(:));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
figure;
ndxCorr = nan(mlb.seqLength);
pvCorr = nan(mlb.seqLength);
pvCorrTime = nan(mlb.seqLength, mlb.seqLength, length(timeVect));
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        sortedPV = sortrows([grpPVsort(:,op1), grpPVmaxFR(:,op2), grpPV(:,:,op2)]);
        subplot(mlb.seqLength,mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
        imagesc(timeVect, 1:sum(sortedPV(:,2)>=frThresh), sortedPV(sortedPV(:,2)>=frThresh,3:end), [0 1]);
%         imagesc(timeVect, 1:sum(sortedPV(:,2)>=frThresh), sortedPV(sortedPV(:,2)>=frThresh,3:end), [0 20]);
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(repmat(mean(pokeOutLats), [1 2]), get(gca, 'ylim'), '-k')
        plot(repmat(mean(pokeOutLats)-std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
        plot(repmat(mean(pokeOutLats)+std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
        
        plot(repmat(mean(rwdDelivLat), [1 2]), get(gca, 'ylim'), '-w')
        plot(repmat(mean(rwdDelivLat)-std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--w');
        plot(repmat(mean(rwdDelivLat)+std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--w');
        set(gca, 'xtick', [0, mean(pokeOutLats), mean(rwdDelivLat)], 'xticklabel', [{'Poke'}, {'Exit'}, {'Rwd'}], 'xticklabelrotation', 45);
        
        title(sprintf('Pos:%.0f, Sort:%.0f', op1, op2)); 
        
        ndxCorr(op1,op2) = corr(grpPVsort(:,op1), grpPVsort(:,op2), 'rows', 'pairwise');
        tempX = grpPV(grpPVmaxFR(:,op2)>=frThresh,:,op1);
        tempY = grpPV(grpPVmaxFR(:,op2)>=frThresh,:,op2);
        pvCorr(op1,op2) = corr(tempX(:), tempY(:), 'Rows', 'pairwise');
        for t = 1:length(timeVect)
            pvCorrTime(op1,op2,t) = corr(tempX(:,t), tempY(:,t), 'Rows', 'pairwise');
        end
    end
end
linkaxes
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Population Vectors Sortings:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

%% Plot Ensemble Similarity Summaries
pvCorrByPos = nan(mlb.seqLength, mlb.seqLength, size(grpPV,1));
for uni = 1:size(grpPV,1)
    for op1 = 1:mlb.seqLength
        for op2 = 1:mlb.seqLength
            pvCorrByPos(op2,op1,uni) = corr(grpPV(uni,:,op1)', grpPV(uni,:,op2)');
%             pvCorrByPos(op2,op1,uni) = 1-pdist([grpPV(uni,:,op1); grpPV(uni,:,op2)], 'cosine');
        end
    end
end
figure; 
subplot(1,2,1);
imagesc(mean(pvCorrByPos,3, 'omitnan'), [0 0.5]);
title('Individual Cell FR Similarity');
colorbar;
subplot(1,2,2);
imagesc(pvCorr, [0 0.75]);
title('Population Vector Similarity');
colorbar;
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Trial Position Correlations:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Plot PV Correlation Over Time
figure;
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
        plot(timeVect, reshape(pvCorrTime(op1,op2,:), [1,size(pvCorrTime,3)]));
        axis tight;
        set(gca, 'ylim', [0 1]);
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(repmat(mean(pokeOutLats), [1 2]), get(gca, 'ylim'), '-k');
        plot(repmat(mean(pokeOutLats)-std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
        plot(repmat(mean(pokeOutLats)+std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
        plot(repmat(mean(rwdDelivLat), [1 2]), get(gca, 'ylim'), '-r');
        plot(repmat(mean(rwdDelivLat)-std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--r');
        plot(repmat(mean(rwdDelivLat)+std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--r');
        set(gca, 'xtick', [0, mean(pokeOutLats), mean(rwdDelivLat)], 'xticklabel', [{'Poke'}, {'Exit'}, {'Rwd'}], 'xticklabelrotation', 45);
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("PV Correlation Across Time:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Plot Temporal Remapping
figure;
propStable = cell(mlb.seqLength);
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
%         corrScatPlot(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op1), grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op2), op1,op2,[]);
        [tempCounts, xedge,yedge] = histcounts2(timeVect(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op1)), timeVect(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op2)),downsample(timeVect,binSize/dsRate),downsample(timeVect,binSize/dsRate));
        for t = 1:size(tempCounts,1)
            tempCounts(:,t) = tempCounts(:,t)/sum(tempCounts(:,t));
        end
        imagesc(xedge(2:end)-mode(diff(xedge))/2, yedge(2:end)-mode(diff(yedge))/2, tempCounts, [0 0.5]);
        propStable{op2,op1} = tempCounts(logical(eye(size(tempCounts))));
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');                                                plot(get(gca, 'xlim'), [0 0], '-k');                
        plot(repmat(mean(pokeOutLats), [1 2]), get(gca, 'ylim'), '-k');                     plot(get(gca, 'xlim'), repmat(mean(pokeOutLats), [1 2]), '-k')
        plot(repmat(mean(pokeOutLats)-std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');   plot(get(gca, 'xlim'),repmat(mean(pokeOutLats)-std(pokeOutLats), [1,2]), '--k');
        plot(repmat(mean(pokeOutLats)+std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');   plot(get(gca, 'xlim'),repmat(mean(pokeOutLats)+std(pokeOutLats), [1,2]), '--k');
        
        plot(repmat(mean(rwdDelivLat), [1 2]), get(gca, 'ylim'), '-w');                     plot(get(gca, 'xlim'),repmat(mean(rwdDelivLat), [1 2]), '-w');
        plot(repmat(mean(rwdDelivLat)-std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--w');   plot(get(gca, 'xlim'),repmat(mean(rwdDelivLat)-std(rwdDelivLat), [1,2]), '--w');
        plot(repmat(mean(rwdDelivLat)+std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--w');   plot(get(gca, 'xlim'),repmat(mean(rwdDelivLat)+std(rwdDelivLat), [1,2]), '--w');
        set(gca, 'xtick', [0, mean(pokeOutLats), mean(rwdDelivLat)], 'xticklabel', [{'Poke'}, {'Exit'}, {'Rwd'}], 'xticklabelrotation', 20,...
            'ytick', [0, mean(pokeOutLats), mean(rwdDelivLat)], 'yticklabel', [{'Poke'}, {'Exit'}, {'Rwd'}]);
        
        set(gca, 'ydir','normal');
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Temporal Remapping:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Plot Proportion Stable Cells
permPropStable = cellfun(@(a){a'},propStable);
proStableDiag = permPropStable(~logical(eye(4)));
figure; 
plot(xedge(2:end)-mode(diff(xedge))/2, cell2mat(proStableDiag)');
hold on;
plot([0 0], get(gca, 'ylim'), '-k');
plot(repmat(mean(pokeOutLats), [1 2]), get(gca, 'ylim'), '-k');
plot(repmat(mean(pokeOutLats)-std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
plot(repmat(mean(pokeOutLats)+std(pokeOutLats), [1,2]), get(gca, 'ylim'), '--k');
plot(repmat(mean(rwdDelivLat), [1 2]), get(gca, 'ylim'), '-r');
plot(repmat(mean(rwdDelivLat)-std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--r');
plot(repmat(mean(rwdDelivLat)+std(rwdDelivLat), [1,2]), get(gca, 'ylim'), '--r');


% Plot Max Rate & Time Remapping
rateTimeRemap = nan(size(grpPVmaxFR,1),2);
timeRateRemap = nan(size(grpPVmaxFR,1),2);
meanMaxFR = mean(grpPVmaxFR,2);
for u = 1:size(grpPVmaxFR,1)
    diffMtx = squareform(pdist(grpPVmaxFR(u,:)'));
    rateTimeRemap(u,1) = max(max(diffMtx));
    [r,c] = find(diffMtx==max(max(diffMtx)) & triu(true(mlb.seqLength)));
    rateTimeRemap(u,2) = abs(grpPVsort(u,r(1))-grpPVsort(u,c(1)));
    
    diffMtx = squareform(pdist(grpPVsort(u,:)'));
    timeRateRemap(u,1) = max(max(diffMtx));
    [r,c] = find(diffMtx==max(max(diffMtx)) & triu(true(mlb.seqLength)));
    timeRateRemap(u,2) = abs(grpPVmaxFR(u,r(1))-grpPVmaxFR(u,c(1)));
end
figure; 
subplot(2,2,1)
corrScatPlot(rateTimeRemap(:,1), rateTimeRemap(:,2), 'Max \DeltaFR', 'Max \DeltaLatency', []);
title('Latency Remapping at Max Rate Remap')
subplot(2,2,2)
corrScatPlot(timeRateRemap(:,1), timeRateRemap(:,2), 'Max \DeltaLatency', 'Max \DeltaFR', []);
title('Rate Remapping at Max Latency Remap');
subplot(2,2,3)
corrScatPlot(meanMaxFR, rateTimeRemap(:,1), 'Mean Max FR', 'Max \DeltaFR', []);
title('Mean FR vs Rate Remap');
subplot(2,2,4)
corrScatPlot(meanMaxFR, timeRateRemap(:,1), 'Mean Max FR', 'Max \DeltaLatency',[]);
title('Mean FR vs Latency Remap');
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Rate vs Time Remapping:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

% Plot Rate & Time Remapping Across Positions
figure;
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        if op1 ~= op2
            subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
            corrScatPlot(abs(grpPVmaxFR(grpPVmaxFR(:,op2)>=frThresh,op1)-grpPVmaxFR(grpPVmaxFR(:,op2)>=frThresh,op2)),abs(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op1)-grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op2)), '\Delta Rate', '\Delta Time',[]);
            set(gca, 'ydir','normal');
            title(sprintf('Pos %i vs %i', op1,op2));
        end
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Rate vs Temporal Remapping:\n     binType = %s, binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binType, binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');



    