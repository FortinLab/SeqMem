% PFC_Behavior

%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];
ssnTtl = 'PFC Well-Trained';

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
% ssnTtl = 'CA1 Well-Trained';

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
ssnTtl = 'PFC Well-Trained';

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
% ssnTtl = 'PFC Dual-List';

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
respMtx = cell2mat(permute(cellfun(@(a){a.responseMatrix}, mlb), [1,3,2]));
respMtxSFP = cell2mat(permute(cellfun(@(a){a.responseMatrixSFP}, mlb), [1,3,2]));

%% Transposition Matrix
tempTransMatAcc = permute(cellfun(@(a){a.transMatPerf}, mlb), [1,3,2]);
transMatC = sum(cell2mat(cellfun(@(a){a(:,:,1)}, tempTransMatAcc)),3);
transMatIC = sum(cell2mat(cellfun(@(a){a(:,:,2)}, tempTransMatAcc)),3);
transMatAcc = transMatC./(transMatC+transMatIC);

figure;
imagesc(transMatAcc, [0 1]);
for r = 1:mlb{ani}.seqLength
    for c = 1:mlb{ani}.seqLength
        if isnan(transMatAcc(r,c))
            patch([c-0.5 c+0.5 c+0.5 c-0.5], [r-0.5 r-0.5 r+0.5 r+0.5], 'w', 'linestyle', 'none');
        end
    end
end
title('Trial Accuracy (Across Animals)');
xlabel('Positions');
ylabel('Odor');
set(gca,'xtick', 1:mlb{ani}.seqLength, 'ytick', 1:mlb{ani}.seqLength, 'yticklabels', mlb{ani}.Rosetta(1:mlb{ani}.seqLength));
colormap(cMap);
%% Evaluate Accuracy & Latency by Lag
binNum = 20;
latThresh = 2; % There's an error somewhere in my latency calculation code where it miscalculates the hold duration for trials. By my count it only affects three trials in the data set, only one aggregiously so. Sorry, I'm too lazy to track down a bug for such a small subset of trials.

tempLagLatAll = permute(cellfun(@(a){permute(a.lagLatVectRaw, [3,2,1])}, mlb), [1,3,2]);
tempLagLatAll = reshape([tempLagLatAll{:}], [size(tempLagLatAll{1},1), size(tempLagLatAll{1},2), length(fileDirs)]);
grpLagLat = cell(size(tempLagLatAll(:,:,1)));
figure;
ndx = 1:3:size(grpLagLat,2)*3;
subplot(2,2,1);
aniLagAcc = nan(length(fileDirs),size(grpLagLat,2));
for ani = 1:length(fileDirs)
    tempAccVect = cell2mat(cellfun(@(a){length(a)}, tempLagLatAll(:,:,ani)));
    aniLagAcc(ani,:) = tempAccVect(1,:)./sum(tempAccVect,1);
end
for n = 2:length(ndx)
    mlb{ani}.PlotMeanVarSwarmBar(ndx(n), aniLagAcc(:,n),1,0.05, 'k', 'error', 'SEM', 'binNum', 5);
end
set(gca,'ylim', [0 1.2], 'xlim', [min(ndx)-1, max(ndx)+1], 'xtick', ndx, 'xticklabels', (1:size(grpLagLat,2))-mlb{ani}.seqLength);
title('Accuracy by Lag (Mean +/- SEM across Animals)');
xlabel('Lag');
ylabel('Accuracy (% correct)');

subplot(2,2,3);
for r = 1:size(grpLagLat,1)
    for c = 1:size(grpLagLat,2)
        grpLagLat{r,c} = cell2mat(squeeze(tempLagLatAll(r,c,:)));
        grpLagLat{r,c}(grpLagLat{r,c}>latThresh) = [];
        if ~isempty(grpLagLat{r,c})
            if r ==1
                if c==mlb{ani}.seqLength
                    corr = mlb{ani}.PlotMeanVarViolin(ndx(c)-0.5, grpLagLat{r,c}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
                else
                    mlb{ani}.PlotMeanVarViolin(ndx(c)-0.5, grpLagLat{r,c}, 1, 0.05, 'k', 'filled','error', 'SEM', 'binNum', binNum);
                end                    
            else
                if c==mlb{ani}.seqLength
                    inCorr = mlb{ani}.PlotMeanVarViolin(ndx(c)+0.5, grpLagLat{r,c}, 1, 0.05, 'r', 'filled','error', 'SEM', 'binNum', binNum);
                else
                    mlb{ani}.PlotMeanVarViolin(ndx(c)+0.5, grpLagLat{r,c}, 1, 0.05, 'r', 'filled','error', 'SEM', 'binNum', binNum);
                end
            end
        end
    end
end
legend([corr, inCorr], {'Correct', 'Incorrect'});
set(gca,'ylim', [0 2], 'xlim', [min(ndx)-1, max(ndx)+1], 'xtick', ndx, 'xticklabels', (1:size(grpLagLat,2))-mlb{ani}.seqLength);
title('Response Latency by Lag (All Trials)');
xlabel('Lag');
ylabel('Response Latency (s)');

tempLagLatSFP = permute(cellfun(@(a){permute(a.lagLatVectSFPraw, [3,2,1])}, mlb), [1,3,2]);
tempLagLatSFP = reshape([tempLagLatSFP{:}], [size(tempLagLatSFP{1},1), size(tempLagLatSFP{1},2), length(fileDirs)]);
grpLagLat = cell(size(tempLagLatSFP(:,:,1)));
ndx = 1:3:size(grpLagLat,2)*3;
subplot(2,2,2);
aniLagAcc = nan(length(fileDirs),size(grpLagLat,2));
for ani = 1:length(fileDirs)
    tempAccVect = cell2mat(cellfun(@(a){length(a)}, tempLagLatSFP(:,:,ani)));
    aniLagAcc(ani,:) = tempAccVect(1,:)./sum(tempAccVect,1);
end
for n = 2:length(ndx)
    mlb{ani}.PlotMeanVarSwarmBar(ndx(n), aniLagAcc(:,n),1,0.05, 'k', 'error', 'SEM', 'binNum', 5);
end
set(gca,'ylim', [0 1.2], 'xlim', [min(ndx)-1, max(ndx)+1], 'xtick', ndx, 'xticklabels', (1:size(grpLagLat,2))-mlb{ani}.seqLength);
title('Accuracy by Lag SFP (Mean +/- SEM across Animals)');
xlabel('Lag');
ylabel('Accuracy (% correct)');

subplot(2,2,4);
for r = 1:size(grpLagLat,1)
    for c = 1:size(grpLagLat,2)
        grpLagLat{r,c} = cell2mat(squeeze(tempLagLatSFP(r,c,:)));
        grpLagLat{r,c}(grpLagLat{r,c}>latThresh) = [];
        if ~isempty(grpLagLat{r,c})
            if r ==1
                if c==mlb{ani}.seqLength
                    corr = mlb{ani}.PlotMeanVarViolin(ndx(c)-0.5, grpLagLat{r,c}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
                else
                    mlb{ani}.PlotMeanVarViolin(ndx(c)-0.5, grpLagLat{r,c}, 1, 0.05, 'k', 'filled','error', 'SEM', 'binNum', binNum);
                end                    
            else
                if c==mlb{ani}.seqLength
                    inCorr = mlb{ani}.PlotMeanVarViolin(ndx(c)+0.5, grpLagLat{r,c}, 1, 0.05, 'r', 'filled','error', 'SEM', 'binNum', binNum);
                else
                    mlb{ani}.PlotMeanVarViolin(ndx(c)+0.5, grpLagLat{r,c}, 1, 0.05, 'r', 'filled','error', 'SEM', 'binNum', binNum);
                end
            end
        end
    end
end
legend([corr, inCorr], {'Correct', 'Incorrect'});
set(gca,'ylim', [0 2], 'xlim', [min(ndx)-1, max(ndx)+1], 'xtick', ndx, 'xticklabels', (1:size(grpLagLat,2))-mlb{ani}.seqLength);
title('Response Latency by Lag (SFP)');
xlabel('Lag');
ylabel('Response Latency (s)');

%% Histogram of Response Latencies
bins = 0:0.1:2;
respLats = cell(2);
for ani = 1:length(fileDirs)
    respLats{1,1} = [respLats{1,1}; tempLagLatAll{1,mlb{ani}.seqLength,ani}];
    respLats{1,1}(respLats{1,1}>latThresh) = [];
    respLats{1,2} = [respLats{1,latThresh}; tempLagLatAll{2,mlb{ani}.seqLength,ani}];
    respLats{1,2}(respLats{1,2}>latThresh) = [];
    respLats{2,1} = [respLats{2,1}; cell2mat(reshape(tempLagLatAll(2,(1:size(tempLagLatAll,2))~=mlb{ani}.seqLength,ani), [size(tempLagLatAll,2)-1,1]))];
    respLats{2,1}(respLats{2,1}>latThresh) = [];
    respLats{2,2} = [respLats{2,2}; cell2mat(reshape(tempLagLatAll(1,(1:size(tempLagLatAll,2))~=mlb{ani}.seqLength,ani), [size(tempLagLatAll,2)-1,1]))];
    respLats{2,2}(respLats{2,2}>latThresh) = [];
end
figure;
subplot(1,3,1:2);
for n = 1:numel(respLats)
    histogram(respLats{n}, bins, 'facealpha', 0.25);
    hold on;
end
subplot(1,3,3);
mlb{ani}.PlotMeanVarViolin(1,respLats{1,1}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
mlb{ani}.PlotMeanVarViolin(2,respLats{1,2}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
mlb{ani}.PlotMeanVarViolin(3,respLats{2,1}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
mlb{ani}.PlotMeanVarViolin(4,respLats{2,2}, 1, 0.05, 'k', 'filled', 'error', 'SEM', 'binNum', binNum);
set(gca, 'xlim', [0 5], 'xtick', 1:4, 'xticklabel', [{'ISC'}, {'ISI'}, {'OSI'}, {'OSC'}]);

anovaLog = cell2mat(reshape(cellfun(@(a,b){ones(size(a))+b}, respLats, reshape(num2cell(0:3), [2,2])), [4,1]));
[p,tbl,stats] = anovan(cell2mat(respLats(:)), anovaLog, 'display', 'off');
figure;
mc = multcompare(tbl, 'cType', 'bonferroni');
%% Response Bias Measures
