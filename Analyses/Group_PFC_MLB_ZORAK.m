fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
aniIDs = [{'GE11'}, {'GE13'}, {'GE14'}, {'GE17'}, {'GE24'}];
binSize = 200;
dsRate = 50;

%%
tic;
for fl = 1:length(fileDirs)
    pfcMLB(fl) = PFC_TrialEvent_MLB_SM(fileDirs{fl}, binSize, dsRate); %#ok<SAGROW>
%     pfcMLB(fl).binSize = binSize;
%     pfcMLB(fl).dsRate = dsRate;
%     pfcMLB(fl).beginTrialWindow = [-500 1200];
%     pfcMLB(fl).endTrialWindow = [0 0];
%     pfcMLB(fl).RunAnalysis;
end
toc
%%
%% Comment this out if running straight or just run line by line below here if loading the data
% save('PFCmlbData.mat', 'pfcMLB', '-v7.3');
%% Trial After OutSeq Decoding
figure; 
subplot(3,1,1)
taoDecodes = {pfcMLB.taoDecode};
taoPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), taoDecodes, 'uniformoutput', 0));
taoPrevPosVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), taoDecodes, 'uniformoutput',0));
taoPrevOdrVals = cell2mat(cellfun(@(a)nanmean(a==-2,2), taoDecodes, 'uniformoutput', 0));
taoOtherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), taoDecodes, 'uniformoutput', 0));

meanPos = mean(taoPosVals,2);
stdPos = std(taoPosVals,1,2)./sqrt(size(taoPosVals,2)-1);
% stdPos = std(taoPosVals,1,2);
a = plot(meanPos, 'r', 'linewidth', 2);
hold on;
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanPrevPos = mean(taoPrevPosVals,2);
stdPrevPos = std(taoPrevPosVals,1,2)./sqrt(size(taoPrevPosVals,2)-1);
% stdPrevPos = std(taoPrevPosVals,1,2);
b = plot(meanPrevPos, 'b', 'linewidth', 2);
patch([1:length(meanPrevPos), length(meanPrevPos):-1:1], [meanPrevPos+stdPrevPos; flipud(meanPrevPos-stdPrevPos)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanPrevOdr = mean(taoPrevOdrVals,2);
stdPrevOdr = std(taoPrevOdrVals,1,2)./sqrt(size(taoPrevOdrVals,2)-1);
% stdPrevOdr = std(taoPrevOdrVals,1,2);
c = plot(meanPrevOdr, 'g', 'linewidth', 2);
patch([1:length(meanPrevOdr), length(meanPrevOdr):-1:1], [meanPrevOdr+stdPrevOdr; flipud(meanPrevOdr-stdPrevOdr)],...
    'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g');
meanOther = mean(taoOtherVals,2);
stdOther = std(taoOtherVals,1,2)./sqrt(size(taoOtherVals,2)-1);
d = plot(meanOther, 'k', 'linewidth', 2);
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Trial After OutSeq');
legend([a b c d], [{'CurPos'}, {'PrevPos'}, {'PrevOdr'}, {'Other'}]);

subplot(3,1,2)
taoDecodes = {pfcMLB.taSkpDecode};
taoPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), taoDecodes, 'uniformoutput', 0));
taoPrevPosVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), taoDecodes, 'uniformoutput',0));
taoPrevOdrVals = cell2mat(cellfun(@(a)nanmean(a==-2,2), taoDecodes, 'uniformoutput', 0));
taoOtherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), taoDecodes, 'uniformoutput', 0));

meanPos = mean(taoPosVals,2);
stdPos = std(taoPosVals,1,2)./sqrt(size(taoPosVals,2)-1);
% stdPos = std(taoPosVals,1,2);
a = plot(meanPos, 'r', 'linewidth', 2);
hold on;
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanPrevPos = mean(taoPrevPosVals,2);
stdPrevPos = std(taoPrevPosVals,1,2)./sqrt(size(taoPrevPosVals,2)-1);
% stdPrevPos = std(taoPrevPosVals,1,2);
b = plot(meanPrevPos, 'b', 'linewidth', 2);
patch([1:length(meanPrevPos), length(meanPrevPos):-1:1], [meanPrevPos+stdPrevPos; flipud(meanPrevPos-stdPrevPos)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanPrevOdr = mean(taoPrevOdrVals,2);
stdPrevOdr = std(taoPrevOdrVals,1,2)./sqrt(size(taoPrevOdrVals,2)-1);
% stdPrevOdr = std(taoPrevOdrVals,1,2);
c = plot(meanPrevOdr, 'g', 'linewidth', 2);
patch([1:length(meanPrevOdr), length(meanPrevOdr):-1:1], [meanPrevOdr+stdPrevOdr; flipud(meanPrevOdr-stdPrevOdr)],...
    'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g');
meanOther = mean(taoOtherVals,2);
stdOther = std(taoOtherVals,1,2)./sqrt(size(taoOtherVals,2)-1);
d = plot(meanOther, 'k', 'linewidth', 2);
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Trial After Skips');
% legend([a b c d], [{'CurPos'}, {'PrevPos'}, {'PrevOdr'}, {'Other'}]);

subplot(3,1,3)
taoDecodes = {pfcMLB.taRepDecode};
taoPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), taoDecodes, 'uniformoutput', 0));
taoPrevPosVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), taoDecodes, 'uniformoutput',0));
taoPrevOdrVals = cell2mat(cellfun(@(a)nanmean(a==-2,2), taoDecodes, 'uniformoutput', 0));
taoOtherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), taoDecodes, 'uniformoutput', 0));

meanPos = mean(taoPosVals,2);
stdPos = std(taoPosVals,1,2)./sqrt(size(taoPosVals,2)-1);
% stdPos = std(taoPosVals,1,2);
a = plot(meanPos, 'r', 'linewidth', 2);
hold on;
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanPrevPos = mean(taoPrevPosVals,2);
stdPrevPos = std(taoPrevPosVals,1,2)./sqrt(size(taoPrevPosVals,2)-1);
% stdPrevPos = std(taoPrevPosVals,1,2);
b = plot(meanPrevPos, 'b', 'linewidth', 2);
patch([1:length(meanPrevPos), length(meanPrevPos):-1:1], [meanPrevPos+stdPrevPos; flipud(meanPrevPos-stdPrevPos)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanPrevOdr = mean(taoPrevOdrVals,2);
stdPrevOdr = std(taoPrevOdrVals,1,2)./sqrt(size(taoPrevOdrVals,2)-1);
% stdPrevOdr = std(taoPrevOdrVals,1,2);
c = plot(meanPrevOdr, 'g', 'linewidth', 2);
patch([1:length(meanPrevOdr), length(meanPrevOdr):-1:1], [meanPrevOdr+stdPrevOdr; flipud(meanPrevOdr-stdPrevOdr)],...
    'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g');
meanOther = mean(taoOtherVals,2);
stdOther = std(taoOtherVals,1,2)./sqrt(size(taoOtherVals,2)-1);
d = plot(meanOther, 'k', 'linewidth', 2);
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Trial After Repeats');
% legend([a b c d], [{'CurPos'}, {'PrevPos'}, {'PrevOdr'}, {'Other'}]);
%% TAO Pre-Trial Bar Plots
taoDecodes = {pfcMLB.taoDecode};
taoPosValsAni = nan(1,length(fileDirs));
taoPrevPosValsAni = nan(1,length(fileDirs));
taoPrevOdrValsAni = nan(1,length(fileDirs));
for a = 1:length(fileDirs)
    taoPosValsAni(a) = mean(mean(taoDecodes{a}(pfcMLB(a).trialPeriodTimeLog==1,:)==0));
    taoPrevPosValsAni(a) = mean(mean(taoDecodes{a}(pfcMLB(a).trialPeriodTimeLog==1,:)==-1));
    taoPrevOdrValsAni(a) = mean(mean(taoDecodes{a}(pfcMLB(a).trialPeriodTimeLog==1,:)==-2));
end
figure;
BarPlotErrorbars([mean(taoPosValsAni), mean(taoPrevPosValsAni), mean(taoPrevOdrValsAni)],...
    [SEMcalc(taoPosValsAni'), SEMcalc(taoPrevPosValsAni'), SEMcalc(taoPrevOdrValsAni')]);
set(gca, 'xticklabel', [{'Current Position'}, {'Previous Position'}, {'Previous Odor'}], 'xticklabelrotation',45);
ylabel('Decoding Probability');
title('Trial After OutSeq Pre-Trial Decoding');
[p,tbl,stats] = anova1([taoPosValsAni', taoPrevPosValsAni', taoPrevOdrValsAni']);
multcompare(stats, 'CType', 'bonferroni');

%% OutSeq Decoding
figure;
subplot(3,1,1)
osDecodes = {pfcMLB.osDecode};
osPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), osDecodes, 'uniformoutput', 0));
osOdrVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), osDecodes, 'uniformoutput', 0));
otherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), osDecodes, 'uniformoutput', 0));
meanPos = nanmean(osPosVals,2);
stdPos = std(osPosVals,1,2)./sqrt(size(osPosVals,2)-1);
% stdPos = std(osPosVals,1,2);
a = plot(meanPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanOdr = nanmean(osOdrVals,2);
stdOdr = std(osOdrVals,1,2)./sqrt(size(osOdrVals,2)-1);
% stdOdr = std(osOdrVals,1,2);
b = plot(meanOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanOdr), length(meanOdr):-1:1], [meanOdr+stdOdr; flipud(meanOdr-stdOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanOther = nanmean(otherVals,2);
stdOther = std(otherVals,1,2)./sqrt(size(otherVals,2)-1);
% stdOther = std(otherVals,1,2);
c = plot(meanOther, 'k', 'linewidth', 2);
hold on
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
title('All OutSeq');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
legend([a b c], [{'Position'}, {'Odor'}, {'Other'}]);
set(gca, 'ylim', [0 1]);

subplot(3,1,2)
repDecodes = {pfcMLB.repDecode};
repPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), repDecodes, 'uniformoutput', 0));
repOdrVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), repDecodes, 'uniformoutput', 0));
repOtherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), repDecodes, 'uniformoutput', 0));
meanRepPos = mean(repPosVals,2);
stdRepPos = std(repPosVals,1,2)./sqrt(size(repPosVals,2)-1);
% stdPos = std(osPosVals,1,2);
plot(meanRepPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanRepPos), length(meanRepPos):-1:1], [meanRepPos+stdRepPos; flipud(meanRepPos-stdRepPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanRepOdr = mean(repOdrVals,2);
stdRepOdr = std(repOdrVals,1,2)./sqrt(size(repOdrVals,2)-1);
% stdOdr = std(osOdrVals,1,2);
plot(meanRepOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanRepOdr), length(meanRepOdr):-1:1], [meanRepOdr+stdRepOdr; flipud(meanRepOdr-stdRepOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanRepOther = mean(repOtherVals,2);
stdRepOther = std(repOtherVals,1,2)./sqrt(size(repOtherVals,2)-1);
% stdOther = std(otherVals,1,2);
plot(meanRepOther, 'k', 'linewidth', 2);
hold on
patch([1:length(meanRepOther), length(meanRepOther):-1:1], [meanRepOther+stdRepOther; flipud(meanRepOther-stdRepOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Repeats Only');
set(gca, 'ylim', [0 1]);


subplot(3,1,3)
skpDecodes = {pfcMLB.skpDecode};
skpPosVals = cell2mat(cellfun(@(a)nanmean(a==0,2), skpDecodes, 'uniformoutput', 0));
skpOdrVals = cell2mat(cellfun(@(a)nanmean(a==-1,2), skpDecodes, 'uniformoutput', 0));
skpOtherVals = cell2mat(cellfun(@(a)nanmean(a==1,2), skpDecodes, 'uniformoutput', 0));
meanSkpPos = mean(skpPosVals,2);
stdSkpPos = std(skpPosVals,1,2)./sqrt(size(skpPosVals,2)-1);
% stdPos = std(osPosVals,1,2);
plot(meanSkpPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanSkpPos), length(meanSkpPos):-1:1], [meanSkpPos+stdSkpPos; flipud(meanSkpPos-stdSkpPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanSkpOdr = mean(skpOdrVals,2);
stdMeanOdr = std(skpOdrVals,1,2)./sqrt(size(skpOdrVals,2)-1);
% stdOdr = std(osOdrVals,1,2);
plot(meanSkpOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanSkpOdr), length(meanSkpOdr):-1:1], [meanSkpOdr+stdMeanOdr; flipud(meanSkpOdr-stdMeanOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
meanSkpOther = mean(skpOtherVals,2);
stdSkpOther = std(skpOtherVals,1,2)./sqrt(size(skpOtherVals,2)-1);
% stdOther = std(otherVals,1,2);
plot(meanSkpOther, 'k', 'linewidth', 2);
hold on
patch([1:length(meanSkpOther), length(meanSkpOther):-1:1], [meanSkpOther+stdSkpOther; flipud(meanSkpOther-stdSkpOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Skips Only');
set(gca, 'ylim', [0 1]);

%% OutSeq Decoding
figure;
subplot(3,1,1)
a = plot(meanPos, 'k', 'linewidth', 2);
hold on
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
b = plot(meanRepPos, 'b', 'linewidth', 2);
patch([1:length(meanRepPos), length(meanRepPos):-1:1], [meanRepPos+stdRepPos; flipud(meanRepPos-stdRepPos)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
c = plot(meanSkpPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanSkpPos), length(meanSkpPos):-1:1], [meanSkpPos+stdSkpPos; flipud(meanSkpPos-stdSkpPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Position Decoding');
legend([a b c], [{'All OutSeq'}, {'Repeats'}, {'Skips'}]);
set(gca, 'ylim', [0 1]);

subplot(3,1,2)
plot(meanOdr, 'k', 'linewidth', 2);
hold on
patch([1:length(meanOdr), length(meanOdr):-1:1], [meanOdr+stdOdr; flipud(meanOdr-stdOdr)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
plot(meanRepOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanRepOdr), length(meanRepOdr):-1:1], [meanRepOdr+stdRepOdr; flipud(meanRepOdr-stdRepOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
plot(meanSkpOdr, 'r', 'linewidth', 2);
hold on
patch([1:length(meanSkpOdr), length(meanSkpOdr):-1:1], [meanSkpOdr+stdMeanOdr; flipud(meanSkpOdr-stdMeanOdr)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Odor Decoding');
set(gca, 'ylim', [0 1]);

subplot(3,1,3)
plot(meanOther, 'k', 'linewidth', 2);
hold on
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
plot(meanRepOther, 'b', 'linewidth', 2);
hold on
patch([1:length(meanRepOther), length(meanRepOther):-1:1], [meanRepOther+stdRepOther; flipud(meanRepOther-stdRepOther)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
plot(meanSkpOther, 'r', 'linewidth', 2);
hold on
patch([1:length(meanSkpOther), length(meanSkpOther):-1:1], [meanSkpOther+stdSkpOther; flipud(meanSkpOther-stdSkpOther)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Other Decoding');
set(gca, 'ylim', [0 1]);

%% TAO & OutSeq in common figure
figure
sp(1) = subplot(5,4,[1:3, 5:7]);
meanPos = mean(taoPosVals,2);
stdPos = std(taoPosVals,1,2)./sqrt(size(taoPosVals,2)-1);
% stdPos = std(taoPosVals,1,2);
a = plot(meanPos, 'color', [247/255, 148/255, 29/255], 'linewidth', 2);
hold on;
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    [247/255, 148/255, 29/255], 'FaceAlpha', 0.5, 'EdgeColor', [247/255, 148/255, 29/255]);
meanPrevPos = mean(taoPrevPosVals,2);
stdPrevPos = std(taoPrevPosVals,1,2)./sqrt(size(taoPrevPosVals,2)-1);
% stdPrevPos = std(taoPrevPosVals,1,2);
b = plot(meanPrevPos, 'r', 'linewidth', 2);
patch([1:length(meanPrevPos), length(meanPrevPos):-1:1], [meanPrevPos+stdPrevPos; flipud(meanPrevPos-stdPrevPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanPrevOdr = mean(taoPrevOdrVals,2);
stdPrevOdr = std(taoPrevOdrVals,1,2)./sqrt(size(taoPrevOdrVals,2)-1);
% stdPrevOdr = std(taoPrevOdrVals,1,2);
c = plot(meanPrevOdr, 'b', 'linewidth', 2);
patch([1:length(meanPrevOdr), length(meanPrevOdr):-1:1], [meanPrevOdr+stdPrevOdr; flipud(meanPrevOdr-stdPrevOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
axis tight
set(gca, 'ylim', [0 1], 'ytick', 0:0.2:1, 'yticklabel', 0:20:100, 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Trial After OutSeq');
legend([a b c], [{'CurPos'}, {'PrevPos'}, {'PrevOdr'}]);

sp(2) = subplot(5,4,9:11);
plot(zeros(size(meanPos)), 'k');
patch([[size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], [size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1))],...
    [-0.5 0.5 0.5 -0.5], 'w', 'EdgeColor', 'k');
axis off;
axis tight

sp(3) = subplot(5,4,[13:15, 17:19]);
meanPos = nanmean(osPosVals,2);
stdPos = std(osPosVals,1,2)./sqrt(size(osPosVals,2)-1);
% stdPos = std(osPosVals,1,2);
a = plot(meanPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
meanOdr = nanmean(osOdrVals,2);
stdOdr = std(osOdrVals,1,2)./sqrt(size(osOdrVals,2)-1);
% stdOdr = std(osOdrVals,1,2);
b = plot(meanOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanOdr), length(meanOdr):-1:1], [meanOdr+stdOdr; flipud(meanOdr-stdOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
title('All OutSeq');
axis tight
set(gca, 'ylim', [0 1], 'ytick', 0:0.2:1, 'yticklabel', 0:20:100, 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
legend([a b], [{'Position'}, {'Odor'}]);
set(gca, 'ylim', [0 1]);

linkaxes(sp, 'x')

sp2(1) = subplot(5,4,[4,8]);
BarPlotErrorbars([mean(taoPosValsAni), mean(taoPrevPosValsAni), mean(taoPrevOdrValsAni)],...
    [SEMcalc(taoPosValsAni'), SEMcalc(taoPrevPosValsAni'), SEMcalc(taoPrevOdrValsAni')]);
set(gca, 'xticklabel', [{'Current Position'}, {'Previous Position'}, {'Previous Odor'}], 'xticklabelrotation',45, 'tickdir', 'out', 'ylim', [0 0.5], 'ytick', 0:0.25:0.5, 'yticklabel', 0:25:50);
box off
ylabel('Decoding Probability');
title('Trial After OutSeq Pre-Trial Decoding');

sp2(2) = subplot(5,4,[16 20]);
BarPlotErrorbars([mean(osPosValsAni), mean(osOdrValsAni)],...
    [SEMcalc(osPosValsAni'), SEMcalc(osOdrValsAni')]);
set(gca, 'xticklabel', [{'Trial Position'}, {'Trial Odor'}], 'xticklabelrotation',45, 'tickdir', 'out', 'ylim', [0 0.5], 'ytick', 0:0.25:0.5, 'yticklabel', 0:25:50);
box off
ylabel('Decoding Probability');
title('OutSeq Decoding During Odor Sampling');

linkaxes(sp2, 'y');

%% OutSeq During Odor Sampling Bar Plots
osDecodes = {pfcMLB.osDecode};
osPosValsAni = nan(1,length(fileDirs));
osOdrValsAni = nan(1,length(fileDirs));
for a = 1:length(fileDirs)
    osPosValsAni(a) = mean(mean(osDecodes{a}(pfcMLB(a).trialPeriodTimeLog==2 | pfcMLB(a).trialPeriodTimeLog==3,:)==0));
    osOdrValsAni(a) = mean(mean(osDecodes{a}(pfcMLB(a).trialPeriodTimeLog==2 | pfcMLB(a).trialPeriodTimeLog==3,:)==-1));
end
figure;
BarPlotErrorbars([mean(osPosValsAni), mean(osOdrValsAni)],...
    [SEMcalc(osPosValsAni'), SEMcalc(osOdrValsAni')]);
set(gca, 'xticklabel', [{'Trial Position'}, {'Trial Odor'}], 'xticklabelrotation',45);
ylabel('Decoding Probability');
title('OutSeq Decoding During Odor Sampling');
[h,p,ci,stats] = ttest2(osPosValsAni, osOdrValsAni);

%%
odrTimeDecodes = cell2mat(reshape({pfcMLB.fisL1OdecodeOdrTime}, [1 1 length(fileDirs)]));
figure; 
subplot(2,4,1:2);
hold on;
for o = 1:4
    meanDecode = mean(odrTimeDecodes(:,o,:),3);
    stdDecode = std(odrTimeDecodes(:,o,:),1,3)./sqrt(size(odrTimeDecodes,3)-1);
%     stdDecode = std(odrTimeDecodes(:,o,:),1,3);
    plot(meanDecode, 'color', pfcMLB(1).PositionColors(o,:), 'linewidth', 2);
    patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
        pfcMLB(1).PositionColors(o,:), 'FaceAlpha', 0.5, 'EdgeColor', pfcMLB(1).PositionColors(o,:));
end
title('Odor/Position Decoding Over Time');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([find(pfcMLB(1).beginTrialTime==0) find(pfcMLB(1).beginTrialTime==0)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1) find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

odrDecodes = cell2mat(reshape({pfcMLB.fisL1OdecodeOdr}, [1 1 length(fileDirs)]));
subplot(2,4,3);
imagesc(nanmean(odrDecodes,3), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4);
dPrms = cell2mat(cellfun(@(a)norminv(nanmean(a(logical(eye(4)))))-norminv(nanmean(a(logical(abs(eye(4)-1))))), {pfcMLB.fisL1OdecodeOdr}, 'uniformoutput',0));
title([{'Position Decoding'};{sprintf('d'' = %.02f +/- %.02f', mean(dPrms), std(dPrms)/sqrt(length(dPrms)))}]);
odrPosts = cell2mat(reshape({pfcMLB.fisL1OdecodePosts}, [1 1 length(fileDirs)]));
subplot(2,4,4);
imagesc(nanmean(odrPosts,3), [0 0.1]);
title('Posteriors');

timeTimeDecodes = cell2mat(cellfun(@(a)nanmean(a,2), {pfcMLB.fisL1OdecodeTimeTime}, 'uniformoutput', 0));
subplot(2,4,5:6);
hold on;
meanDecode = nanmean(timeTimeDecodes(:,o,:),2);
stdDecode = nanstd(timeTimeDecodes,1,2)./sqrt(sum(~isnan(timeTimeDecodes),2)-1);
% stdDecode = std(timeTimeDecodes,1,2);
plot(meanDecode, 'color', 'k', 'linewidth', 2);
patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
    'k', 'FaceAlpha', 0.75, 'EdgeColor', 'k');
tmeDecodes = cell2mat(reshape({pfcMLB.fisL1OdecodeTime}, [1 1 length(fileDirs)]));
axis tight
line(get(gca, 'xlim'), [0 0], 'linestyle', '--', 'color', 'k', 'linewidth', 1);
set(gca, 'ylim', [-1.2 1.2], 'xtick', []);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([find(pfcMLB(1).beginTrialTime==0) find(pfcMLB(1).beginTrialTime==0)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1) find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end
title('Temporal Discrepancy Over Time');
subplot(2,4,7);
imagesc(nanmean(tmeDecodes,3), [-0.1 0.1]);
colormap(gca, 'bone');
title('Average Temporal Discrepancy')
subplot(2,4,8);
tuLog = repmat(logical(triu(ones(4),1)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(tuLog(:)))*3,tmeDecodes(tuLog));
hold on;
isLog = repmat(logical(eye(4)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(isLog(:)))*2,tmeDecodes(isLog));
tlLog = repmat(logical(tril(ones(4),-1)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(tlLog(:)))*1,tmeDecodes(tlLog));
title('Discrepancy by Lag (+ or -)');

annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Separate Ordinal Decoding Figure
figure; 
sp(1) = subplot(6,1,1:5);
hold on;
for o = 1:4
    meanDecode = mean(odrTimeDecodes(:,o,:),3);
    stdDecode = std(odrTimeDecodes(:,o,:),1,3)./sqrt(size(odrTimeDecodes,3)-1);
%     stdDecode = std(odrTimeDecodes(:,o,:),1,3);
    plot(meanDecode, 'color', pfcMLB(1).PositionColors(o,:), 'linewidth', 2);
    patch([1:length(meanDecode), length(meanDecode):-1:1], [meanDecode+stdDecode; flipud(meanDecode-stdDecode)],...
        pfcMLB(1).PositionColors(o,:), 'FaceAlpha', 0.5, 'EdgeColor', pfcMLB(1).PositionColors(o,:));
end
title('Odor/Position Decoding Over Time');
axis tight
set(gca, 'ylim', [0 1], 'xtick', [], 'tickdir', 'out', 'ytick',0:0.2:1, 'yticklabel', 0:20:100);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([find(pfcMLB(1).beginTrialTime==0) find(pfcMLB(1).beginTrialTime==0)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1) find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

sp(2)= subplot(6,1,6);
plot(zeros(size(meanDecode)), '-k');
hold on;
for op = 1:4
    patch([[find(pfcMLB(1).beginTrialTime==0) find(pfcMLB(1).beginTrialTime==0)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1), [find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1) find(pfcMLB(1).endTrialTime==0)+size(pfcMLB(1).beginTrialTime,1)]+(size(pfcMLB(1).beginTrialTime,1)+size(pfcMLB(1).endTrialTime,1))*(op-1)],...
        [-0.5 0.5 0.5 -0.5], pfcMLB(1).PositionColors(op,:), 'EdgeColor', pfcMLB(1).PositionColors(op,:))
end
axis off
axis tight;
linkaxes(sp, 'x');
%% Plot FIS posteriors
pfcFISposts = nanmean(cell2mat(reshape(cellfun(@(a)nanmean(a,3), {pfcMLB.fisSeqPosts}, 'uniformoutput',0), [1 1 length(fileDirs)])),3);
figure; 
imagesc(pfcFISposts, [0 0.1]);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'w', 'linewidth', 2);  
    line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'w', 'linewidth', 2);
    line([size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'w', 'linewidth', 2);

    line(get(gca, 'xlim'), [size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, 'linestyle', '-', 'color', 'w', 'linewidth', 2);  
    line(get(gca, 'xlim'), [size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), 'linestyle', ':', 'color', 'w', 'linewidth', 2);
    line(get(gca, 'xlim'), [size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), 'linestyle', ':', 'color', 'w', 'linewidth', 2);

end
set(gca, 'ydir', 'normal')
title('Decoding Posteriors');
%%
fisOdrDecodePrd = {pfcMLB.fisL1OdecodeOdr_TrlPrd};
fisTimeDecodePrd = {pfcMLB.fisL1OdecodeTime_TrlPrd};
fisPostsPrd = {pfcMLB.fisL1OdecodePosts_TrlPrd};
figure;
for prd = 1:4
subplot(3,4,sub2ind([4,3], prd,1))
tempOdrDecode = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,prd), fisOdrDecodePrd, 'uniformoutput',0), [1 1 length(fileDirs)])),3);
imagesc(tempOdrDecode, [0 0.35]);
dPrms = cell2mat(cellfun(@(a)norminv(nanmean(a(logical(eye(4)))))-norminv(nanmean(a(logical(abs(eye(4)-1))))), cellfun(@(a)a(:,:,prd), fisOdrDecodePrd, 'uniformoutput',0), 'uniformoutput',0));
switch prd
    case 1
        tit = 'Pre-Trial';
        ylabel('Odor Decoding')
    case 2
        tit = 'Early Trial';
    case 3
        tit = 'Late Trial';
    case 4
        tit = 'Post-Trial';
end
title([{tit};{sprintf('d'' = %.02f +/- %.02f', mean(dPrms), std(dPrms)/sqrt(length(dPrms)))}]);
% title(sprintf('%s, d = %.02f', tit, dPrm));
subplot(3,4,sub2ind([4,3], prd,2))
imagesc(nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,prd), fisTimeDecodePrd, 'uniformoutput',0), [1 1 length(fileDirs)])),3), [-0.1 0.1]);
if prd==1
    ylabel('Time Decoding')
end
colormap(gca, 'bone');
subplot(3,4,sub2ind([4,3], prd,3))
imagesc(nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,prd), fisPostsPrd, 'uniformoutput',0), [1 1 length(fileDirs)])),3), [0 0.5]);
if prd==1
    ylabel('Posteriors')
end
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%%
% tic;
% bs = 10:2:200;
% corVals = nan(length(fileDirs), length(bs));
% for b = 1:length(bs)
%     for fl = 1:length(fileDirs)
%         pfcMLB(fl).dsRate = 50;
%         pfcMLB(fl).binSize = bs(b);
%         pfcMLB(fl).RunAnalysis
% %         tempDecode = pfcMLB(fl).fisL1OdecodeOdr(:);
% %         tempDecode = pfcMLB(fl).fisL1OdecodeOdr(logical(eye(4)));
%         tempDecode = pfcMLB(fl).fisL1OdecodeOdr;
% %         tempTemp = eye(4);
% %         tempTemp(tempTemp==0) = 0.001;
%         tempTemp = ones(1,4);
% %         corVals(fl,b) = sum(arrayfun(@(a,b)log(b/a), tempDecode, tempTemp(:)));
%         corVals(fl,b) = norminv(mean(tempDecode(logical(eye(4))))) - norminv(mean(tempDecode(logical(abs(eye(4)-1)))));
%     end
% %     corVals(:,b) = cell2mat(cellfun(@(a,b)corr(a(:), b(:)), {pfcMLB.fisL1OdecodeOdr}, repmat({eye(4)}, [1 length(pfcMLB)]),  'uniformoutput', 0));
%     if mod(b,10)==0
%         fprintf('%i Completed\n', b);
%     end
% end
% figure;
% plot(bs,corVals');
% legend(aniIDs);
% xlabel('Bin Size');
% ylabel('Decoding Accuracy (d'')');
% title('Odor discriminability improves with longer integration windows');