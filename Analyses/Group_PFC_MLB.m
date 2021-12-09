fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
binSize = 200;
dsRate = 50;

% pfcMLB = nan(size(fileDirs));
%%
for fl = 1:length(fileDirs)
%     pfcMLB(fl) = PFC_TrialEvent_MLB_SM(fileDirs{fl}, binSize, dsRate); %#ok<SAGROW>
    pfcMLB(fl) = PFC_TrialEvent_MLB_SM(fileDirs{fl});
    pfcMLB(fl).binSize = binSize;
    pfcMLB(fl).dsRate = dsRate;
%     pfcMLB(fl).beginTrialWindow = [-800 1200];
%     pfcMLB(fl).endTrialWindow = [0 0];
    pfcMLB(fl).RunAnalysis;
end
%%
%% Comment this out if running straight or just run line by line below here if loading the data
% save('PFCmlbData.mat', 'pfcMLB', '-v7.3');
%%
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
a = plot(meanPos, 'r', 'linewidth', 2);
hold on
patch([1:length(meanPos), length(meanPos):-1:1], [meanPos+stdPos; flipud(meanPos-stdPos)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
b = plot(meanRepPos, 'b', 'linewidth', 2);
patch([1:length(meanRepPos), length(meanRepPos):-1:1], [meanRepPos+stdRepPos; flipud(meanRepPos-stdRepPos)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
c = plot(meanSkpPos, 'k', 'linewidth', 2);
hold on
patch([1:length(meanSkpPos), length(meanSkpPos):-1:1], [meanSkpPos+stdSkpPos; flipud(meanSkpPos-stdSkpPos)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Position Decoding');
legend([a b c], [{'All OutSeq'}, {'Repeats'}, {'Skips'}]);
set(gca, 'ylim', [0 1]);

subplot(3,1,2)
plot(meanOdr, 'r', 'linewidth', 2);
hold on
patch([1:length(meanOdr), length(meanOdr):-1:1], [meanOdr+stdOdr; flipud(meanOdr-stdOdr)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
plot(meanRepOdr, 'b', 'linewidth', 2);
hold on
patch([1:length(meanRepOdr), length(meanRepOdr):-1:1], [meanRepOdr+stdRepOdr; flipud(meanRepOdr-stdRepOdr)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
plot(meanSkpOdr, 'k', 'linewidth', 2);
hold on
patch([1:length(meanSkpOdr), length(meanSkpOdr):-1:1], [meanSkpOdr+stdMeanOdr; flipud(meanSkpOdr-stdMeanOdr)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Odor Decoding');
set(gca, 'ylim', [0 1]);

subplot(3,1,3)
plot(meanOther, 'r', 'linewidth', 2);
hold on
patch([1:length(meanOther), length(meanOther):-1:1], [meanOther+stdOther; flipud(meanOther-stdOther)],...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r');
plot(meanRepOther, 'b', 'linewidth', 2);
hold on
patch([1:length(meanRepOther), length(meanRepOther):-1:1], [meanRepOther+stdRepOther; flipud(meanRepOther-stdRepOther)],...
    'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b');
plot(meanSkpOther, 'k', 'linewidth', 2);
hold on
patch([1:length(meanSkpOther), length(meanSkpOther):-1:1], [meanSkpOther+stdSkpOther; flipud(meanSkpOther-stdSkpOther)],...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('Other Decoding');
set(gca, 'ylim', [0 1]);


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
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

odrDecodes = cell2mat(reshape({pfcMLB.fisL1OdecodeOdr}, [1 1 length(fileDirs)]));
subplot(2,4,3);
imagesc(nanmean(odrDecodes,3), [0 0.5]);
odrPosts = cell2mat(reshape({pfcMLB.fisL1OdecodePosts}, [1 1 length(fileDirs)]));
subplot(2,4,4);
imagesc(nanmean(odrPosts,3), [0 0.5]);

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
    line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end
subplot(2,4,7);
imagesc(nanmean(tmeDecodes,3), [-0.1 0.1]);
colormap(gca, 'bone');
subplot(2,4,8);
tuLog = repmat(logical(triu(ones(4),1)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(tuLog(:)))*3,tmeDecodes(tuLog));
hold on;
isLog = repmat(logical(eye(4)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(isLog(:)))*2,tmeDecodes(isLog));
tlLog = repmat(logical(tril(ones(4),-1)), [1 1 length(fileDirs)]);
swarmchart(ones(1,sum(tlLog(:)))*1,tmeDecodes(tlLog));

annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Plot FIS posteriors
pfcFISposts = nanmean(cell2mat(reshape(cellfun(@(a)nanmean(a,3), {pfcMLB.fisSeqPosts}, 'uniformoutput',0), [1 1 length(fileDirs)])),3);
figure; 
load('batlowW.mat');
colormap(batlowW);
z = colormap;
z = flipud(z);
% timeTicks = pfcMLB(1).beginTrialWindow(1):pfcMLB(1).dsRate:(diff(pfcMLB(1).beginTrialWindow)+diff(pfcMLB(1).endTrialWindow)*4);
% xLbl = repmat([pfcMLB(1).beginTrialTime; pfcMLB(1).endTrialTime], [4,1]);
% reshape(find(xLbl==0), [2,length(find(xLbl==0))/2]);
imagesc(pfcFISposts, [0 0.1]);
for op = 1:4
    line([size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 1);
    line([size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 1);

    line(get(gca, 'xlim'), [size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1), size(pfcMLB(1).beginTrialMtx,1) + size(pfcMLB(1).endTrialMtx,1)]*op, 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line(get(gca, 'xlim'), [size(pfcMLB(1).beginTrialMtx,1)/2 size(pfcMLB(1).beginTrialMtx,1)/2]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), 'linestyle', ':', 'color', 'k', 'linewidth', 1);
    line(get(gca, 'xlim'), [size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1) size(pfcMLB(1).endTrialMtx,1)/2+size(pfcMLB(1).beginTrialMtx,1)]+(size(pfcMLB(1).beginTrialMtx,1)+size(pfcMLB(1).endTrialMtx,1))*(op-1), 'linestyle', ':', 'color', 'k', 'linewidth', 1);

end
colormap(z);
set(gca, 'ydir', 'normal', 'xtick', [], 'ytick', [])
%%
fisOdrDecodePrd = {pfcMLB.fisL1OdecodeOdr_TrlPrd};
fisTimeDecodePrd = {pfcMLB.fisL1OdecodeTime_TrlPrd};
fisPostsPrd = {pfcMLB.fisL1OdecodePosts_TrlPrd};
figure;
for prd = 1:4
subplot(3,4,sub2ind([4,3], prd,1))
tempOdrDecode = nanmean(cell2mat(reshape(cellfun(@(a)a(:,:,prd), fisOdrDecodePrd, 'uniformoutput',0), [1 1 length(fileDirs)])),3);
imagesc(tempOdrDecode, [0 0.35]);
dPrm = norminv(nanmean(tempOdrDecode(logical(eye(4))))) - norminv(nanmean(tempOdrDecode(logical(abs(eye(4)-1)))));
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
title(sprintf('%s, d = %.02f', tit, dPrm));
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
% plot(bs,corVals')
    