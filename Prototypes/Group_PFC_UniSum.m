fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
binSize = 200;
dsRate = 5;

%%
tic;
for fl = 1:length(fileDirs)
    pfcUniSum(fl) = PFC_UniSum_MLB_SM(fileDirs{fl}, binSize, dsRate); %#ok<SAGROW>
%     pfcUniSum(fl).binSize = binSize;
%     pfcUniSum(fl).dsRate = dsRate;
%     pfcUniSum(fl).CalculateEpochCorrMtx;
%     pfcUniSum(fl).TrialPeriodModulation;
end
toc
%%
greatBigANOVAs = {pfcUniSum.trialEpochOdrPosPerfF};
uniSigsVals = cell2mat(cellfun(@(a)reshape(cell2mat(a(2:end-2,end,:)), [size(a,1)-3, size(a,3)]), greatBigANOVAs,'UniformOutput', 0));
uniSigs5 = uniSigsVals<0.05;
uniSigs1 = uniSigsVals<0.01;
prcntSigUnis5 = mean(uniSigs5,2);
prcntSigUnis1 = mean(uniSigs1,2);
figure;
xLabels = greatBigANOVAs{1}(2:end-2,1);
subplot(2,1,1)
bar(1:length(prcntSigUnis5), prcntSigUnis5)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xTickLabel', []);
grid on;
title('Proportion of Cells w/Significant ANOVA F-Vals (@0.05)');
subplot(2,1,2)
bar(1:length(prcntSigUnis1), prcntSigUnis1)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xticklabel', xLabels, 'xticklabelrotation', 45);
grid on;
title('Proportion of Cells w/Significant ANOVA F-Vals (@0.01)');

%%
maxRespPrd = [pfcUniSum.maxResponsePeriod];
figure; 
histogram(maxRespPrd, 'Normalization', 'probability');
set(gca, 'xtick', 1:4, 'xticklabel', [{'Pre'}, {'Early'}, {'Late'}, {'Post'}]);
title('Distribution of Cells by Max Response Periods');

prcntSigUnis5 = nan(4,size(greatBigANOVAs{1},1)-3);
prcntSigUnis1 = nan(4,size(greatBigANOVAs{1},1)-3);
for prd = 1:4
    prcntSigUnis5(prd,:) = mean(uniSigs5(:,maxRespPrd==prd),2);
    prcntSigUnis1(prd,:) = mean(uniSigs1(:,maxRespPrd==prd),2);
end
figure
xLabels = greatBigANOVAs{1}(2:end-2,1);
subplot(2,1,1)
bar(1:length(xLabels), prcntSigUnis5, 1)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xTickLabel', []);
grid on;
title('Proportion of Cells w/Significant ANOVA F-Vals (@0.05)');
subplot(2,1,2)
bar(1:length(xLabels), prcntSigUnis1, 1)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xticklabel', xLabels, 'xticklabelrotation', 45);
grid on;
title('Proportion of Cells w/Significant ANOVA F-Vals (@0.01)');
legend('Pre', 'Early', 'Late', 'Post')

%%
rewardANOVA = {pfcUniSum.trialRewardOdrPosF};
uniSigsVals = cell2mat(cellfun(@(a)reshape(cell2mat(a(2:end-2,end,:)), [size(a,1)-3, size(a,3)]), rewardANOVA,'UniformOutput', 0));
uniSigs5 = uniSigsVals<0.05;
uniSigs1 = uniSigsVals<0.01;

prcntSigUnis5 = mean(uniSigs5,2);
prcntSigUnis1 = mean(uniSigs1,2);
figure;
xLabels = rewardANOVA{1}(2:end-2,1);
subplot(2,1,1)
bar(1:length(prcntSigUnis5), prcntSigUnis5)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xTickLabel', []);
grid on;
title('Proportion of Cells w/Significant Reward ANOVA F-Vals (@0.05)');
subplot(2,1,2)
bar(1:length(prcntSigUnis1), prcntSigUnis1)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xticklabel', xLabels, 'xticklabelrotation', 45);
grid on;
title('Proportion of Cells w/Significant Reward ANOVA F-Vals (@0.01)');

% sum(isnan(uniSigsVals),2)
%%
errorANOVA = {pfcUniSum.trialErrorOdrPosF};
uniSigsVals = cell2mat(cellfun(@(a)reshape(cell2mat(a(2:end-2,end,:)), [size(a,1)-3, size(a,3)]), errorANOVA,'UniformOutput', 0));
uniSigs5 = uniSigsVals<0.05;
uniSigs1 = uniSigsVals<0.01;

prcntSigUnis5 = mean(uniSigs5,2);
prcntSigUnis1 = mean(uniSigs1,2);
figure;
xLabels = errorANOVA{1}(2:end-2,1);
subplot(2,1,1)
bar(1:length(prcntSigUnis5), prcntSigUnis5)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xTickLabel', []);
grid on;
title('Proportion of Cells w/Significant Error ANOVA F-Vals (@0.05)');
subplot(2,1,2)
bar(1:length(prcntSigUnis1), prcntSigUnis1)
set(gca, 'ylim', [0 1], 'xtick', 1:length(xLabels), 'xticklabel', xLabels, 'xticklabelrotation', 45);
grid on;
title('Proportion of Cells w/Significant Error ANOVA F-Vals (@0.01)');

% sum(isnan(uniSigsVals),2)