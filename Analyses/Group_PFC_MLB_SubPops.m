% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
binSize = 200;
dsRate = 5;
sigThresh = 0.05;
% selContrast = 'Position';
% selContrast = 'Epoch*Position';
selContrast = 'Odor';
% selContrast = 'Position*Odor';
% selContrast = 'Performance';

%%
odorDecodeSig = cell(1,length(fileDirs));
odorDecodeNonSig = cell(1,length(fileDirs));
numUnis = nan(2,length(fileDirs));
tic
for fl = 1:length(fileDirs)
    tempUniSum = PFC_UniSum_MLB_SM(fileDirs{fl}, binSize, dsRate);
    tempSigValMatrix = reshape(cell2mat(tempUniSum.trialEpochOdrPosPerfF(2:end-2,end,:)), [size(tempUniSum.trialEpochOdrPosPerfF,1)-3, size(tempUniSum.trialEpochOdrPosPerfF,3)]);
    contrasts = tempUniSum.trialEpochOdrPosPerfF(2:end-2,1);
    rowLog = strcmp(contrasts, selContrast);
    mlb = PFC_TrialEvent_MLB_SM(fileDirs{fl});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.popVectIncludeLog = tempSigValMatrix(rowLog,:)<sigThresh;
    numUnis(1,fl) = sum(mlb.popVectIncludeLog);
    if  numUnis(1,fl) ~= 0
        numUnis(1,fl) = sum(mlb.popVectIncludeLog);
        mlb.RunAnalysis;
        odorDecodeSig{fl} = mlb.fisL1OdecodeOdr;
    else
        odorDecodeSig{fl} = nan;
    end
    numUnis(2,fl) = sum(~mlb.popVectIncludeLog);
    if numUnis(2,fl) ~= 0    
        mlb.popVectIncludeLog = tempSigValMatrix(rowLog,:)>sigThresh;
        mlb.RunAnalysis;
        odorDecodeNonSig{fl} = mlb.fisL1OdecodeOdr;
    else
        odorDecodeNonSig{fl} = nan;
    end
    clear tempUniSum mlb
end
toc

%%
figure;
subPopDprm = nan(1,length(fileDirs));
for sp = 1:5
    subplot(2,5,sp)
    imagesc(odorDecodeSig{sp}, [0 0.5]);
    subPopDprm(sp) = norminv(nanmean(odorDecodeSig{sp}(logical(eye(4))))) - norminv(nanmean(odorDecodeSig{sp}(logical(abs(eye(4)-1)))));
    title([{sprintf('%i/%i', numUnis(1,sp), sum(numUnis(:,sp)))};...
        {sprintf('d'' = %.02f', subPopDprm(sp))}]);
end
nonSubPopDprm = nan(1,length(fileDirs));
for sp = 6:10
    subplot(2,5,sp)
    imagesc(odorDecodeNonSig{sp-5}, [0 0.5]);
    nonSubPopDprm(sp) = norminv(nanmean(odorDecodeNonSig{sp-5}(logical(eye(4))))) - norminv(nanmean(odorDecodeNonSig{sp-5}(logical(abs(eye(4)-1)))));
    title(sprintf('d'' = %.02f', nonSubPopDprm(sp)));
end

annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, Sig = %.02f, Contrast = %s', binSize, dsRate, sigThresh, selContrast),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

[h,p,ci,stats] = ttest2(subPopDprm, nonSubPopDprm);
annotation(gcf,'textbox', [0 0.5 1 0.05],'String', sprintf('t(%i)=%.02f, p=%.05f',stats.df, stats.tstat, p),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
