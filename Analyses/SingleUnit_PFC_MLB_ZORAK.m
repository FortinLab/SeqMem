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

%%
odorDecodeSig = cell(1,length(fileDirs));
odorDecodeNonSig = cell(1,length(fileDirs));
numUnis = nan(2,length(fileDirs));
tic
for fl = 1:length(fileDirs)
    mlb = PFC_TrialEvent_MLB_SM(fileDirs{fl});
    figure;
    sps = ceil(sqrt(length(mlb.ensembleMatrixColIDs)));
    for u = 1:length(mlb.ensembleMatrixColIDs)
        subplot(sps,sps,u);
        mlb.binSize = binSize;
        mlb.dsRate = dsRate;
        mlb.popVectIncludeLog = false(size(mlb.ensembleMatrixColIDs));
        mlb.popVectIncludeLog(u) = true;
        mlb.RunAnalysis;
        imagesc(mlb.fisL1OdecodeOdr, [0 1]);
        title(sprintf('%.02f', norminv(nanmean(mlb.fisL1OdecodeOdr(logical(eye(4)))))-norminv(nanmean(mlb.fisL1OdecodeOdr(logical(abs(eye(4)-1)))))));
        drawnow
    end
end
toc