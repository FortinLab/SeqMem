fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

trialWindow = [-500 1200];
alignment = 'PokeIn';
binSize = 100;

%%
figure;
sp = nan(length(fileDirs),4);
for ani = 1:length(fileDirs)
    data = MLB_SM(fileDirs{ani});
    data.binSize = binSize;
    data.dsRate = 1;
    [spkMtx, timeMtx] = data.PP_TrialMatrix_Spiking(trialWindow, alignment);
    binarySpkMtx = spkMtx ~= 0;
    odorID = [data.trialInfo.Odor];
    posID = [data.trialInfo.Position];
    perfLog = [data.trialInfo.Performance];
    for p = 1:max(posID)
        sp(ani,p) = subplot(length(fileDirs), 4, sub2ind([4,length(fileDirs)], p,ani));
        tempISC = mean(binarySpkMtx(:,:,odorID==p & posID==p & perfLog),2);
        plot(timeMtx, mean(tempISC,3), '-k');
        hold on;
        plot(timeMtx, mean(tempISC,3)+((std(tempISC,1,3)./sqrt(size(tempISC,3)-1))*2), '--k');
        plot(timeMtx, mean(tempISC,3)-((std(tempISC,1,3)./sqrt(size(tempISC,3)-1))*2), '--k');
        drawnow;
    end
    data = [];
end
linkaxes(sp, 'xy');
for spNum = 1:numel(sp)
    plot(sp(spNum),repmat(find(timeMtx==0), [1,2]), get(gca, 'ylim'), '-k');
end