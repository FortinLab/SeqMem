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

binSize = 200;
dsRate = 1;
trlWindow = {[-800 1500]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};

frThresh = 1/(binSize/1000);
% frThresh = 0.001;
%% Data Vectors
popVects = cell(length(fileDirs),1);
popVectsSortVect = cell(length(fileDirs),1);
popVectsThreshVect = cell(length(fileDirs),1);
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.binType = 'gauss';
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
        end
        nsmblSpks(:,:,op) = tempNsmbl';
    end
    popVects{ani} = nsmblSpks;
    popVectsSortVect{ani} = peakNdx;
    popVectsThreshVect{ani} = peakRt;
end

%% Plot PopVects
grpPV = cell2mat(popVects);
% grpPVsort = timeVect(cell2mat(popVectsSortVect));
grpPVsort = cell2mat(popVectsSortVect);
% grpPVmaxFR(sum(grpPVmaxFR>=frThresh,2)>=1,:) = 100;
% grpPVmaxFR(grpPVmaxFR~=100)=0;
grpPVmaxFR = cell2mat(popVectsThreshVect);

figure;
ndxCorr = nan(mlb.seqLength);
pvCorr = nan(mlb.seqLength);
pvCorrTime = nan(mlb.seqLength, mlb.seqLength, length(timeVect));
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        sortedPV = sortrows([grpPVsort(:,op1), grpPVmaxFR(:,op2), grpPV(:,:,op2)]);
        subplot(mlb.seqLength,mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
        imagesc(timeVect, 1:sum(sortedPV(:,2)>=frThresh), sortedPV(sortedPV(:,2)>=frThresh,3:end), [0 1]);
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
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Population Vectors Sortings:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');


figure; 
subplot(1,2,1);
imagesc(ndxCorr, [0 1]);
title('Max FR Index Similarity');
subplot(1,2,2);
imagesc(pvCorr, [0 1]);
title('Population Vector Similarity');
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Trial Position Correlations:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

figure;
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
        plot(timeVect, reshape(pvCorrTime(op1,op2,:), [1,size(pvCorrTime,3)]));
        axis tight;
        set(gca, 'ylim', [0 1]);
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("PV Correlation Across Time:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

figure;
for op1 = 1:mlb.seqLength
    for op2 = 1:mlb.seqLength
        subplot(mlb.seqLength, mlb.seqLength, sub2ind([mlb.seqLength, mlb.seqLength], op2, op1));
%         corrScatPlot(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op1), grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op2), op1,op2,[]);
        [tempCounts, xedge,yedge] = histcounts2(timeVect(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op1)), timeVect(grpPVsort(grpPVmaxFR(:,op2)>=frThresh,op2)),downsample(timeVect,binSize/dsRate),downsample(timeVect,binSize/dsRate));
        for t = 1:size(tempCounts,1)
            tempCounts(:,t) = tempCounts(:,t)/sum(tempCounts(:,t));
        end
        imagesc(xedge(2:end)-mode(diff(xedge))/2, yedge(2:end)-mode(diff(yedge))/2, tempCounts, [0 0.5]);
        set(gca, 'ydir','normal');
    end
end
annotation(gcf,'textbox', [0.05 0.95 0.8 0.05],...
    'String', sprintf("Temporal Remapping:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');


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
    'String', sprintf("Rate vs Time Remapping:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');

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
    'String', sprintf("Rate vs Temporal Remapping:\n     binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), FR Threshold = %.02f spk/s",...
    binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, frThresh),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');



    