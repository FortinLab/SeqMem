% PFC_SingleUnit_Modulation
%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];

% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

%% Colormap setup
% cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);
% cMap2 = load('batlow.mat');
% cMap2 = load('batlowK.mat');
% cMap2 = load('batlowW.mat');
% cMap2 = load('davos.mat');
% cMap2 = load('grayC.mat');
cMap2 = load('roma.mat'); %%
cMap2 = cMap2.(cell2mat(fieldnames(cMap2)));
cMap2 = flipud(cMap2);
%% Compile SU Objects
su = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    su{ani} = SingleUnit_SM(fileDirs{ani});
end

%% Plot Single Unit Summaries
for ani = 1:length(fileDirs)
    su{ani}.gaussWinDur = 200;
    su{ani}.binSize = 200;
    su{ani}.dsRate = 10;
%     su{ani}.savePlots = true;
    su{ani}.PlotUnitSummary;
    su{ani}.savePlots = false;
    close all
end

%% Plot Position Ensembles organized by different positions
frThresh = 0.05;
window = [-1200 2000];
meanPopVecPos = cell(1,length(fileDirs));
iscPokeDur = cell(1,length(fileDirs));
iscRwdLat = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    su{ani}.gaussWinDur = 200;
    su{ani}.gaussWinWidth = 2.5;
    su{ani}.binSize = 200;
            
    isLog = [su{ani}.trialInfo.TranspositionDistance]==0;
    perfLog = [su{ani}.trialInfo.Performance]==1;
    iscPokeDur{ani} = [su{ani}.trialInfo(isLog & perfLog).PokeDuration];
    iscRwdLat{ani} = [su{ani}.trialInfo(isLog & perfLog).RewardIndex] - [su{ani}.trialInfo(isLog & perfLog).PokeInIndex];
        
    eventSpikes = su{ani}.ExtractTrialMatrix(su{ani}.ensembleMatrix, [window(1)-su{ani}.gaussWinDur/2, window(2)+su{ani}.gaussWinDur/2], 'PokeIn');
    [evtSpks, trlIDs] = su{ani}.ExtractTrialSpikes(eventSpikes, 'isc');
    
    instFRgauss = gausswin(su{ani}.gaussWinDur,su{ani}.gaussWinWidth);
    instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(su{ani}.tsVect)));
    
    instFR = nan(size(evtSpks));
    for trl = 1:size(instFR,3)
        for uni = 1:size(instFR,2)
            instFR(:,uni,trl) = conv(evtSpks(:,uni,trl), instFRgauss, 'same');
        end
    end
    instFR = instFR(su{ani}.gaussWinDur/2+1:size(instFR,1)-(su{ani}.gaussWinDur/2),:,:);
    tempPopV = nan(size(instFR,1), size(instFR,2), su{ani}.seqLength);
    for pos = 1:su{ani}.seqLength
        curPosTrls = instFR(:,:,trlIDs(2,:)==pos);
        tempPopV(:,:,pos) = mean(curPosTrls,3,'omitnan');
    end
    meanPopVecPos{ani} = tempPopV;
end
meanPopVecPos = permute(cell2mat(meanPopVecPos), [2,1,3]);
normPopVecPos = nan(size(meanPopVecPos));
maxFR = nan(size(meanPopVecPos,1), size(meanPopVecPos,3));
peakNdx = nan(size(meanPopVecPos,1), size(meanPopVecPos,3));
for pos = 1:size(meanPopVecPos,3)
    for uni = 1:size(meanPopVecPos,1)
        normPopVecPos(uni,:,pos) = meanPopVecPos(uni,:,pos)./max(meanPopVecPos(uni,:,pos));
        maxFR(uni,pos) = max(meanPopVecPos(uni,:,pos));
        peakNdx(uni,pos) = find(meanPopVecPos(uni,:,pos)==maxFR(uni,pos), 1, 'first');
    end
end

figure;
for pos = 1:su{ani}.seqLength
    curVect = normPopVecPos(:,:,pos);
    for sort = 1:su{ani}.seqLength
        sortPopV = sortrows([peakNdx(:,sort),maxFR(:,pos), curVect]);
        subplot(su{ani}.seqLength, su{ani}.seqLength, sub2ind(repmat(su{ani}.seqLength,[1,2]), sort, pos));        
        imagesc(window(1):window(2), 1:sum(sortPopV(:,2)>frThresh), sortPopV(sortPopV(:,2)>frThresh,2:end), [0 1]);
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(repmat(median(cell2mat(iscPokeDur))*1000,[1,2]), get(gca, 'ylim'), '-k');
        plot(repmat(median(cell2mat(iscRwdLat)),[1,2]), get(gca, 'ylim'), ':k');
    end
end
colormap(cMap2);
colorbar;

%% Evaluate Position Info coding w/in ensemble
window = [-1200 2000];
posInfo = cell(1,length(fileDirs));
posInfoSANSA = cell(1,length(fileDirs));
iscSpks = cell(1,length(fileDirs));
iscSpksSANSA = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    su{ani}.dsRate = 10;
    su{ani}.binSize = 200;
    [posInfo{ani}, tsVect, tempSpks, tempTrlIDs] = su{ani}.QuantPosInfo(window, 'PokeIn', 'pos', 'corr');
    posInfoSANSA{ani} = su{ani}.QuantPosInfo(window, 'PokeIn', 'pos', 'corrSANSA');
    
    iscSpks{ani} = mean(tempSpks(:,:,tempTrlIDs(5,:)==1),3);
    iscSpksSANSA{ani} = mean(tempSpks(:,:,tempTrlIDs(5,:)==1 & tempTrlIDs(2,:)~=1),3);
end
%
ensmblPosInfo = cell2mat(posInfo);
ensmblPosInfoSANSA = cell2mat(posInfoSANSA);
ensmblSpks = cell2mat(iscSpks);
ensmblSpksSANSA = cell2mat(iscSpksSANSA);
normPI = nan(size(ensmblPosInfo));
loc = nan(size(ensmblPosInfo,2),1);
normPIsansa = nan(size(ensmblPosInfo));
locSANSA = nan(size(ensmblPosInfo,2),1);
uniCorSimAsansA = nan(size(ensmblPosInfo,2),1);
spkNorm = nan(size(ensmblPosInfo));
spkLoc = nan(size(ensmblPosInfo,2),1);
spkNormSANSA = nan(size(ensmblPosInfo));
spkLocSANSA = nan(size(ensmblPosInfo,2),1);
for u = 1:size(ensmblPosInfo,2)
    normPI(:,u) = ensmblPosInfo(:,u)./max(ensmblPosInfo(:,u));
    loc(u) = find(normPI(:,u)==1,1,'first');
    normPIsansa(:,u) = ensmblPosInfoSANSA(:,u)./max(ensmblPosInfoSANSA(:,u));
    locSANSA(u) = find(normPIsansa(:,u)==1,1,'first');
    spkNorm(:,u) = ensmblSpks(:,u)./max(ensmblSpks(:,u));
    spkLoc(u) = find(spkNorm(:,u)==1,1,'first');
    spkNormSANSA(:,u) = ensmblSpksSANSA(:,u)./max(ensmblSpksSANSA(:,u));
    if sum(spkNormSANSA(:,u))~=1 && ~isnan(sum(spkNormSANSA(:,u)))
        spkLocSANSA(u) = find(spkNormSANSA(:,u)==1,1,'first');
    end
    nanLog = sum(isnan([ensmblPosInfo(:,u), ensmblPosInfoSANSA(:,u)]),2)==0;
    uniCorSimAsansA(u) = corr(ensmblPosInfo(nanLog,u), ensmblPosInfoSANSA(nanLog,u));
end

sortNormPInorm = sortrows([loc, normPI']);
sortNormPInorm = sortNormPInorm(:,2:end);
sortNormPIsansa = sortrows([locSANSA, normPI']);
sortNormPIsansa = sortNormPIsansa(:,2:end);
sortNormPIspks = sortrows([spkLoc, normPI']);
sortNormPIspks = sortNormPIspks(:,2:end);
sortNormPIspksSANSA = sortrows([spkLocSANSA, normPI']);
sortNormPIspksSANSA = sortNormPIspksSANSA(:,2:end);
sortSANSAnormPI = sortrows([loc, normPIsansa']);
sortSANSAnormPI = sortSANSAnormPI(:,2:end);
sortSANSAsansa = sortrows([locSANSA, normPIsansa']);
sortSANSAsansa = sortSANSAsansa(:,2:end);
sortSANSAspks = sortrows([spkLoc, normPIsansa']);
sortSANSAspks = sortSANSAspks(:,2:end);
sortSANSAspksSANSA = sortrows([spkLocSANSA, normPIsansa']);
sortSANSAspksSANSA = sortSANSAspksSANSA(:,2:end);
sortSPKnorm = sortrows([loc, spkNorm']);
sortSPKnorm = sortSPKnorm(:,2:end);
sortSPKsansa = sortrows([locSANSA, spkNorm']);
sortSPKsansa = sortSPKsansa(:,2:end);
sortSPKspk = sortrows([spkLoc, spkNorm']);
sortSPKspk = sortSPKspk(:,2:end);
sortSPKspkSANSA = sortrows([spkLocSANSA, spkNorm']);
sortSPKspkSANSA = sortSPKspkSANSA(:,2:end);
sortSPKsansaNORM = sortrows([loc, spkNormSANSA']);
sortSPKsansaNORM = sortSPKsansaNORM(:,2:end);
sortSPKsansaSANSA = sortrows([locSANSA, spkNormSANSA']);
sortSPKsansaSANSA = sortSPKsansaSANSA(:,2:end);
sortSPKsansaSPK = sortrows([spkLoc, spkNormSANSA']);
sortSPKsansaSPK = sortSPKsansaSPK(:,2:end);
sortSPKsansaSPKsansa = sortrows([spkLocSANSA, spkNormSANSA']);
sortSPKsansaSPKsansa = sortSPKsansaSPKsansa(:,2:end);

figure; 
oiSP(1) = subplot(4,4,1);
imagesc(tsVect, 1:size(sortNormPInorm,1), sortNormPInorm, [0 1]);
title('All Trials OrdInfo (sort by all trials)');
oiSP(2) = subplot(4,4,2);
imagesc(tsVect, 1:size(sortNormPIsansa,1), sortNormPIsansa, [0 1]);
title('All Trials OrdInfo (sort by SANS1)');
oiSP(3) = subplot(4,4,3);
imagesc(tsVect, 1:size(sortNormPIspks,1), sortNormPIspks, [0 1]);
title('All Trials (sort by mean SpkR)');
oiSP(4) = subplot(4,4,4);
imagesc(tsVect, 1:size(sortNormPIspksSANSA,1), sortNormPIspksSANSA, [0 1]);
title('All Trials (sort by mean SpkR SANS1)');

oiSP(5) = subplot(4,4,5);
imagesc(tsVect, 1:size(sortSANSAnormPI,1), sortSANSAnormPI, [0 1]);
title('SANS1 (sort by all)');
oiSP(6) = subplot(4,4,6);
imagesc(tsVect, 1:size(sortSANSAsansa,1), sortNormPInorm, [0 1]);
title('SANS1 (sort by SANS1)');
oiSP(7) = subplot(4,4,7);
imagesc(tsVect, 1:size(sortSANSAspks,1), sortSANSAspks, [0 1]);
title('SANS1 (sort by mean SpkR)');
oiSP(8) = subplot(4,4,8);
imagesc(tsVect, 1:size(sortSANSAspksSANSA,1), sortSANSAspksSANSA, [0 1]);
title('SANS1 (sort by mean SpkR SANS1)');

srSP(1) = subplot(4,4,9);
imagesc(tsVect, 1:size(sortSPKnorm,1), sortSPKnorm, [0 1]);
title('Mean SpkR (sort by norm PI)');
srSP(2) = subplot(4,4,10);
imagesc(tsVect, 1:size(sortSPKsansa,1), sortSPKsansa, [0 1]);
title('Mean SpkR (sort by SANS1)');
srSP(3) = subplot(4,4,11);
imagesc(tsVect, 1:size(sortSPKspk,1), sortSPKspk, [0 1]);
title('Mean SpkR (sort by mean SpkR)');
srSP(4) = subplot(4,4,12);
imagesc(tsVect, 1:size(sortSPKspkSANSA,1), sortSPKspkSANSA, [0 1]);
title('Mean SpkR (sort by mean SpkR SANS1)');

srSP(5) = subplot(4,4,13);
imagesc(tsVect, 1:size(sortSPKsansaNORM,1), sortSPKsansaNORM, [0 1]);
title('Mean SpkR SANS1 (sort by PI)');
srSP(6) = subplot(4,4,14);
imagesc(tsVect, 1:size(sortSPKsansaSANSA,1), sortSPKsansaSANSA, [0 1]);
title('Mean SpkR SANS1 (sort by PI SANS1)');
srSP(7) = subplot(4,4,15);
imagesc(tsVect, 1:size(sortSPKsansaSPK,1), sortSPKsansaSPK, [0 1]);
title('Mean SpkR SANS1 (sort by mean SpkR)');
srSP(8) = subplot(4,4,16);
imagesc(tsVect, 1:size(sortSPKsansaSPKsansa,1), sortSPKsansaSPKsansa, [0 1]);
title('Mean SpkR SANS1 (sort by mean SpkR SANS1)');

arrayfun(@(a) colormap(a, cMap2), oiSP);
% arrayfun(@(a) colormap(a, cMap), srSP);

figure; 
subplot(2,1,1)
scatter(tsVect(loc), uniCorSimAsansA);
subplot(2,1,2)
scatter(tsVect(locSANSA), uniCorSimAsansA);
%% Evaluate Properties of Firing Rate change during Ordinal Info Peaks
timeBins = -1200:200:1600;
peakBins = 0:5:50;
widthBins = 0:10:100;
corrBins = -1:0.1:1;
polyRatBins = 0.95:0.05:1.5;

peakInfo = [];
for ani = 1:length(fileDirs)
%     peakInfo = [peakInfo, su{ani}.InterrogatePeakInfo([min(timeBins) max(timeBins)], 'PokeIn', 'pos', 'iscSANSA', 3.3)];  %#ok<AGROW>
    peakInfo = [peakInfo, su{ani}.InterrogatePeakInfo([min(timeBins) max(timeBins)], 'PokeIn', 'pos', 'isc', 3.3)];  %#ok<AGROW>
end
locs = cell2mat({peakInfo.Locations})';
peaks = cell2mat({peakInfo.Peaks}');
widths = cell2mat({peakInfo.Widths}');
proms = cell2mat({peakInfo.Prominance}');
fits = cell2mat({peakInfo.LinearFit}');
polyRat = cell2mat({peakInfo.PolyRatio}');

tallPeakLocs = cell2mat(arrayfun(@(a){a.Locations(find(a.Peaks==max(a.Peaks),1,'first'))},peakInfo));
tallPeakPeaks = cell2mat(arrayfun(@(a){a.Peaks(find(a.Peaks==max(a.Peaks),1,'first'))},peakInfo)');
tallPeakWidths = cell2mat(arrayfun(@(a){a.Widths(find(a.Peaks==max(a.Peaks),1,'first'))},peakInfo)');
tallPeakFits = cell2mat(arrayfun(@(a){a.LinearFit(find(a.Peaks==max(a.Peaks),1,'first'))},peakInfo)');
tallPeakPolyRat = cell2mat(arrayfun(@(a){a.PolyRatio(find(a.Peaks==max(a.Peaks),1,'first'))},peakInfo)');
%
sp(1) = subplot(4,1,1);
imagesc(timeBins(1:end-1)+(mode(diff(timeBins))/2), peakBins(1:end-1)+(mode(diff(peakBins))/2), histcounts2(locs, peaks, timeBins, peakBins)');
hold on;
scatter(locs, peaks, 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.25);
title('F-Values');
sp(2) = subplot(4,1,2);
imagesc(timeBins(1:end-1)+(mode(diff(timeBins))/2), widthBins(1:end-1)+(mode(diff(widthBins))/2), histcounts2(locs, widths, timeBins, widthBins)');
hold on;
scatter(locs, widths, 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.25);
title('Peak Width');
sp(3) = subplot(4,1,3);
imagesc(timeBins(1:end-1)+(mode(diff(timeBins))/2), corrBins(1:end-1)+(mode(diff(corrBins))/2), histcounts2(locs, fits, timeBins, corrBins)');
hold on;
scatter(locs, fits, 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.25);
title('Linear Fit');
sp(4) = subplot(4,1,4);
imagesc(timeBins(1:end-1)+(mode(diff(timeBins))/2), polyRatBins(1:end-1)+(mode(diff(polyRatBins))/2), histcounts2(locs, polyRat , timeBins, polyRatBins)');
hold on;
scatter(locs, polyRat , 'markeredgecolor', 'none', 'markerfacecolor', 'k', 'markerfacealpha', 0.25);
title('Poly Ratio (MSE-poly1/MSE-poly2)');

%%  Firing Rate Change over time across time
winStepSize = 50;
horizon = 10;
trlType = 'isc';
% trlType = 'iscSANSA';
figure;
spkDiffTime = cell(1,1,length(fileDirs));
for ani = 1:length(fileDirs)
%     [spkDiffTime{ani},tsVect] = su{ani}.SpikeDiffAcrossTime([-1200 2000],'PokeIn',trlType,winStepSize,horizon);
    [spkDiffTime{ani},tsVect] = su{ani}.SpikeDiffAcrossTime([-2000 1500],'PokeOut',trlType,winStepSize,horizon);
    iscPokeDur = [su{ani}.trialInfo([su{ani}.trialInfo.TranspositionDistance]==0 & [su{ani}.trialInfo.Performance]==1).PokeDuration];
    windowEnd = floor((mean(iscPokeDur))*100)/100;
    subplot(1,length(fileDirs)+1,ani);
    imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect,mean(spkDiffTime{ani},3,'omitnan'),[-0.5 0.5]);
%     imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect,abs(mean(spkDiffTime{ani},3,'omitnan')),[0 0.5]);
    hold on; 
%     plot(get(gca, 'xlim'), repmat(windowEnd*1000,[1,2]), '-k', 'linewidth', 2);
    plot(get(gca, 'xlim'), repmat(windowEnd*-1000,[1,2]), '-k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
end 
grpSpkDiffTime = cell2mat(spkDiffTime);
subplot(1,length(fileDirs)+1, ani+1);
imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect, mean(grpSpkDiffTime,3), [-0.5 0.5]);
set(gca, 'ydir', 'normal');
colormap(cMap);
figure;
imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect, mean(abs(grpSpkDiffTime),3), [0 2]);
set(gca, 'ydir', 'normal');
% colormap(cMap);

%% Ensemble similarity over time across time
winStepSize = 50;
horizon = 10;
trlType = 'isc';
% trlType = 'iscSANSA';
figure;
stateDiffTime = cell(1,1,length(fileDirs));
for ani = 1:length(fileDirs)
    [stateDiffTime{ani},tsVect] = su{ani}.EnsembleStateDiffAcrossTime([-1200 2000],'PokeIn',trlType,winStepSize,horizon);
%     [stateDiffTime{ani},tsVect] = su{ani}.EnsembleStateDiffAcrossTime([-2000 1500],'PokeOut',trlType,winStepSize,horizon);
    iscPokeDur = [su{ani}.trialInfo([su{ani}.trialInfo.TranspositionDistance]==0 & [su{ani}.trialInfo.Performance]==1).PokeDuration];
    windowEnd = floor((mean(iscPokeDur))*100)/100;
    subplot(1,length(fileDirs)+1,ani);
    imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect,stateDiffTime{ani});
    hold on; 
    plot(get(gca, 'xlim'), repmat(windowEnd*1000,[1,2]), '-k', 'linewidth', 2);
%     plot(get(gca, 'xlim'), repmat(windowEnd*-1000,[1,2]), '-k', 'linewidth', 2);
    set(gca, 'ydir', 'normal');
end 
subplot(1,length(fileDirs)+1, ani+1);
imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect, mean(cell2mat(stateDiffTime),3));
set(gca, 'ydir', 'normal');
colormap(cMap);

figure; 
imagesc(winStepSize:winStepSize:(winStepSize*horizon),tsVect, mean(cell2mat(stateDiffTime),3));
set(gca, 'ydir', 'normal');
colormap(cMap);


%% Pos Info vs Spike Rate
aniPosSpkCorr = cell(length(fileDirs), 1);
aniUniPrd = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    iscPokeDur = [su{ani}.trialInfo([su{ani}.trialInfo.TranspositionDistance]==0 & [su{ani}.trialInfo.Performance]==1).PokeDuration];
%     windowEnd = floor((mean(iscPokeDur)+0.6)*100)/100;
    windowEnd = floor((mean(iscPokeDur))*100)/100;
    su{ani}.binSize = 200;
    su{ani}.dsRate = 50;
    [aniPosSpkCorr{ani}, posInfo, varSpks, tsVect] = su{ani}.QuantPosSpkCorr([-1200 windowEnd*1000],'PokeIn','pos','iscSANSA');
    posSpkNdxs = nan(2,size(varSpks,2));
    for u = 2:size(varSpks,2)
        if sum(isnan(posInfo(:,u)))~=size(posInfo,1)
            posSpkNdxs(1,u) = tsVect(find(posInfo(:,u)==max(posInfo(:,u)),1,'first'));
            posSpkNdxs(2,u) = tsVect(find(varSpks(:,u,end)==max(varSpks(:,u,end)),1,'first'));
        end
    end
    aniUniPrd{ani} = posSpkNdxs;
end

corrVals = cell2mat(aniPosSpkCorr)-1;
peakNdxs = cell2mat(aniUniPrd);
%
figure; 
subplot(1,2,1)
su{ani}.PlotMeanVarSwarmBar(1,corrVals(peakNdxs(1,:)<0,end),1,0.05,'k');
su{ani}.PlotMeanVarSwarmBar(2,corrVals(peakNdxs(1,:)>=0,end),1,0.05,'r');
subplot(1,2,2)
su{ani}.PlotMeanVarSwarmBar(1,corrVals(peakNdxs(2,:)<0,end),1,0.05,'k');
su{ani}.PlotMeanVarSwarmBar(2,corrVals(peakNdxs(2,:)>=0,end),1,0.05,'r');

[h,p,ci,stats] = ttest2(corrVals(peakNdxs(2,:)<0,end), corrVals(peakNdxs(2,:)>0,end))

figure;
for p = 1:size(corrVals,2)
    su{ani}.PlotMeanVarSwarmBar(p,corrVals(:,p),1,0.05,'k');
end
%
figure;
subplot(2,1,1)
% scatter(peakNdxs(1,:),corrVals(:,end), 'k', 'filled', 'markerfacealpha', 0.5);
[~,~,scatPlot] = corrScatPlot(peakNdxs(1,:)',corrVals(:,end));
set(scatPlot, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
xlabel('Peak SpkR Index');
grid on;
subplot(2,1,2)
% scatter(peakNdxs(2,:),corrVals(:,end), 'k', 'filled', 'markerfacealpha', 0.5);
[~,~,scatPlot] = corrScatPlot(peakNdxs(2,:)',corrVals(:,end));
set(scatPlot, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
xlabel('Peak OrdInfo Index');
grid on;
linkaxes;
        
%% Analyses of Session Split for "time cell" plots
% "Time Cell" Plots
% - Ensemble Peak Latency Correlations
%   - By Position 
%   - Mean Response
% - Distribution of Individual cell mean SpkR/Time correlations
%
% Compare the ensemble and individual cell distributions w/in (session split) vs across positions
%   ... If temporal structure is a "backbone" then it should be fairly similar... or scale a bit with position
%   ... Regardless... the w/in position should be higher than across positions (or similar, it shouldn't be worse but
%   the power loss may kill it)



%%
% trlTAB = cell(1,length(fileDirs));
% errTAB = cell(1,length(fileDirs));
% rwdTAB = cell(1,length(fileDirs));
% trlSTATS = cell(1,length(fileDirs));
% errSTATS = cell(1,length(fileDirs));
% rwdSTATS = cell(1,length(fileDirs));
% trlDATA = cell(1,length(fileDirs));
% errDATA = cell(1,length(fileDirs));
% rwdDATA = cell(1,length(fileDirs));
% trlIDS = cell(1,length(fileDirs));
% errIDS = cell(1,length(fileDirs));
% rwdIDS = cell(1,length(fileDirs));
% errT = cell(1,length(fileDirs));
% errTABseq = cell(1,length(fileDirs));
% errSTATSseq = cell(1,length(fileDirs));
% for ani = 1:length(fileDirs)   
%     [trlTAB{ani}, trlSTATS{ani}, trlDATA{ani}, trlIDS{ani}] = su{ani}.TrialPeriodSpiking;
%     [errTAB{ani}, errSTATS{ani}, errDATA{ani}, errIDS{ani}] = su{ani}.ErrorSpiking;
%     [rwdTAB{ani}, rwdSTATS{ani}, rwdDATA{ani}, rwdIDS{ani}] = su{ani}.RewardSpiking;
%     [errT{ani}, errTABseq{ani}, errSTATSseq{ani}] = su{ani}.ErrorSpikingSequential;
% end
%%
% trlPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},trlTAB));
% errPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},errTAB));
% rwdPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},rwdTAB));
% errTPvals = cell2mat(cellfun(@(a){[a.p]}, errT));
% errSeqPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2,end,:)))},errTABseq)')';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% DEAD ENDS %%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate Inter-trial vs Trial Correlations (dead end)
% trlSpks = cell(length(fileDirs),1);
% itdSpks = cell(length(fileDirs),1);
% itdTrlCorr = cell(1,length(fileDirs));
% for ani = 1:length(fileDirs)
%     iscPokeDur = [su{ani}.trialInfo([su{ani}.trialInfo.TranspositionDistance]==0 & [su{ani}.trialInfo.Performance]==1).PokeDuration];
%     windowEnd = floor((mean(iscPokeDur))*100)/100;
%     trlSpksBin = su{ani}.BinTrialEventSpikes([0 windowEnd*1000], 'PokeIn');
%     [tempTrlSpks, trlIDs] = su{ani}.ExtractTrialSpikes(squeeze(trlSpksBin./windowEnd), 'isc');
%     trlSpks{ani} = tempTrlSpks;
%     itdSpksBin = su{ani}.BinTrialEventSpikes([-800 0], 'PokeIn');
%     tempITDspks = su{ani}.ExtractTrialSpikes(squeeze(itdSpksBin./0.8), 'isc');
%     itdSpks{ani} = tempITDspks;
%     tempCorr = nan(1,size(tempTrlSpks,1));
%     for u = 1:size(tempTrlSpks,1)
% %         tempCorr(u) = corr(tempTrlSpks(u,:)', tempITDspks(u,:)');
%         tempCorr(u) = mean(tempTrlSpks(u,:)- tempITDspks(u,:));
%     end
%     itdTrlCorr{ani} = tempCorr;
% end
% meanTrlSpks = cell2mat(cellfun(@(a){mean(a,2)},trlSpks));
% meanItdSpks = cell2mat(cellfun(@(a){mean(a,2)},itdSpks));
% meanSpks = mean([meanTrlSpks,meanItdSpks],2);

%% Change in firing rate relative to max firing rate (dead end, ignore)
% winStepSize = 50;
% horizon = 10;
% trlType = 'isc';
% % trlType = 'iscSANSA';
% changeVals = cell(length(fileDirs),1);
% ndxs = cell(length(fileDirs),1);
% maxs = cell(length(fileDirs),1);
% for ani = 1:length(fileDirs)
%     [changeVals{ani},ndxs{ani},maxs{ani},windowTSvect] = su{ani}.SpikeDiffWindowedOnMaxRate([-1200 2000],'PokeIn',trlType,winStepSize,horizon);
% end
% grpChangeVals = cell2mat(changeVals);
% grpNdxs = cell2mat(ndxs);
% grpMaxs = cell2mat(maxs);
% % grpChangeVals(grpMaxs<10,:) = [];
% % grpNdxs(grpMaxs<10,:) = [];
% sortGroup = sortrows([grpNdxs,grpChangeVals]);
% figure;
% imagesc(windowTSvect,1:size(grpChangeVals,1), sortGroup(:,2:end));%, [0 0.5]);
% set(gca,'ydir', 'normal');
% colormap(cMap)
    