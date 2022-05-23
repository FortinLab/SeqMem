% PFC_SingleUnit_Modulation
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

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

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
%
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

%% Colormap setup
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);
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
    su{ani}.savePlots = true;
    su{ani}.PlotUnitSummary;
    su{ani}.savePlots = false;
    close all
end

%% Evaluate Position Info coding w/in ensemble
posInfo = cell(1,length(fileDirs));
posInfoSANSA = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    su{ani}.dsRate = 10;
    su{ani}.binSize = 200;
    [posInfo{ani}, tsVect] = su{ani}.QuantPosInfo([-800 2000], 'PokeIn', 'pos', 'corr');
    posInfoSANSA{ani} = su{ani}.QuantPosInfo([-800 2000], 'PokeIn', 'pos', 'corrSANSA');
end
%
ensmblPosInfo = cell2mat(posInfo);
ensmblPosInfoSANSA = cell2mat(posInfoSANSA);
normPI = nan(size(ensmblPosInfo));
normPIsansa = nan(size(ensmblPosInfo));
for u = 1:size(ensmblPosInfo,2)
    normPI(:,u) = ensmblPosInfo(:,u)./max(ensmblPosInfo(:,u));
    normPIsansa(:,u) = ensmblPosInfoSANSA(:,u)./max(ensmblPosInfoSANSA(:,u));
end
normPI = normPI';
normPIsansa = normPIsansa';
loc = nan(size(normPI,1),1);
locSANSA = nan(size(normPI,1),1);
for u = 1:size(normPI,1)
    loc(u) = find(normPI(u,:)==1,1,'first');
    locSANSA(u) = find(normPIsansa(u,:)==1,1,'first');
end
sortNormPI = sortrows([loc, normPI]);
sortNormPI = sortNormPI(:,2:end);
altSortNormPI = sortrows([locSANSA, normPI]);
altSortNormPI = altSortNormPI(:,2:end);
sortNormPIsansa = sortrows([locSANSA, normPIsansa]);
sortNormPIsansa = sortNormPIsansa(:,2:end);
altSortNormPIsansa = sortrows([loc, normPIsansa]);
altSortNormPIsansa = altSortNormPIsansa(:,2:end);

figure; 
subplot(2,2,1);
imagesc(tsVect, 1:size(sortNormPI,1), sortNormPI);
title('All Trials (sort by all)');
subplot(2,2,2);
imagesc(tsVect, 1:size(altSortNormPI,1), altSortNormPI);
title('All Trials (sort by SANSA)');
subplot(2,2,3);
imagesc(tsVect, 1:size(altSortNormPIsansa,1), altSortNormPIsansa);
title('SANSA (sort by all)');
subplot(2,2,4);
imagesc(tsVect, 1:size(sortNormPIsansa,1), sortNormPI);
title('SANSA (sort by SANSA)');
colormap(cMap);

%% Pos Info vs Spike Rate
aniPosSpkCorr = cell(length(fileDirs), 1);
aniUniPrd = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    iscPokeDur = [su{ani}.trialInfo([su{ani}.trialInfo.TranspositionDistance]==0 & [su{ani}.trialInfo.Performance]==1).PokeDuration];
    windowEnd = floor((mean(iscPokeDur)+0.6)*100)/100;
%     windowEnd = floor((mean(iscPokeDur))*100)/100;
    su{ani}.dsRate = 10;
    su{ani}.binSize = 200;
    [aniPosSpkCorr{ani}, posInfo, varSpks, tsVect] = su{ani}.QuantPosSpkCorr([-1200 windowEnd*1000],'PokeIn','pos','corrSANSA');
    posSpkNdxs = nan(2,size(varSpks,2));
    for u = 1:size(varSpks,2)
        posSpkNdxs(1,u) = tsVect(find(posInfo(:,u)==max(posInfo(:,u)),1,'first'));
        posSpkNdxs(2,u) = tsVect(find(varSpks(:,u,end)==max(varSpks(:,u,end)),1,'first'));
    end
    aniUniPrd{ani} = posSpkNdxs;
end

corrVals = cell2mat(aniPosSpkCorr);
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