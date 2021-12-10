fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

binSize = 200;
dsRate = 150;
stWin = [-500 500];
ndWin = [500 500];
templateProportion = 0.5;
numPerms = 10;

%% Data Outputs
grpPosts = cell(1,1,length(fileDirs));
grpTimeDecode = cell(1,1,length(fileDirs));
grpOdrDecode = cell(1,1,length(fileDirs));
grpPosDecode = cell(1,1,length(fileDirs));

%%
for fl = 1:length(fileDirs)
    %% Create the MLB object
    pfcMLB = MLB_SM(fileDirs{fl});
    pfcMLB.binSize = binSize;
    pfcMLB.dsRate = dsRate;
    
    %% Extract trial start & end trial spike matrices
    [stTrlSpikes, stTrlTimeVect] = pfcMLB.PP_TrialMatrix_Spiking(stWin, 'PokeIn'); 
    [ndTrlSpikes, ndTrlTimeVect] = pfcMLB.PP_TrialMatrix_Spiking(ndWin, 'PokeOut'); 
    trlTimeVect = [stTrlTimeVect; ndTrlTimeVect+1.001];
    
    %% Create Likelihood & Observation vectors based on trial counts
    % Identify ISC trials and Odor
    iscTrialLog = [pfcMLB.trialInfo.ItemItemDistance]==1 & [pfcMLB.trialInfo.Performance]==1;
    odrVect = [pfcMLB.trialInfo(iscTrialLog).Odor];
    [counts,binEdges] = histcounts(odrVect);
    oipID = binEdges(2:end)-0.5;
    oipID(counts==0) = [];
    posID = arrayfun(@(a)mode([pfcMLB.trialInfo([pfcMLB.trialInfo.Odor]==a & iscTrialLog).Position]), oipID);
    counts(counts==0) = [];
    % Identify likelihood size based on proportion factor defined above and ISC trial #s
    likeliSize = floor(min(counts)*templateProportion);
    randPool = nan(likeliSize, length(counts), numPerms);
    % Create likelihood & observation data pools
    odrTrlObsvs = cell(length(counts),numPerms);
    odrTrlLikes = cell(length(counts),numPerms);
    % Compile the likelihoods & observations across repetitions
    for oip = 1:length(counts)
        tempPool = randi(counts(oip), [likeliSize,numPerms]);
        for r = 1:numPerms
            while length(unique(tempPool(:,r)))~=size(tempPool,1)
                tempPool(:,r) = [unique(tempPool(:,r)); randi(counts(oip), [size(tempPool,1)-length(unique(tempPool(:,r))),1])];
            end
        end
        while size(unique(tempPool', 'rows'),1)~=numPerms
            tempPool = [unique(tempPool', 'rows')', randi(counts(oip), [likeliSize,numPerms-size(unique(tempPool', 'rows'),1)])];
            for r = 1:numPerms
                while length(unique(tempPool(:,r)))~=size(tempPool,1)
                    tempPool(:,r) = [unique(tempPool(:,r)); randi(counts(oip), [size(tempPool,1)-length(unique(tempPool(:,r))),1])];
                end
            end
        end
        if size(unique(tempPool', 'rows'),1)~=numPerms
            error('Apparently I need to code for this eventuality too!?');
        end
%         tempPool = nchoosek(1:counts(oipID), likeliSize);
%         randSel = randi(size(tempPool,1), [1,numPerms]);
%         randPool(:,oipID,:) = tempPool(randSel,:)';
        randPool(:,oip,:) = tempPool;
        tempTrlSpikesST = stTrlSpikes(:,:,[pfcMLB.trialInfo.Odor]==oipID(oip) & iscTrialLog);
        tempTrlSpikesND = ndTrlSpikes(:,:,[pfcMLB.trialInfo.Odor]==oipID(oip) & iscTrialLog);
        for perm = 1:numPerms
            tempTempTrlSpikesST = tempTrlSpikesST;
            tempTempTrlSpikesND = tempTrlSpikesND;
            odrTrlLikes{oip,perm} = [tempTempTrlSpikesST(:,:,randPool(:,oip,perm));tempTempTrlSpikesND(:,:,randPool(:,oip,perm))];
            tempTempTrlSpikesST(:,:,randPool(:,oip,perm)) = [];
            tempTempTrlSpikesND(:,:,randPool(:,oip,perm)) = [];
            odrTrlObsvs{oip,perm} = [tempTempTrlSpikesST; tempTempTrlSpikesND];
        end
    end
    % Create posterior identity decoding vectors
    timeIDvect = repmat(trlTimeVect, [length(oipID),1]);
    oipIDvect = cell2mat(arrayfun(@(a)repmat(a,size(trlTimeVect)), oipID, 'uniformoutput', 0)');
    posIDvect = cell2mat(arrayfun(@(a)repmat(a,size(trlTimeVect)), posID, 'uniformoutput', 0)');
    %% Decode things!
    permPosts = nan(length(timeIDvect), length(timeIDvect), numPerms);
    permTimeDecode = cell(length(oipID), numPerms);
    permOdrDecode = cell(length(oipID), numPerms);
    permPosDecode = cell(length(oipID), numPerms);
    for perm = 1:numPerms
        tempLike = mean(cell2mat(odrTrlLikes(:,perm)),3,'omitnan');
        tempObsvs = odrTrlObsvs(:,perm);
        tempPosts = cell(size(tempObsvs));
        tempOdrPostRaw = cell(size(tempObsvs));
        for oip = 1:length(tempObsvs)
            tempOdrPost = pfcMLB.CalcStaticBayesPost(tempLike, tempObsvs{oip});
            permTimeDecode{oip,perm} = pfcMLB.DecodeBayesPost(tempOdrPost, timeIDvect) - repmat(trlTimeVect, [1,size(tempObsvs{oip},3)]);
            permOdrDecode{oip,perm} = pfcMLB.DecodeBayesPost(tempOdrPost, oipIDvect);
            permPosDecode{oip,perm} = pfcMLB.DecodeBayesPost(tempOdrPost, posIDvect);
            tempOdrPostRaw{oip} = reshape(permute(tempOdrPost, [2,1,3]), [size(tempOdrPost,2), numel(tempOdrPost)/size(tempOdrPost,2)])';
            tempPosts{oip} = mean(tempOdrPost,3,'omitnan');
        end
        permPosts(:,:,perm) = cell2mat(tempPosts);
    end
    grpPosts{fl} = permPosts;
    grpTimeDecode{fl} = permTimeDecode;
    grpOdrDecode{fl} = permOdrDecode;
    grpPosDecode{fl} = permPosDecode;
end
posID = unique(posID);
trlPrdID = [(stTrlTimeVect<0)*1 + (stTrlTimeVect>0)*2;(ndTrlTimeVect<0)*3 + (ndTrlTimeVect>0)*4];

%% Plot Group Posteriors
itemPosBounds = [find(diff(oipIDvect)~=0); length(oipIDvect)];
posTrlCenters = itemPosBounds-mode(diff(itemPosBounds))/2;
pokeInNdxs = find(timeIDvect==0);
pokeOutNdxs = find(timeIDvect==1.001);

figure;
imgPlot = subplot(7,4,1:16);
imagesc(mean(cell2mat(cellfun(@(a)mean(a,3), grpPosts, 'uniformoutput', 0)),3,'omitnan'),[0 0.015]);
set(gca, 'ydir', 'normal',...
    'xtick', posTrlCenters, 'xticklabel', arrayfun(@(a)sprintf('%02.0f',a),oipID,'UniformOutput',false),...
    'ytick', posTrlCenters, 'yticklabel', arrayfun(@(a)sprintf('%02.0f',a),oipID,'UniformOutput',false));
xlabel('Observed Time');
ylabel('Decoded Time');
title('Dual List Decoding Posteriors');
hold on;
plot([length(oipIDvect)/2+0.5, length(oipIDvect)/2+0.5], [0, length(oipIDvect)], '-k', 'linewidth', 2);
plot([0, length(oipIDvect)], [length(oipIDvect)/2+0.5, length(oipIDvect)/2+0.5], '-k', 'linewidth', 2);

load('batlowW.mat');
colormap(batlowW);
z = colormap;
% z = flipud(z);
% z(1,:) = [1 1 1];

aniTimeDecode = cell(length(oipID),length(fileDirs));
aniOdrDecode = cell(length(oipID),1,length(fileDirs));
aniOipTrlPrdDecode = cell(1,4,length(fileDirs));
aniPosDecode = cell(length(oipID),1,length(fileDirs));
for fl = 1:length(fileDirs)
    tempTimeDecode = grpTimeDecode{fl};
    tempOdrDecode = grpOdrDecode{fl};
    tempPosDecode = grpPosDecode{fl};
    tempTrlPrdDecode = repmat({nan(length(oipID))}, [1,4]);
    for oip = 1:length(oipID)
        % Temporal Decoding
        tempPosTimeDecode = cell2mat(tempTimeDecode(oip,:));
        aniTimeDecode{oip,fl} = mean(tempPosTimeDecode,2);
        % Odor Decoding
        tempTempOdrDecode = cell2mat(tempOdrDecode(oip,:));
        odrDecodes = nan(size(tempTempOdrDecode,1),length(oipID));
        for o = 1:length(oipID)
            odrDecodes(:,o) = mean(tempTempOdrDecode==oipID(o),2);
            for trlPrd = 1:4
                tempTrlPrdDecode{trlPrd}(oip,o) = sum(sum(tempTempOdrDecode==oipID(o) & trlPrdID==trlPrd))/(sum(trlPrdID==trlPrd)*size(tempTempOdrDecode,2));
            end
        end
        aniOdrDecode{oip,1,fl} = odrDecodes;
        % Position Decoding
        tempTempPosDecode = cell2mat(tempPosDecode(oip,:));
        posDecodes = nan(size(tempTempPosDecode,1),length(posID));
        for p = 1:length(posID)
            posDecodes(:,p) = mean(tempTempPosDecode==posID(p),2);
        end
        aniPosDecode{oip,1,fl} = posDecodes;
    end
    aniOipTrlPrdDecode(:,:,fl) = tempTrlPrdDecode;
end
lnPlot(1) = subplot(7,4,17:20);
timeDecode = cell2mat(aniTimeDecode);
meanTmDec = mean(timeDecode,2);
semTmDec = pfcMLB.SEMcalc(timeDecode')';
plot(1:length(timeIDvect), meanTmDec, '-k');
hold on;
patch('XData', [1:length(timeIDvect), length(timeIDvect):-1:1],...
    'YData', [(meanTmDec+semTmDec)', flipud(meanTmDec-semTmDec)'], 'edgecolor', 'k',...
    'facecolor', 'k', 'facealpha', 0.5);
for oip = 1:length(pokeInNdxs)
    plot(repmat(pokeInNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(pokeOutNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(itemPosBounds(oip), [1,2]), get(gca, 'ylim'), '--k');
end
box off;
set(gca, 'color', 'none', 'XAxisLocation', 'origin', 'xtick', posTrlCenters, 'xticklabel', []);
title('Temporal Discrepancy');
axis tight

lnPlot(2) = subplot(7,4,21:24);
odrDecode = cell2mat(aniOdrDecode);
for oip = 1:size(odrDecode,2)
    tempOIP = reshape(odrDecode(:,oip,:), [size(odrDecode,1), size(odrDecode,3)]);
    tempMeanOdrDec = mean(tempOIP,2);
    tempSemOdrDec = pfcMLB.SEMcalc(tempOIP')';
    if oipID(oip)<10
        plot(1:length(oipIDvect),tempMeanOdrDec, 'color', pfcMLB.PositionColors(oipID(oip),:), 'linestyle', '-');
        hold on;
        patch('XData', [1:length(oipIDvect), length(oipIDvect):-1:1],...
            'YData', [(tempMeanOdrDec+tempSemOdrDec)', flipud(tempMeanOdrDec-tempSemOdrDec)'], 'edgecolor', pfcMLB.PositionColors(oipID(oip),:),...
            'facecolor', pfcMLB.PositionColors(oipID(oip),:), 'facealpha', 0.5);
    else
        plot(1:length(oipIDvect),tempMeanOdrDec, 'color', pfcMLB.PositionColors(oipID(oip)-10,:), 'linestyle', '--');
        hold on;
        patch('XData', [1:length(oipIDvect), length(oipIDvect):-1:1],...
            'YData', [(tempMeanOdrDec+tempSemOdrDec)', flipud(tempMeanOdrDec-tempSemOdrDec)'], 'edgecolor', pfcMLB.PositionColors(oipID(oip)-10,:),...
            'linestyle', '--', 'facecolor', pfcMLB.PositionColors(oipID(oip)-10,:), 'facealpha', 0.5);
    end
end
for oip = 1:length(pokeInNdxs)
    plot(repmat(pokeInNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(pokeOutNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(itemPosBounds(oip), [1,2]), get(gca, 'ylim'), '--k');
end
set(gca, 'xtick', posTrlCenters, 'xticklabel', []);
title('Odor-In-Position Decoding');

lnPlot(3) = subplot(7,4,25:28);
posDecode = cell2mat(aniPosDecode);
for pos = 1:size(posDecode,2)
    tempPOS = reshape(posDecode(:,pos,:), [size(posDecode,1), size(posDecode,3)]);
    tempMeanPosDec = mean(tempPOS,2);
    tempSemPosDec = pfcMLB.SEMcalc(tempPOS')';
    plot(1:length(posIDvect), tempMeanPosDec, 'color', pfcMLB.PositionColors(posID(pos),:), 'linestyle', '-');
    hold on;
    patch('XData', [1:length(posIDvect), length(posIDvect):-1:1],...
        'YData', [(tempMeanPosDec+tempSemPosDec)', flipud(tempMeanPosDec-tempSemPosDec)'], 'edgecolor', pfcMLB.PositionColors(posID(pos),:),...
        'facecolor', pfcMLB.PositionColors(oipID(pos),:), 'facealpha', 0.5);
end
for oip = 1:length(pokeInNdxs)
    plot(repmat(pokeInNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(pokeOutNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
    plot(repmat(itemPosBounds(oip), [1,2]), get(gca, 'ylim'), '--k');
end
set(gca, 'xtick', posTrlCenters, 'xticklabel',  arrayfun(@(a)sprintf('%02.0f',a),oipID,'UniformOutput',false));
title('Positional Decoding')
    

axis(lnPlot, 'tight');        
linkaxes([imgPlot, lnPlot], 'x');
            
annotation(gcf,'textbox', [0.1 0.95 0.8 0.05],...
    'String', sprintf('Dual List Group:\n binSize = %.0fms, dsRate = %.0fms, SartWindow = (%.0fms:%.0fms to pokeIn), EndWindow = (%.0fms:%.0fms to pokeOut), Likelihood Porportion = %.02f, NumPerms = %.0f',...
    binSize, dsRate, stWin(1), stWin(2), ndWin(1), ndWin(2), templateProportion, numPerms),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%%
if sum(oipID>10)>=1
    %%
    figure;
    for sp = 1:4
        subplot(2,2,sp)
        imagesc(median(cell2mat(aniOipTrlPrdDecode(1,sp,:)),3,'omitnan'), [0 0.5]);
        set(gca, 'xtick', 1:length(oipID), 'xticklabel', arrayfun(@(a)sprintf('%02.0f',a),oipID,'UniformOutput',false),...
            'ytick', 1:length(oipID), 'yticklabel', arrayfun(@(a)sprintf('%02.0f',a),oipID,'UniformOutput',false));
        title(sp);
    end
    %%
    odrDiffs = nan(size(odrDecode,1), length(fileDirs), length(posID));
    seqDiffs = nan(size(odrDecode,1)/length(oipID), length(fileDirs), 2);
    for ani = 1:length(fileDirs)
        tempSeqDiff = nan(size(odrDecode,1)/length(oipID), max(posID), 2);
        for pos = 1:max(posID)
            odrDiffs(:,ani,pos) = odrDecode(:,pos,ani)-odrDecode(:,pos+max(posID),ani);
            odrDiffs(posIDvect~=pos,ani,pos) = 0;
            tempSeqDiff(:,pos,1) = odrDiffs(oipIDvect==pos,ani,pos);
            tempSeqDiff(:,pos,2) = odrDiffs(oipIDvect==(pos+10),ani,pos);
        end
        seqDiffs(:,ani,1) = mean(tempSeqDiff(:,:,1),2);
        seqDiffs(:,ani,2) = mean(tempSeqDiff(:,:,2),2);
    end
    figure;
    subplot(1,3,1:2)
    for pos = 1:max(posID)
        tempMeanOdrDiff = mean(odrDiffs(:,:,pos),2);
        tempSemOdrDiff = pfcMLB.SEMcalc(odrDiffs(:,:,pos)')';
        plot(tempMeanOdrDiff, 'color', pfcMLB.PositionColors(pos,:));
        hold on;
        patch('XData', [1:size(odrDiffs,1), size(odrDiffs,1):-1:1],...
        'YData', [(tempMeanOdrDiff+tempSemOdrDiff)', flipud(tempMeanOdrDiff-tempSemOdrDiff)'], 'edgecolor', pfcMLB.PositionColors(pos,:),...
        'facecolor', pfcMLB.PositionColors(pos,:), 'facealpha', 0.5);
    end
    for oip = 1:length(pokeInNdxs)
        plot(repmat(pokeInNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
        plot(repmat(pokeOutNdxs(oip), [1,2]), get(gca, 'ylim'), ':k');
        plot(repmat(itemPosBounds(oip), [1,2]), get(gca, 'ylim'), '--k');
    end
    box off;
    set(gca, 'color', 'none', 'XAxisLocation', 'origin', 'xtick', posTrlCenters, 'xticklabel', []);
    subplot(1,3,3)
    seq1SeqDiffMean = mean(seqDiffs(:,:,1),2);
    seq1SeqDiffSEM = pfcMLB.SEMcalc(seqDiffs(:,:,1)')';
    plot(trlTimeVect,seq1SeqDiffMean, '-k');
    hold on;
    patch('XData', [trlTimeVect;flipud(trlTimeVect)],...
        'YData', [seq1SeqDiffMean+seq1SeqDiffSEM ; flipud(seq1SeqDiffMean-seq1SeqDiffSEM)], 'edgecolor', 'k',...
        'facecolor', 'k', 'facealpha', 0.5);
    seq2SeqDiffMean = mean(seqDiffs(:,:,2),2);
    seq2SeqDiffSEM = pfcMLB.SEMcalc(seqDiffs(:,:,2)')';
    plot(trlTimeVect,seq2SeqDiffMean, '--k');
    hold on;
    patch('XData', [trlTimeVect; flipud(trlTimeVect)],...
        'YData', [seq2SeqDiffMean+seq2SeqDiffSEM; flipud(seq2SeqDiffMean - seq2SeqDiffSEM)], 'edgecolor', 'k', 'linestyle', '--',...
        'facecolor', 'k', 'facealpha', 0.5);
    set(gca, 'color', 'none', 'XAxisLocation', 'origin');
    box off;
    axis tight
end

    