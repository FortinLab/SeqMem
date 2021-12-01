% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];

binSize = 200;
dsRate = 50;
trlWindow = [-800 1200];

PositionColors = [44/255, 168/255, 224/255;...
    154/255, 133/255, 122/255;...
    9/255, 161/255, 74/255;...
    128/255, 66/255, 151/255;...
    241/255, 103/255, 36/255];

odrPosts = cell(1,length(fileDirs));
odrDecodes = cell(1,length(fileDirs));
trlInfo = odrPosts;
fiscs = odrPosts;
meanPosts = cell(5,5,length(fileDirs));
% meanPosts = cell(4,4,length(fileDirs));
meanDecodes = meanPosts;
lfpPhase = odrPosts;
lfpPower = odrPosts;
for fl = 1:length(fileDirs)
    tic;
    [odrPosts{fl}, odrDecodes{fl}, trlInfo{fl}, fiscs{fl}, meanPosts(:,:,fl), meanDecodes(:,:,fl), lfpPhase{fl}, lfpPower{fl}, trlTimeVect]...
        = TemporalInvariance_MLB(fileDirs{fl}, binSize, dsRate, trlWindow);
    toc
    
end
%
close all;
%%
% figure;
% for p = 1:4
%     for pp = 1:4
%         subplot(4,4,sub2ind([4,4],pp,p));        
%         imagesc(trlTimeVect,trlTimeVect, median(cell2mat(meanPosts(p,pp,:)),3,'omitnan'), [0 0.5]);
%         title(sprintf('Post %i During %i', pp, p));
%         xlabel('Template Time');
%         ylabel('Decode Time');
%         set(gca, 'ydir', 'normal');
%         drawnow
%     end
% end
figure;
load('broc.mat');
colormap(broc);
z = colormap;
% z = flipud(z);
% z(1,:) = [1 1 1];
for p = 1:size(fiscs{1},1)
    for pp = 1:size(fiscs{1},1)
        subplot(size(fiscs{1},1),size(fiscs{1},1),sub2ind([size(fiscs{1},1),size(fiscs{1},1)],pp,p));
        threshMean = mean(cell2mat(meanDecodes(p,pp,:)),3,'omitnan');
%         tempMean(tempMean<0.3) = nan;
        imagesc(trlTimeVect,trlTimeVect, threshMean, [0 0.5]);
        hold on;
        plot(trlTimeVect,trlTimeVect, '--k');
        colormap(z);
        title(sprintf('Decode %i During %i', pp, p));
        xlabel('Template Time');
        ylabel('Trial Time');
        set(gca, 'ydir', 'normal');
    end
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Mean Animal; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Look at difference between Hits & False Alarms

iscDiffs = cell(1,size(fiscs{1},1));
for p = 1:size(fiscs{1},1)
    for a = 1:length(fileDirs)
        tempHit = meanDecodes{p,p,a};
        tempFA = median(cell2mat(reshape(meanDecodes(1:size(fiscs{1},1)~=p,p,a), [1,1,size(fiscs{1},1)-1])),3, 'omitnan');
        iscDiffs{1,p,a} = tempHit-tempFA;
    end
end

figure;
for p = 1:size(fiscs{1},1)
    subplot(1,size(fiscs{1},1),p);
    tempMean = median(cell2mat(iscDiffs(1,p,:)),3);
    imagesc(trlTimeVect,trlTimeVect, tempMean, [0 0.5]);
    hold on;
    plot(trlTimeVect,trlTimeVect, '--k');
    colormap(z);
    title(sprintf('Decode %i During %i', p, p));
    xlabel('Template Time');
    ylabel('Trial Time');
    set(gca, 'ydir', 'normal');
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Hit Rate - False Alarm; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Calculate timecourse of LFP power & Compare with HR-FA Decoding
iscPwr = cell(1,size(fiscs{1},1));
for ani = 1:length(fileDirs)
    tempISlog = [trlInfo{ani}.Odor]==[trlInfo{ani}.Position];
    tempPerfLog = [trlInfo{ani}.Performance];
    tempPwr = lfpPower{ani};
    nonAfiscTrls = fiscs{ani}(2:end,:);
    iscPwr{ani} = mean(tempPwr(:,:,nonAfiscTrls(:)),3,'omitnan');
end
iscPwr = cell2mat(iscPwr);

nonA = mean(cell2mat(reshape(iscDiffs(:,2:end,:), [1,1,numel(iscDiffs(:,2:end,:))])),3);
figure; 
sp1 = subplot(5,4,1:16);
imagesc(trlTimeVect, trlTimeVect, nonA, [0 0.5]);
hold on;
plot(trlTimeVect, trlTimeVect, '--k');
colormap(z);
title('Matched Decoding (non A1 Trials)');
ylabel('Trial Time');
set(gca, 'ydir', 'normal');
sp2 = subplot(5,4,17:20);
plot(trlTimeVect, mean(iscPwr,2), '-k');
patch('XData', [trlTimeVect; flipud(trlTimeVect)], 'YData', [mean(iscPwr,2)-SEMcalc(iscPwr')'; flipud(mean(iscPwr,2)+SEMcalc(iscPwr')')],...
    'facealpha', 0.5);
linkaxes([sp1 sp2], 'x');
title('Theta Power (animal mean +/- SEM)');
ylabel('Power');
xlabel('Template Time');
axis tight

axes(gcf, 'position',[0.7 0.3 0.2 0.2]);
corrScatPlot(trapz(nonA)', mean(iscPwr,2));
% set(axes, 'fontcolor', 'w')
        
%%
% histogram(cell2mat(cellfun(@(a)a(:), cellfun(@(a)mean(a), lfpPower, 'uniformoutput', 0), 'uniformoutput', 0)'))
%%
isc = cell(1,size(fiscs{1},1));
iscPwr = cell(1,size(fiscs{1},1));
for ani = 1:length(fileDirs)
    tempISlog = [trlInfo{ani}.Odor]==[trlInfo{ani}.Position];
    tempPerfLog = [trlInfo{ani}.Performance];
    tempPwr = trapz(lfpPower{ani});
    tempPwr = tempPwr(:);
    for op = 1:size(fiscs{1},1)
%         tempTempISlog = tempISlog - ismember([trlInfo{ani}.TrialNum], fiscs{ani}(op,:));
%         tempTempISlog = ismember([trlInfo{ani}.TrialNum], fiscs{ani}(op,:));
        tempTempISlog = tempISlog;
        isc{op} = [isc{op}; odrPosts{ani}(tempTempISlog & tempPerfLog & [trlInfo{ani}.Odor]==op)];
        iscPwr{op} = [iscPwr{op}; tempPwr(tempTempISlog & tempPerfLog & [trlInfo{ani}.Odor]==op)];
    end
end

pwr=cell2mat(iscPwr(:));

figure;
for p = 1:size(fiscs{1},1)
    tempTrls = isc{p};
    tempPwr = iscPwr{p};
    tempTrls(tempPwr>median(pwr)) = [];
    for pp = 1:size(fiscs{1},1)
        subplot(size(fiscs{1},1),size(fiscs{1},1),sub2ind([size(fiscs{1},1),size(fiscs{1},1)],pp,p));
        trlMean = mean(cell2mat(reshape(cellfun(@(a){a(:,:,pp)}, tempTrls), [1,1,size(tempTrls,1)])),3,'omitnan');
%         trlMean(trlMean<0.3)=nan;
        threshMean = trlMean;
        threshMean(threshMean<0.35) = 0;
        threshMean(threshMean>0) = 1;
        imagesc(trlTimeVect, trlTimeVect, trlMean, [0 0.5]);
        hold on;
        plot(trlTimeVect,trlTimeVect, '--k');
        bounds = bwboundaries(threshMean);
%         for b = 1:length(bounds)
%             plot(trlTimeVect(bounds{b}(:,2)), trlTimeVect(bounds{b}(:,1)), 'k', 'linewidth', 5);
%         end
        colormap(z);
        title(sprintf('Decode %i During %i', pp, p));
        xlabel('Template Time');
        ylabel('Trial Time');
        set(gca, 'ydir', 'normal');
    end
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Mean Trials; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
%%
lagVect = (size(fiscs{1},1)-2)*-1:(size(fiscs{1},1)-2);
trlLag = cell(size(lagVect));
for p = 2:size(fiscs{1},1)
    tempTrls = isc{p};
    for pp = 2:size(fiscs{1},1)
        curLag = lagVect==pp-p;
        trlLag{curLag} = [trlLag{curLag}; cellfun(@(a){a(:,:,pp)}, tempTrls)];
    end
end
trlLagMean = cell(size(lagVect));
figure 
for lag = 1:length(lagVect)
    subplot(length(lagVect),1, lag);
    curLagMean = mean(cell2mat(reshape(trlLag{lag}, [1,1,length(trlLag{lag})])),3,'omitnan');
%     curLagMean(curLagMean<0.3)=nan;
    imagesc(trlTimeVect, trlTimeVect, curLagMean, [0 0.5]);
    hold on;
    plot(trlTimeVect,trlTimeVect, '--k');
    title(sprintf('Decode Lag%i', lagVect(lag)));
    xlabel('Template Time');
    ylabel('Trial Time');
    set(gca, 'ydir', 'normal');
    colormap(z);
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Decoding by Lag Collapsed; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%%

figure 
for lag = 1:length(lagVect)
    subplot(length(lagVect),1, lag);
    curLagMean = mean(cell2mat(reshape(trlLag{lag}, [1,1,length(trlLag{lag})])),3,'omitnan');
    tempLagMean = curLagMean;
    curLagMean(curLagMean<0.35) = 0;
    curLagMean(curLagMean>0) = 1;
    imagesc(trlTimeVect, trlTimeVect, tempLagMean, [0 0.5]);
    hold on;
    
    bounds = bwboundaries(curLagMean);
    for b = 1:length(bounds)
        plot(trlTimeVect(bounds{b}(:,2)), trlTimeVect(bounds{b}(:,1)), 'k', 'linewidth', 5);
    end
    set(gca, 'ydir', 'normal');
    colormap(z);
end
    
%%
objPosDecode = cell(1,size(fiscs{1},1)+1);
for p = 1:size(fiscs{1},1)
    tempTrls = isc{p};
    for pp = 1:size(fiscs{1},1)
        if p==pp
            if p ~= 1 && p~=4
                objPosDecode{end} = [objPosDecode{end}; cellfun(@(a){a(:,:,pp)}, tempTrls)];
            end
        else
            objPosDecode{pp} = [objPosDecode{pp}; cellfun(@(a){a(:,:,pp)},tempTrls)];
        end
    end
end
figure;
if size(fiscs{1},1) == 4
    subplots = [1, 3, 5, 7, 2];
    ttls = [{'pos1'}, {'pos2'}, {'pos3'}, {'pos4'}, {'Matched'}];
elseif size(fiscs{1},1) == 5
    subplots = [1, 3, 5, 7, 9, 2];
    ttls = [{'pos1'}, {'pos2'}, {'pos3'}, {'pos4'}, {'pos5'}, {'Matched'}];
end

for sp = 1:length(subplots)
    subplot(size(fiscs{1},1),2,subplots(sp));
    curMean = mean(cell2mat(reshape(objPosDecode{sp}, [1,1,length(objPosDecode{sp})])),3,'omitnan');
%     curMean(curMean<0.3)=nan;
    imagesc(trlTimeVect, trlTimeVect, curMean, [0 0.5]);
    hold on;
    plot(trlTimeVect,trlTimeVect, '--k');
    title(ttls(sp));
    xlabel('Template Time');
    ylabel('Trial Time');
    set(gca, 'ydir', 'normal');
    colormap(z);
end

annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Non-Match Decoding; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
%%
timeWindow = [-0.6 -0.4];
figure;
for sp = 1:length(subplots)
    curMean = mean(cell2mat(reshape(objPosDecode{sp}, [1,1,length(objPosDecode{sp})])),3,'omitnan');
    windowMean = mean(curMean(:,trlTimeVect>timeWindow(1) & trlTimeVect<timeWindow(2)),2);
%     windowMean = windowMean./max(windowMean);
    if sum(sp == 1:4)==1
        plot(trlTimeVect,windowMean, 'color', PositionColors(sp,:), 'linewidth', 2);
    else
        plot(trlTimeVect,windowMean, 'color', 'k', 'linewidth', 2);
    end
    hold on;
end
legend([{'1'}, {'2'}, {'3'}, {'4'}, {'Match'}]);
title(sprintf('Template between %ims and %ims', timeWindow(1)*1000, timeWindow(2)*1000));
xlabel('Trial Time');
ylabel('Probability (% Trials decoded)');


annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Time-Slice Match/Non-Match Decoding; BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
%%
for t = 1:length(trlTimeVect)
    if trlTimeVect(t)-binSize/2000 <= trlTimeVect(1)
        rlyWindow = trlTimeVect(1);
    else
        rlyWindow = trlTimeVect(t)-binSize/2000;
    end
    if trlTimeVect(t)+binSize/2000 >= trlTimeVect(end)
        latWindow = trlTimeVect(end);
    else
        latWindow = trlTimeVect(t)+binSize/2000;
    end
    timeWindow = [rlyWindow latWindow];
%     timeWindow = [0.49 0.21];
    figure;
    pPlt = nan(1,size(fiscs{1},1));
    sp = nan(1,size(fiscs{1},1));
    for p = 1:size(fiscs{1},1)
        sp(p) = subplot(size(fiscs{1},1),1,p);
        patch('XData', [timeWindow, fliplr(timeWindow)], 'YData', [0 0 1 1], 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.2);
        hold on;
        tempTrials = isc{p};
        for pp = 1:size(fiscs{1},1)
            tempDecode = cell2mat(cellfun(@(a){mean(a(:,trlTimeVect>timeWindow(1) & trlTimeVect<timeWindow(2),pp),2)}, tempTrials)');
            meanDecode = mean(tempDecode,2);
            semDecode = SEMcalc(tempDecode,0,2);
            patch('XData', [trlTimeVect', fliplr(trlTimeVect')], 'YData', [(meanDecode-semDecode)', fliplr((meanDecode+semDecode)')],...
                'EdgeColor', PositionColors(pp,:), 'FaceColor', PositionColors(pp,:), 'FaceAlpha', 0.5);
            pPlt(pp) = plot(trlTimeVect, meanDecode, 'color', PositionColors(pp,:), 'linewidth', 2);
        end
        title(sprintf('Position %i', p));
    end
    if size(fiscs{1},1) == 4
        legend(pPlt,[{'1'}, {'2'}, {'3'}, {'4'}], 'location', 'north', 'orientation', 'horizontal');
    elseif size(fiscs{1},1) == 5
        legend(pPlt,[{'1'}, {'2'}, {'3'}, {'4'}, {'5'}], 'location', 'north', 'orientation', 'horizontal');
    end
    linkaxes(sp, 'xy');
    axis tight
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('Time-Slice Decoding; Centered @%.00fms; BinSize = %i, DSrate = %i',trlTimeVect(t)*1000, binSize, dsRate),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    orient(gcf, 'tall')
%     print(gcf, '-dpdf', sprintf('plot%i.pdf', t));
    close(gcf);
end
    

