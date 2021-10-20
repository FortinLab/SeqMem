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

binSize = 100;
dsRate = 50;

odrPosts = cell(1,length(fileDirs));
odrDecodes = cell(1,length(fileDirs));
trlInfo = odrPosts;
fiscs = odrPosts;
meanPosts = cell(4,4,length(fileDirs));
meanDecodes = meanPosts;
lfpPhase = odrPosts;
lfpPower = odrPosts;
for fl = 1:length(fileDirs)
    tic;
    [odrPosts{fl}, odrDecodes{fl}, trlInfo{fl}, fiscs{fl}, meanPosts(:,:,fl), meanDecodes(:,:,fl), lfpPhase{fl}, lfpPower{fl}, trlTimeVect]...
        = TemporalInvariance_MLB(fileDirs{fl}, binSize, dsRate);
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
colormap(hot);
z = colormap;
z = flipud(z);
% z(1,:) = [1 1 1];
for p = 1:4
    for pp = 1:4
        subplot(4,4,sub2ind([4,4],p,pp));
        tempMean = mean(cell2mat(meanDecodes(p,pp,:)),3,'omitnan');
        tempMean(tempMean<0.3) = nan;
        imagesc(trlTimeVect,trlTimeVect, tempMean, [0 0.5]);
        hold on;
        plot(trlTimeVect,trlTimeVect, '--k');
        colormap(z);
        title(sprintf('Decode %i During %i', pp, p));
        xlabel('Template Time');
        ylabel('Decode Time');
        set(gca, 'ydir', 'normal');
    end
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%%
% histogram(cell2mat(cellfun(@(a)a(:), cellfun(@(a)mean(a), lfpPower, 'uniformoutput', 0), 'uniformoutput', 0)'))
%%
isc = cell(1,4);
iscPwr = cell(1,4);
for ani = 1:length(fileDirs)
    tempISlog = [trlInfo{ani}.Odor]==[trlInfo{ani}.Position];
    tempPerfLog = [trlInfo{ani}.Performance];
    tempPwr = mean(lfpPower{ani});
    tempPwr = tempPwr(:);
    for op = 1:4
        isc{op} = [isc{op}; odrPosts{ani}(tempISlog & tempPerfLog & [trlInfo{ani}.Odor]==op)];
        iscPwr{op} = [iscPwr{op}; tempPwr(tempISlog & tempPerfLog & [trlInfo{ani}.Odor]==op)];
    end
end

figure;
for p = 1:4
    tempTrls = isc{p};
    tempPwr = iscPwr{p};
%     tempTrls(tempPwr<0) = [];
    for pp = 1:4
        subplot(4,4,sub2ind([4,4],p,pp));
        trlMean = mean(cell2mat(reshape(cellfun(@(a){a(:,:,pp)}, tempTrls), [1,1,size(tempTrls,1)])),3,'omitnan');
        trlMean(trlMean<0.3)=nan;
        imagesc(trlTimeVect, trlTimeVect, trlMean, [0 0.5]);
        hold on;
        plot(trlTimeVect,trlTimeVect, '--k');
        colormap(z);
        title(sprintf('Decode %i During %i', pp, p));
        xlabel('Template Time');
        ylabel('Decode Time');
        set(gca, 'ydir', 'normal');
    end
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
%%
lagVect = -3:3;
trlLag = cell(size(lagVect));
for p = 1:4
    tempTrls = isc{p};
    for pp = 1:4
        curLag = lagVect==pp-p;
        trlLag{curLag} = [trlLag{curLag}; cellfun(@(a){a(:,:,pp)}, tempTrls)];
    end
end
trlLagMean = cell(size(lagVect));
figure 
for lag = 1:length(lagVect)
    subplot(length(lagVect),1, lag);
    curLagMean = mean(cell2mat(reshape(trlLag{lag}, [1,1,length(trlLag{lag})])),3,'omitnan');
    curLagMean(curLagMean<0.3)=nan;
    imagesc(trlTimeVect, trlTimeVect, curLagMean, [0 0.5]);
    hold on;
    plot(trlTimeVect,trlTimeVect, '--k');
    title(sprintf('Decode Lag%i', lagVect(lag)));
    xlabel('Template Time');
    ylabel('Decode Time');
    set(gca, 'ydir', 'normal');
    colormap(z);
end
annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i', binSize, dsRate),...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
%%
objPosDecode = cell(1,5);
for p = 1:4
    tempTrls = isc{p};
    for pp = 1:4
        if p==pp
            objPosDecode{5} = [objPosDecode{5}; cellfun(@(a){a(:,:,pp)}, tempTrls)];
        else
            objPosDecode{pp} = [objPosDecode{pp}; cellfun(@(a){a(:,:,pp)},tempTrls)];
        end
    end
end
figure;
subplots = [1, 3, 7, 9, 5];
ttls = [{'pos1'}, {'pos2'}, {'pos3'}, {'pos4'}, {'Matched'}];
for sp = 1:length(subplots)
    subplot(3,3,subplots(sp));
    curMean = mean(cell2mat(reshape(objPosDecode{sp}, [1,1,length(objPosDecode{sp})])),3,'omitnan');
    curMean(curMean<0.3)=nan;
    imagesc(trlTimeVect, trlTimeVect, curMean, [0 0.5]);
    hold on;
    plot(trlTimeVect,trlTimeVect, '--k');
    title(ttls(sp));
    xlabel('Template Time');
    ylabel('Decode Time');
    set(gca, 'ydir', 'normal');
    colormap(z);
end