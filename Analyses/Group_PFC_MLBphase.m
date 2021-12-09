% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = {'D:\WorkBigDataFiles\CA3 Data\HC01_071217\Cut3\'};

fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
    {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];

binSize = 100;
dsRate = 50;
stWin = [0 500];
dnWin = [-500 0];

phaseBins = -pi:pi/4:pi;
powerBins = -2:0.5:6;
odorBins = 0.5:4.5;
timeBins = -0.5:0.05:0.5;

powerThresh = -4;

% nbp = 1;
% nbp = 0;
nbp = nan;

%%
%%
iscDecodesOdor = cell(4,length(fileDirs));
iscDecodesTime = cell(4,length(fileDirs));
thetaPhase = cell(4,length(fileDirs));
thetaPower = cell(4,length(fileDirs));
betaPhase = cell(4,length(fileDirs));
betaPower = cell(4,length(fileDirs));
for fl = 1:length(fileDirs)
    pfcMLB = PFC_TrialEvent_MLB_SM(fileDirs{fl});
    pfcMLB.binSize = binSize;
    pfcMLB.dsRate = dsRate;
    pfcMLB.beginTrialWindow = stWin;
    pfcMLB.endTrialWindow = dnWin;
%     pfcMLB.lfpRefTet = 17;
    pfcMLB.RunAnalysis;
    o = figure('Position', [88 350 560 420]);
    t = figure('Position', [1248 350 560 420]);
    for op = 1:4
        iscLog = ([pfcMLB.trialInfo.Odor] == [pfcMLB.trialInfo.Position]) &...
            ([pfcMLB.trialInfo.Performance] == 1) &...
            ([pfcMLB.trialInfo.Odor]==op);
        tempISCposts = pfcMLB.trlPosts(:,:,iscLog);
        [tempISCdecodeOdor, ~] = pfcMLB.DecodeBayesPost(tempISCposts, pfcMLB.fisSeqSpikeOdorLog);
        iscDecodesOdor{op,fl} = tempISCdecodeOdor;
        tempISCdecodeTime = pfcMLB.DecodeBayesPost(tempISCposts, pfcMLB.fisSeqSpikeTimeLog) - pfcMLB.trialTimeLog;
        iscDecodesTime{op,fl} = tempISCdecodeTime;
        tempThetaPhase = pfcMLB.trialThetaPhaseMtx(:,iscLog);
        thetaPhase{op,fl} = tempThetaPhase;
        tempThetaPower = pfcMLB.trialThetaPowerMtx(:,iscLog);
        thetaPower{op,fl} = tempThetaPower;
        tempBetaPhase = pfcMLB.trialBetaPhaseMtx(:,iscLog);
        betaPhase{op,fl} = tempBetaPhase;
        tempBetaPower = pfcMLB.trialBetaPowerMtx(:,iscLog);
        betaPower{op,fl} = tempBetaPower;
        
        thetaPowerLog = tempThetaPower>powerThresh;
        betaPowerLog = tempBetaPower>powerThresh;
        
        figure(o);
        thetaDecodePhaseOdor = histcounts2(tempThetaPhase(thetaPowerLog), tempISCdecodeOdor(thetaPowerLog), phaseBins, odorBins);
        if nbp == 1
            thetaDecodePhaseOdor = thetaDecodePhaseOdor./repmat(max(thetaDecodePhaseOdor,[],2), [1,4]);
        elseif nbp == 0
            thetaDecodePhaseOdor = thetaDecodePhaseOdor./repmat(max(thetaDecodePhaseOdor), [size(thetaDecodePhaseOdor,1), 1]);
        end
        subplot(4,4,op)
        imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, odorBins(2:end)-0.5, thetaDecodePhaseOdor');
        set(gca, 'ydir', 'normal');
        ylabel('Odor');
        xlabel('\theta Phase');
        title(op);
        thetaDecodePowerOdor = histcounts2(tempThetaPower(:), tempISCdecodeOdor(:), powerBins, odorBins);
        if nbp == 1
            thetaDecodePowerOdor = thetaDecodePowerOdor./repmat(max(thetaDecodePowerOdor,[],2), [1,4]);
        elseif nbp == 0
            thetaDecodePowerOdor = thetaDecodePowerOdor./repmat(max(thetaDecodePowerOdor), [size(thetaDecodePowerOdor,1), 1]);
        end
        subplot(4,4,op+4)
        imagesc(powerBins(2:end)-mode(diff(powerBins))/2, odorBins(2:end)-0.5, thetaDecodePowerOdor');
        set(gca, 'ydir', 'normal');
        ylabel('Odor');
        xlabel('\theta Power');
        betaDecodePhaseOdor = histcounts2(tempBetaPhase(betaPowerLog), tempISCdecodeOdor(betaPowerLog), phaseBins, odorBins);
        if nbp == 1
            betaDecodePhaseOdor = betaDecodePhaseOdor./repmat(max(betaDecodePhaseOdor,[],2), [1,4]);
        elseif nbp == 0
            betaDecodePhaseOdor = betaDecodePhaseOdor./repmat(max(betaDecodePhaseOdor), [size(betaDecodePhaseOdor,1), 1]);
        end
        subplot(4,4,op+8)
        imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, odorBins(2:end)-0.5, betaDecodePhaseOdor');
        set(gca, 'ydir', 'normal');
        ylabel('Odor');
        xlabel('\beta Phase');
        betaDecodePowerOdor = histcounts2(tempBetaPower(:), tempISCdecodeOdor(:), powerBins, odorBins);
        if nbp == 1
            betaDecodePowerOdor = betaDecodePowerOdor./repmat(max(betaDecodePowerOdor,[],2), [1,4]);
        elseif nbp == 0
            betaDecodePowerOdor = betaDecodePowerOdor./repmat(max(betaDecodePowerOdor), [size(betaDecodePowerOdor,1), 1]);
        end
        subplot(4,4,op+12)
        imagesc(powerBins(2:end)-mode(diff(powerBins))/2, odorBins(2:end)-0.5, betaDecodePowerOdor');
        set(gca, 'ydir', 'normal');
        ylabel('Odor');
        xlabel('\beta Power');
        drawnow;
        
        figure(t);
        thetaDecodePhaseTime = histcounts2(tempThetaPhase(thetaPowerLog), tempISCdecodeTime(thetaPowerLog), phaseBins, timeBins);
        if nbp == 1
            thetaDecodePhaseTime = thetaDecodePhaseTime./repmat(max(thetaDecodePhaseTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
        elseif nbp == 0
            thetaDecodePhaseTime = thetaDecodePhaseTime./repmat(max(thetaDecodePhaseTime), [size(thetaDecodePhaseTime,1), 1]);
        end
        subplot(4,4,op)
        imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, thetaDecodePhaseTime');
        set(gca, 'ydir', 'normal');
        ylabel('Time');
        xlabel('\theta Phase');
        title(op);
        thetaDecodePowerTime = histcounts2(tempThetaPower(:), tempISCdecodeTime(:), powerBins, timeBins);
        if nbp == 1
            thetaDecodePowerTime = thetaDecodePowerTime./repmat(max(thetaDecodePowerTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
        elseif nbp == 0
            thetaDecodePowerTime = thetaDecodePowerTime./repmat(max(thetaDecodePowerTime), [size(thetaDecodePowerTime,1), 1]);
        end
        subplot(4,4,op+4)
        imagesc(powerBins(2:end)-mode(diff(powerBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, thetaDecodePowerTime');
        set(gca, 'ydir', 'normal');
        ylabel('Time');
        xlabel('\theta Power');
        betaDecodePhaseTime = histcounts2(tempBetaPhase(betaPowerLog), tempISCdecodeTime(betaPowerLog), phaseBins, timeBins);
        if nbp == 1
            betaDecodePhaseTime = betaDecodePhaseTime./repmat(max(betaDecodePhaseTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
        elseif nbp == 0
            betaDecodePhaseTime = betaDecodePhaseTime./repmat(max(betaDecodePhaseTime), [size(betaDecodePhaseTime,1), 1]);
        end
        subplot(4,4,op+8)
        imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, betaDecodePhaseTime');
        set(gca, 'ydir', 'normal');
        ylabel('Time');
        xlabel('\beta Phase');
        betaDecodePowerTime = histcounts2(tempBetaPower(:), tempISCdecodeTime(:), powerBins, timeBins);
        if nbp == 1
            betaDecodePowerTime = betaDecodePowerTime./repmat(max(betaDecodePowerTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
        elseif nbp == 0
            betaDecodePowerTime = betaDecodePowerTime./repmat(max(betaDecodePowerTime), [size(betaDecodePowerTime,1), 1]);
        end
        subplot(4,4,op+12)
        imagesc(powerBins(2:end)-mode(diff(powerBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, betaDecodePowerTime');
        set(gca, 'ydir', 'normal');
        ylabel('Time');
        xlabel('\beta Power');
        drawnow;
    end
    if nbp == 1
        figure(o);
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Power/Phase', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
        figure(t);
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Power/Phase', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    elseif nbp == 0
        figure(o);
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Odor', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
        figure(t)
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Odor', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    else
        figure(o);
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, No Norm', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
        figure(t)
        annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, No Norm', binSize, dsRate, powerThresh),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    end
    figure(o);
    annotation(gcf,'textbox', [0 0.05 1 0.05],'String', fileDirs{fl},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    figure(t);
    annotation(gcf,'textbox', [0 0.05 1 0.05],'String', fileDirs{fl},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     pfcMLB = [];
end

%%
o = figure('Position', [88 350 560 420]);
t = figure('Position', [1248 350 560 420]);
for op = 1:4
    grpOdor = cell2mat(iscDecodesOdor(op,:));
    grpTime = cell2mat(iscDecodesTime(op,:));
    grpThetaPhase = cell2mat(thetaPhase(op,:));
    grpThetaPower = cell2mat(thetaPower(op,:));
    grpBetaPhase = cell2mat(betaPhase(op,:));
    grpBetaPower = cell2mat(betaPower(op,:));
    
    thetaPowerLog = grpThetaPower>powerThresh;
    betaPowerLog = grpBetaPower>powerThresh;
    
    figure(o);
    thetaDecodePhaseOdor = histcounts2(grpThetaPhase(thetaPowerLog), grpOdor(thetaPowerLog), phaseBins, odorBins);
    if nbp == 1
        thetaDecodePhaseOdor = thetaDecodePhaseOdor./repmat(max(thetaDecodePhaseOdor,[],2), [1,4]);
    elseif nbp == 0
        thetaDecodePhaseOdor = thetaDecodePhaseOdor./repmat(max(thetaDecodePhaseOdor), [size(thetaDecodePhaseOdor,1), 1]);
    end
    subplot(4,4,op)
    imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, odorBins(2:end)-0.5, thetaDecodePhaseOdor');
    set(gca, 'ydir', 'normal');
    ylabel('Odor');
    xlabel('\theta Phase');
    title(op);
    thetaDecodePowerOdor = histcounts2(grpThetaPower(:), grpOdor(:), powerBins, odorBins);
    if nbp == 1
        thetaDecodePowerOdor = thetaDecodePowerOdor./repmat(max(thetaDecodePowerOdor,[],2), [1,4]);
    elseif nbp == 0
        thetaDecodePowerOdor = thetaDecodePowerOdor./repmat(max(thetaDecodePowerOdor), [size(thetaDecodePowerOdor,1), 1]);
    end
    subplot(4,4,op+4)
    imagesc(powerBins(2:end)-mode(diff(powerBins))/2, odorBins(2:end)-0.5, thetaDecodePowerOdor');
    set(gca, 'ydir', 'normal');
    ylabel('Odor');
    xlabel('\theta Power');
    betaDecodePhaseOdor = histcounts2(grpBetaPhase(betaPowerLog), grpOdor(betaPowerLog), phaseBins, odorBins);
    if nbp == 1
        betaDecodePhaseOdor = betaDecodePhaseOdor./repmat(max(betaDecodePhaseOdor,[],2), [1,4]);
    elseif nbp == 0
        betaDecodePhaseOdor = betaDecodePhaseOdor./repmat(max(betaDecodePhaseOdor), [size(betaDecodePhaseOdor,1), 1]);
    end
    subplot(4,4,op+8)
    imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, odorBins(2:end)-0.5, betaDecodePhaseOdor');
    set(gca, 'ydir', 'normal');
    ylabel('Odor');
    xlabel('\beta Phase');
    betaDecodePowerOdor = histcounts2(grpBetaPower(:), grpOdor(:), powerBins, odorBins);
    if nbp == 1
        betaDecodePowerOdor = betaDecodePowerOdor./repmat(max(betaDecodePowerOdor,[],2), [1,4]);
    elseif nbp == 0
        betaDecodePowerOdor = betaDecodePowerOdor./repmat(max(betaDecodePowerOdor), [size(betaDecodePowerOdor,1), 1]);
    end
    subplot(4,4,op+12)
    imagesc(powerBins(2:end)-mode(diff(powerBins))/2, odorBins(2:end)-0.5, betaDecodePowerOdor');
    set(gca, 'ydir', 'normal');
    ylabel('Odor');
    xlabel('\beta Power');
    drawnow;
    
    figure(t);
    thetaDecodePhaseTime = histcounts2(grpThetaPhase(thetaPowerLog), grpTime(thetaPowerLog), phaseBins, timeBins);
    if nbp == 1
        thetaDecodePhaseTime = thetaDecodePhaseTime./repmat(max(thetaDecodePhaseTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
    elseif nbp == 0
        thetaDecodePhaseTime = thetaDecodePhaseTime./repmat(max(thetaDecodePhaseTime), [size(thetaDecodePhaseTime,1), 1]);
    end
    subplot(4,4,op)
    imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, thetaDecodePhaseTime');
    set(gca, 'ydir', 'normal');
    ylabel('Time');
    xlabel('\theta Phase');
    title(op);
    thetaDecodePowerTime = histcounts2(grpThetaPower(:), grpTime(:), powerBins, timeBins);
    if nbp == 1
        thetaDecodePowerTime = thetaDecodePowerTime./repmat(max(thetaDecodePowerTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
    elseif nbp == 0
        thetaDecodePowerTime = thetaDecodePowerTime./repmat(max(thetaDecodePowerTime), [size(thetaDecodePowerTime,1), 1]);
    end
    subplot(4,4,op+4)
    imagesc(powerBins(2:end)-mode(diff(powerBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, thetaDecodePowerTime');
    set(gca, 'ydir', 'normal');
    ylabel('Time');
    xlabel('\theta Power');
    betaDecodePhaseTime = histcounts2(grpBetaPhase(betaPowerLog), grpTime(betaPowerLog), phaseBins, timeBins);
    if nbp == 1
        betaDecodePhaseTime = betaDecodePhaseTime./repmat(max(betaDecodePhaseTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
    elseif nbp == 0
        betaDecodePhaseTime = betaDecodePhaseTime./repmat(max(betaDecodePhaseTime), [size(betaDecodePhaseTime,1), 1]);
    end
    subplot(4,4,op+8)
    imagesc(phaseBins(2:end)-mode(diff(phaseBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, betaDecodePhaseTime');
    set(gca, 'ydir', 'normal');
    ylabel('Time');
    xlabel('\beta Phase');
    betaDecodePowerTime = histcounts2(grpBetaPower(:), grpTime(:), powerBins, timeBins);
    if nbp == 1
        betaDecodePowerTime = betaDecodePowerTime./repmat(max(betaDecodePowerTime,[],2), [1,size(thetaDecodePhaseTime,2)]);
    elseif nbp == 0
        betaDecodePowerTime = betaDecodePowerTime./repmat(max(betaDecodePowerTime), [size(betaDecodePowerTime,1), 1]);
    end
    subplot(4,4,op+12)
    imagesc(powerBins(2:end)-mode(diff(powerBins))/2, timeBins(2:end)-mode(diff(timeBins))/2, betaDecodePowerTime');
    set(gca, 'ydir', 'normal');
    ylabel('Time');
    xlabel('\beta Power');
    drawnow;
end
if nbp == 1
    figure(o);
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Power/Phase', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    figure(t);
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Power/Phase', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
elseif nbp == 0
    figure(o);
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Odor', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    figure(t)
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, Norm by Odor', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
else
    figure(o);
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, No Norm', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    figure(t)
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('BinSize = %i, DSrate = %i, PowerThresh = %i, No Norm', binSize, dsRate, powerThresh),...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end
figure(o);
annotation(gcf,'textbox', [0 0.05 1 0.05],'String', 'Group Data',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
figure(t);
annotation(gcf,'textbox', [0 0.05 1 0.05],'String', 'Group Data',...
    'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

