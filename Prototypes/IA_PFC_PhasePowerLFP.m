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

trialWindow = [0 1200];
alignment = 'PokeIn';

phaseBins = -pi:pi/8:pi;
powerBins = -2:0.5:6;

%#ok<*AGROW>
%%
for ani = 1:length(fileDirs)
    data = MLB_SM(fileDirs{ani});
    tetWMU = find(data.numUniPerTet==max(data.numUniPerTet), 1, 'first');
    
    sessionWindow = [data.trialInfo(1).PokeInIndex-1000; data.trialInfo(end).PokeOutIndex+2000]; % To remove the pre/post session recordings from merged files that were cut to have unit IDs across recordings
    
    [thetaFilt, thetaPhase, thetaPwr] = data.SimpleFilter(data.lfpMatrix(sessionWindow(1):sessionWindow(2),tetWMU), [4 12]);
    thetaFilt = [zeros(sessionWindow(1)-1,1); thetaFilt; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)]; 
    thetaPhase = [zeros(sessionWindow(1)-1,1); thetaPhase; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    thetaPwr = [zeros(sessionWindow(1)-1,1); thetaPwr; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    [betaFilt, betaPhase, betaPwr] = data.SimpleFilter(data.lfpMatrix(sessionWindow(1):sessionWindow(2),tetWMU), [16 32]);
    betaFilt = [zeros(sessionWindow(1)-1,1); betaFilt; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    betaPhase = [zeros(sessionWindow(1)-1,1); betaPhase; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    betaPwr = [zeros(sessionWindow(1)-1,1); betaPwr; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    [gammaFilt, gammaPhase, gammaPwr] = data.SimpleFilter(data.lfpMatrix(sessionWindow(1):sessionWindow(2),tetWMU), [40 100]);
    gammaFilt = [zeros(sessionWindow(1)-1,1); gammaFilt; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    gammaPhase = [zeros(sessionWindow(1)-1,1); gammaPhase; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    gammaPwr = [zeros(sessionWindow(1)-1,1); gammaPwr; zeros(size(data.lfpMatrix,1)-sessionWindow(2),1)];
    
    thetaFiltTrl = data.ExtractTrialMatrix(thetaFilt, trialWindow, alignment);
    thetaPhaseTrl = data.ExtractTrialMatrix(thetaPhase, trialWindow, alignment);
    thetaPwrTrl = data.ExtractTrialMatrix(thetaPwr, trialWindow, alignment);
    
    betaFiltTrl = data.ExtractTrialMatrix(betaFilt, trialWindow, alignment);
    betaPhaseTrl = data.ExtractTrialMatrix(betaPhase, trialWindow, alignment);
    betaPwrTrl = data.ExtractTrialMatrix(betaPwr, trialWindow, alignment);
    
    gammaFiltTrl = data.ExtractTrialMatrix(gammaFilt, trialWindow, alignment);
    gammaPhaseTrl = data.ExtractTrialMatrix(gammaPhase, trialWindow, alignment);
    gammaPwrTrl = data.ExtractTrialMatrix(gammaPwr, trialWindow, alignment);
    
    trialLFP = data.ExtractTrialMatrix(data.lfpMatrix(:,tetWMU), trialWindow, alignment);
    %%
%     trl = 0;
%     trl = trl+1;
%     tempRaw = trialLFP(:,1,trl);
%     tempBetaFilt = betaFiltTrl(:,1,trl);
%     tempBetaPhase = betaPhaseTrl(:,1,trl)./pi;
%     tempBetaPwr = betaPwrTrl(:,1,trl);
%     tempUnwrappedBetaPhase = unwrap(betaPhaseTrl(:,1,trl));
%     
%     tempThetaFilt = thetaFiltTrl(:,1,trl);
%     tempThetaPhase = thetaPhaseTrl(:,1,trl)./pi;
%     tempThetaPwr = thetaPwrTrl(:,1,trl);
%     tempUnwrappedThetaPhase = unwrap(thetaPhaseTrl(:,1,trl));
%     
%     tempGammaFilt = gammaFiltTrl(:,1,trl);
%     tempGammaPhase = gammaPhaseTrl(:,1,trl)./pi;
%     tempGammaPwr = gammaPwrTrl(:,1,trl);
%     tempUnwrappedGammaPhase = unwrap(gammaPhaseTrl(:,1,trl));
    %%
%     figure;
%     subplot(2,1,1);
%     plot(tempRaw, 'ydatasource', 'tempRaw', 'color', 'k');
%     hold on;
%     plot(tempThetaFilt, 'ydatasource', 'tempThetaFilt');
%     plot(tempThetaPwr, 'ydatasource', 'tempThetaPwr', 'linestyle', '-.');
%     % plot(tempThetaPhase, 'ydatasource', 'tempThetaPhase', 'linestyle', ':')
%     plot(tempBetaFilt, 'ydatasource', 'tempBetaFilt');
%     plot(tempBetaPwr, 'ydatasource', 'tempBetaPwr', 'linestyle', '-.');
%     % plot(tempBetaPhase, 'ydatasource', 'tempBetaPhase', 'linestyle', ':');
%     plot(tempGammaFilt, 'ydatasource', 'tempGammaFilt');
%     plot(tempGammaPwr, 'ydatasource', 'tempGammaPwr', 'linestyle', '-.');
%     % plot(tempGammaPhase, 'ydatasource', 'tempGammaPhase', 'linestyle', ':')
%     legend('Raw', '\theta', '\theta-Power', '\beta', '\beta-Power', '\gamma', '\gamma-Power');
%     % legend('Raw', '\theta', '\theta-Power', '\theta-Phase', '\beta', '\beta-Power',  '\beta-Phase', '\gamma', '\gamma-Power', '\gamma-Phase');
%     subplot(2,1,2);
%     plot(tempUnwrappedThetaPhase, 'ydatasource', 'tempUnwrappedThetaPhase');
%     hold on;
%     plot(tempUnwrappedBetaPhase, 'ydatasource', 'tempUnwrappedBetaPhase');
%     plot(tempUnwrappedGammaPhase, 'ydatasource', 'tempUnwrappedGammaPhase');
%     legend('Unwrapped \theta', 'Unwrapped \beta', 'Unwrapped \gamma');
%     linkaxes([subplot(2,1,1), subplot(2,1,2)], 'x');
%     
%     annotation(gcf,'textbox', [0 0.95 1 0.05],'String', fileDirs{ani},...
%         'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     
    %%
    iscLog = ([data.trialInfo.Odor] == [data.trialInfo.Position]) &...
        ([data.trialInfo.Performance] == 1);
    
    iscThetaPhase = thetaPhaseTrl(:,:,iscLog);
    iscBetaPhase = betaPhaseTrl(:,:,iscLog);
    iscGammaPhase = gammaPhaseTrl(:,:,iscLog);
    
    iscThetaPower = thetaPwrTrl(:,:,iscLog);
    iscBetaPower = betaPwrTrl(:,:,iscLog);
    iscGammaPower = gammaPwrTrl(:,:,iscLog);
    
    numThetaCycsPrTrl = nan(1,sum(iscLog));
    numBetaCycsPrTrl = nan(1,sum(iscLog));
    numGammaCycsPrTrl = nan(1,sum(iscLog));
    
    thetaTrlCycPwr = nan(30,sum(iscLog));
    betaTrlCycPwr = nan(75,sum(iscLog));
    gammaTrlCycPwr = nan(150,sum(iscLog));
    
    for trl = 1:sum(iscLog)
        curTrlThetaPhase = iscThetaPhase(:,:,trl);
        [~, thetaPhaseBounds] = findpeaks(diff(curTrlThetaPhase)*-1, 'MinPeakProminence', 5);
        numThetaCycsPrTrl(trl) = length(thetaPhaseBounds);
        thetaPhaseBounds = [1; thetaPhaseBounds];
        for tPB = 2:length(thetaPhaseBounds)
            thetaTrlCycPwr(tPB-1,trl) = mean(iscThetaPower(thetaPhaseBounds(tPB-1):thetaPhaseBounds(tPB),:,trl));
        end
        
        curTrlBetaPhase = iscBetaPhase(:,:,trl);
        [~, betaPhaseBounds] = findpeaks(diff(curTrlBetaPhase)*-1, 'MinPeakProminence', 5);
        numBetaCycsPrTrl(trl) = length(betaPhaseBounds);
        betaPhaseBounds = [1; betaPhaseBounds];
        for bPB = 2:length(betaPhaseBounds)
            betaTrlCycPwr(bPB-1,trl) = mean(iscBetaPower(betaPhaseBounds(bPB-1):betaPhaseBounds(bPB),:,trl));
        end
        
        curTrlGammaPhase = iscGammaPhase(:,:,trl);
        [~, gammaPhaseBounds] = findpeaks(diff(curTrlGammaPhase)*-1, 'MinPeakProminence', 5);
        numGammaCycsPrTrl(trl) = length(gammaPhaseBounds);
        gammaPhaseBounds = [1; gammaPhaseBounds];
        for gPB = 2:length(gammaPhaseBounds)
            gammaTrlCycPwr(gPB-1,trl) = mean(iscGammaPower(gammaPhaseBounds(gPB-1):gammaPhaseBounds(gPB),:,trl));
        end        
    end
    figure;
    subplot(2,3,1)
    histogram(numThetaCycsPrTrl, 'Normalization', 'pdf');
    subplot(2,3,2)
    histogram(numBetaCycsPrTrl, 'Normalization', 'pdf');
    subplot(2,3,3)
    histogram(numGammaCycsPrTrl, 'Normalization', 'pdf');
    subplot(2,3,4)
    imagesc(thetaTrlCycPwr, [-4 4]); set(gca, 'ydir', 'normal'); xlabel('Trial'); ylabel('Cycle');
    subplot(2,3,5)
    imagesc(betaTrlCycPwr, [-4 4]); set(gca, 'ydir', 'normal'); xlabel('Trial'); ylabel('Cycle');
    subplot(2,3,6)
    imagesc(gammaTrlCycPwr, [-4 4]); set(gca, 'ydir', 'normal'); xlabel('Trial'); ylabel('Cycle');

    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    %%
    trlThetaPwr = reshape(thetaPwrTrl(:,1,data.fiscTrials(:)), [size(thetaPwrTrl,1),numel(data.fiscTrials)]);
    trlThetaPhase = reshape(thetaPhaseTrl(:,1,data.fiscTrials(:)), [size(thetaPhaseTrl,1),numel(data.fiscTrials)]);
    
    trlBetaPwr = reshape(betaPwrTrl(:,1,data.fiscTrials(:)), [size(betaPwrTrl,1),numel(data.fiscTrials)]);
    trlBetaPhase = reshape(betaPhaseTrl(:,1,data.fiscTrials(:)), [size(betaPhaseTrl,1),numel(data.fiscTrials)]);
    
    trlGammaPwr = reshape(gammaPwrTrl(:,1,data.fiscTrials(:)), [size(gammaPwrTrl,1),numel(data.fiscTrials)]);
    trlGammaPhase = reshape(gammaPhaseTrl(:,1,data.fiscTrials(:)), [size(gammaPhaseTrl,1),numel(data.fiscTrials)]);
    
    figure;
    subplot(3,3,2)
    [tB, tBx, tBy] = histcounts2(trlThetaPhase(:), trlBetaPwr(:), phaseBins, powerBins);
    imagesc(tBx(1:end-1)-pi/80, tBy(1:end-1)-0.25, tB');
    set(gca, 'ydir', 'normal');
    xlabel('\theta Phase');
    ylabel('\beta Power');
    subplot(3,3,3)
    [tG, tGx, tGy] = histcounts2(trlThetaPhase(:), trlGammaPwr(:), phaseBins, powerBins);
    imagesc(tGx(1:end-1)-pi/80, tGy(1:end-1)-0.25, tG');
    set(gca, 'ydir', 'normal');
    xlabel('\theta Phase');
    ylabel('\gamma Power');
    
    subplot(3,3,4)
    [bT, bTx, bTy] = histcounts2(trlBetaPhase(:), trlThetaPwr(:), phaseBins, powerBins);
    imagesc(bTx(1:end-1)-pi/80, bTy(1:end-1)-0.25, bT');
    set(gca, 'ydir', 'normal');
    xlabel('\beta Phase');
    ylabel('\theta Power');
    subplot(3,3,6)
    [bG, bGx, bGy] = histcounts2(trlBetaPhase(:), trlGammaPwr(:), phaseBins, powerBins);
    imagesc(bGx(1:end-1)-pi/80, bGy(1:end-1)-0.25, bG');
    set(gca, 'ydir', 'normal');
    xlabel('\beta Phase');
    ylabel('\gamma Power')
    
    subplot(3,3,7)
    [gT, gTx, gTy] = histcounts2(trlGammaPhase(:), trlThetaPwr(:), phaseBins, powerBins);
    imagesc(gTx(1:end-1)-pi/80, gTy(1:end-1)-0.25, gT');
    set(gca, 'ydir', 'normal');
    xlabel('\gamma Phase');
    ylabel('\theta Power');
    subplot(3,3,8)
    [gB, gBx, gBy] = histcounts2(trlGammaPhase(:), trlBetaPwr(:), phaseBins, powerBins);
    imagesc(gBx(1:end-1)-pi/80, gBy(1:end-1)-0.25, gB');
    set(gca, 'ydir', 'normal');
    xlabel('\gamma Phase');
    ylabel('\beta Power');
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    %
    figure;
    subplot(3,3,1);
    histogram(trlThetaPwr(:));
    xlabel('\theta Power');
    subplot(3,3,2)
    [tB, tBx, tBy] = histcounts2(trlThetaPwr(:), trlBetaPwr(:), powerBins, powerBins);
    imagesc(tBx(1:end-1)-pi/80, tBy(1:end-1)-0.25, tB');
    set(gca, 'ydir', 'normal');
    xlabel('\theta Power');
    ylabel('\beta Power');
    subplot(3,3,3)
    [tG, tGx, tGy] = histcounts2(trlThetaPwr(:), trlGammaPwr(:), powerBins, powerBins);
    imagesc(tGx(1:end-1)-pi/80, tGy(1:end-1)-0.25, tG');
    set(gca, 'ydir', 'normal');
    xlabel('\theta Power');
    ylabel('\gamma Power');
    
    subplot(3,3,4)
    [bT, bTx, bTy] = histcounts2(trlBetaPwr(:), trlThetaPwr(:), powerBins, powerBins);
    imagesc(bTx(1:end-1)-pi/80, bTy(1:end-1)-0.25, bT');
    set(gca, 'ydir', 'normal');
    xlabel('\beta Power');
    ylabel('\theta Power');
    subplot(3,3,5);
    histogram(trlBetaPwr(:));
    xlabel('\beta Power');
    subplot(3,3,6)
    [bG, bGx, bGy] = histcounts2(trlBetaPwr(:), trlGammaPwr(:), powerBins, powerBins);
    imagesc(bGx(1:end-1)-pi/80, bGy(1:end-1)-0.25, bG');
    set(gca, 'ydir', 'normal');
    xlabel('\beta Power');
    ylabel('\gamma Power')
    
    subplot(3,3,7)
    [gT, gTx, gTy] = histcounts2(trlGammaPwr(:), trlThetaPwr(:), powerBins, powerBins);
    imagesc(gTx(1:end-1)-pi/80, gTy(1:end-1)-0.25, gT');
    set(gca, 'ydir', 'normal');
    xlabel('\gamma Power');
    ylabel('\theta Power');
    subplot(3,3,8)
    [gB, gBx, gBy] = histcounts2(trlGammaPwr(:), trlBetaPwr(:), powerBins, powerBins);
    imagesc(gBx(1:end-1)-pi/80, gBy(1:end-1)-0.25, gB');
    set(gca, 'ydir', 'normal');
    xlabel('\gamma Power');
    ylabel('\beta Power');
    subplot(3,3,9);
    histogram(trlGammaPwr(:));
    xlabel('\gamma Power');
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', fileDirs{ani},...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    drawnow
end
