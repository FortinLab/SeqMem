classdef SingleUnit_SM < SeqMem
    % Single unit analyses from the Sequence Memory project
    properties % General Variables
        gaussWinDur = 100
    end
    %% Creation Method
    methods
        %% Main Creation Method
        function obj = SingleUnit_SM(fileDir)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@SeqMem(fileDir);        
        end
    end
    %% General Method
    methods
        %% Bin trial event spiking
        function binnedSpikes = BinTrialEventSpikes(obj,window,alignment)
            trialSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix, window, alignment);
            binnedSpikes = sum(trialSpikes,1);
        end
        %% Extract Trial Types
        function [trlSpikes, trlIDs] = ExtractTrialSpikes(obj, spiking, trlType, varargin)
            trlDim = find(size(spiking)==length(obj.trialInfo));
            trlNum = [obj.trialInfo.TrialNum];
            posVect = [obj.trialInfo.Position];
            odrVect = [obj.trialInfo.Odor];
            perfLog = [obj.trialInfo.Performance]==1;
            isLog = [obj.trialInfo.TranspositionDistance]==0;
            taoLog = [obj.trialInfo.TranspositionDistance]==0 & [obj.trialInfo.ItemItemDistance]~=1;
            ssnLogs = [trlNum;posVect;odrVect;perfLog;isLog;taoLog];
            if strcmp(trlType,'all')
                trlSpikes = spiking;
                trlIDs = ssnLogs;
            elseif strcmp(trlType,'isc')
                iscLog = isLog & perfLog;
                if trlDim == 1
                    trlSpikes = spiking(iscLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,iscLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,iscLog);
                end
                trlIDs = ssnLogs(:,iscLog);
            elseif strcmp(trlType,'osc')
                oscLog = ~isLog & perfLog;
                if trlDim == 1
                    trlSpikes = spiking(oscLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,oscLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,oscLog);
                end
                trlIDs = ssnLogs(:,oscLog);
            elseif strcmp(trlType, 'isi')
                isiLog = isLog & ~perfLog;
                if trlDim == 1
                    trlSpikes = spiking(isiLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,isiLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,isiLog);
                end
                trlIDs = ssnLogs(:,isiLog);
            elseif strcmp(trlType, 'corr')
                if trlDim == 1
                    trlSpikes = spiking(perfLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,perfLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,perfLog);
                end
                trlIDs = ssnLogs(:,perfLog);
            elseif strcmp(trlType, 'incorr')
                if trlDim == 1
                    trlSpikes = spiking(~perfLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,~perfLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,~perfLog);
                end
                trlIDs = ssnLogs(:,~perfLog);
            end
        end
    end
    %% Specific Analyses
    methods
        %% Compare ISC trial period spiking
        function [anovaTable, anovaStats, anovaData, groupIDs] = TrialPeriodSpiking(obj)
%             ssnLogs = [trlNum;posVect;odrVect;perfLog;isLog;taoLog];
            preTrialSpikes = obj.BinTrialEventSpikes([-500 0], 'PokeIn');
            [iscPRE, trlIDs] = obj.ExtractTrialSpikes(squeeze(preTrialSpikes), 'isc');
            rlyTrialSpikes = obj.BinTrialEventSpikes([0 500], 'PokeIn');
            iscRLY = obj.ExtractTrialSpikes(squeeze(rlyTrialSpikes), 'isc');
            latTrialSpikes = obj.BinTrialEventSpikes([-500 0], 'PokeOut');
            iscLAT = obj.ExtractTrialSpikes(squeeze(latTrialSpikes), 'isc');
            pstTrialSpikes = obj.BinTrialEventSpikes([0 500], 'PokeOut');
            iscPST = obj.ExtractTrialSpikes(squeeze(pstTrialSpikes), 'isc');

            posIDs = repmat(trlIDs(2,:)',[4,1]);
            timeIDs = [ones(size(iscPRE,2),1); ones(size(iscRLY,2),1)+1; ones(size(iscLAT,2),1)+2; ones(size(iscPST,2),1)+3];
            groupIDs = {timeIDs, posIDs};
            anovaTable = cell(6,7,length(obj.ensembleMatrixColIDs));
            anovaStats = cell(1,length(obj.ensembleMatrixColIDs));
            anovaData = cell(1,length(obj.ensembleMatrixColIDs));
            for uni = 1:length(obj.ensembleMatrixColIDs)
                tempAnovaData = [iscPRE(uni,:)';iscRLY(uni,:)';iscLAT(uni,:)';iscPST(uni,:)'];
                anovaData{uni} = tempAnovaData;
                [~,anovaTable(:,:,uni), anovaStats{uni}] = anovan(tempAnovaData, groupIDs, 'model', 'interaction', 'varnames', {'Time', 'Position'}, 'display', 'off');
            end
        end
        %% Compare Reward vs Error Spiking
        function [anovaTable, anovaStats, anovaData, groupIDs] = FeedbackSpiking(obj)
            rwdSpikes = obj.BinTrialEventSpikes([-250 250], 'FrontReward');
            [iscRWD, rwdTrls] = obj.ExtractTrialSpikes(squeeze(rwdSpikes), 'corr');
            errSpikes = obj.BinTrialEventSpikes([-250 250], 'ErrorSignal');
            [iscERR, errTrls] = obj.ExtractTrialSpikes(squeeze(errSpikes), 'incorr');
            
            posIDs = [rwdTrls(2,:)'; errTrls(2,:)'];
            eventIDs = [ones(size(iscRWD,2),1); ones(size(iscERR,2),1)+1];
            groupIDs = {eventIDs, posIDs};
            anovaTable = cell(6,7,length(obj.ensembleMatrixColIDs));
            anovaStats = cell(1,length(obj.ensembleMatrixColIDs));
            anovaData = cell(1,length(obj.ensembleMatrixColIDs));
            for uni = 1:length(obj.ensembleMatrixColIDs)
                tempAnovaData = [iscRWD(uni,:)';iscERR(uni,:)'];
                anovaData{uni} = tempAnovaData;
                [~,anovaTable(:,:,uni), anovaStats{uni}] = anovan(tempAnovaData, groupIDs, 'model', 'interaction', 'varnames', {'Event','Position'}, 'display', 'off');
            end
        end
        %% Evaluate Error Responses
        function [anovaTable, anovaStats, anovaData, groupIDs] = ErrorSpiking(obj)
            preErrSigSpikes = obj.BinTrialEventSpikes([-250 0], 'ErrorSignal');
            [preERR, trialIDs] = obj.ExtractTrialSpikes(squeeze(preErrSigSpikes), 'incorr');
            pstErrSigSpikes = obj.BinTrialEventSpikes([0 250], 'ErrorSignal');
            pstERR = obj.ExtractTrialSpikes(squeeze(pstErrSigSpikes), 'incorr');
            
            posIDs = repmat(trialIDs(2,:)', [2,1]);
            timeIDs = [ones(size(preERR,2),1); ones(size(pstERR,2),1)+1];
            groupIDs = {timeIDs, posIDs};
            anovaTable = cell(6,7,length(obj.ensembleMatrixColIDs));
            anovaStats = cell(1,length(obj.ensembleMatrixColIDs));
            anovaData = cell(1,length(obj.ensembleMatrixColIDs));
            for uni = 1:length(obj.ensembleMatrixColIDs)
                tempAnovaData = [preERR(uni,:)';pstERR(uni,:)'];
                anovaData{uni} = tempAnovaData;
                [~,anovaTable(:,:,uni), anovaStats{uni}] = anovan(tempAnovaData, groupIDs, 'model', 'interaction', 'varnames', {'Time','Position'}, 'display', 'off');
            end
        end
        %% Evaluate Error Response Effects Sequentially
        function [tStats, anovaTable, anovaStats] = ErrorSpikingSequential(obj)
            preErrSigSpikes = obj.BinTrialEventSpikes([-250 0], 'ErrorSignal');
            [preERR, trialIDs] = obj.ExtractTrialSpikes(squeeze(preErrSigSpikes), 'incorr');
            pstErrSigSpikes = obj.BinTrialEventSpikes([0 250], 'ErrorSignal');
            pstERR = obj.ExtractTrialSpikes(squeeze(pstErrSigSpikes), 'incorr');
            
            tp = cell(1,length(obj.ensembleMatrixColIDs));
            ts = cell(1,length(obj.ensembleMatrixColIDs));
            anovaTable = cell(4,6,length(obj.ensembleMatrixColIDs));
            anovaStats = cell(1,length(obj.ensembleMatrixColIDs));
            for uni = 1:length(obj.ensembleMatrixColIDs)
                [~,tp{uni},~,ts{uni}] = ttest(preERR(uni,:),pstERR(uni,:));
                [~,anovaTable(:,:,uni), anovaStats{uni}] = anova1(pstERR(uni,:)', trialIDs(2,:)', 'off');
            end
            tStats = struct('p', tp, 'stats', ts);
        end
        %% Evaluate Reward Responses
        function [anovaTable, anovaStats, anovaData, groupIDs] = RewardSpiking(obj)
            preRwdSigSpikes = obj.BinTrialEventSpikes([-250 0], 'FrontReward');
            [preRWD, trialIDs] = obj.ExtractTrialSpikes(squeeze(preRwdSigSpikes), 'corr');
            pstRwdSigSpikes = obj.BinTrialEventSpikes([0 250], 'FrontReward');
            pstRWD = obj.ExtractTrialSpikes(squeeze(pstRwdSigSpikes), 'corr');
            
            posIDs = repmat(trialIDs(2,:)', [2,1]);
            timeIDs = [ones(size(preRWD,2),1); ones(size(pstRWD,2),1)+1];
            groupIDs = {timeIDs, posIDs};
            anovaTable = cell(6,7,length(obj.ensembleMatrixColIDs));
            anovaStats = cell(1,length(obj.ensembleMatrixColIDs));
            anovaData = cell(1,length(obj.ensembleMatrixColIDs));
            for uni = 1:length(obj.ensembleMatrixColIDs)
                tempAnovaData = [preRWD(uni,:)';pstRWD(uni,:)'];
                anovaData{uni} = tempAnovaData;
                [~,anovaTable(:,:,uni), anovaStats{uni}] = anovan(tempAnovaData, groupIDs, 'model', 'interaction', 'varnames', {'Time','Position'}, 'display', 'off');
            end
        end
    end
    %% Visualizations
    methods
        %% Plot Unit Summary
        function PlotTrialSummary(obj, uni, window, alignment)
            if nargin==2
                window = [-1500 2000];
                alignment = 'PokeIn';
                uni = {uni};
            elseif nargin==1
                window = [-1500 2000];
                alignment = 'PokeIn';
                uni = obj.ensembleMatrixColIDs;
            else
                uni = {uni};
            end
            for u = 1:length(uni)
                figure;
                sp(1) = subplot(4,1,1);
                obj.PlotRastersByPos(window, alignment, uni{u}, 'isc');
                axis tight;
                sp(2) = subplot(4,1,2:4);
                obj.PlotGaussFireByPos(window, alignment, uni{u}, 'isc');
                axis tight;
                linkaxes(sp, 'x');
                drawnow;
            end
        end
        %% Plot Error Summary
        function PlotErrorSummary(obj, uni, window, alignment)
            if nargin==2
                window = [-250 250];
                alignment = 'ErrorSignal';
                uni = {uni};
            elseif nargin==1
                window = [-250 250];
                alignment = 'ErrorSignal';
                uni = obj.ensembleMatrixColIDs;
            else
                uni = {uni};
            end
            for u = 1:length(uni)
                figure;
                sp(1) = subplot(4,1,1);
                obj.PlotRastersByPos(window, alignment, uni{u}, 'incorr');
                axis tight;
                sp(2) = subplot(4,1,2:4);
                obj.PlotGaussFireByPos(window, alignment, uni{u}, 'incorr');
                axis tight;
                linkaxes(sp, 'x');
                drawnow;
            end
        end
        %% Plot Reward Summary
        function PlotRewardSummary(obj, uni, window, alignment)
            if nargin==2
                window = [-250 250];
                alignment = 'FrontReward';
                uni = {uni};
            elseif nargin==1
                window = [-250 250];
                alignment = 'FrontReward';
                uni = obj.ensembleMatrixColIDs;
            else
                uni = {uni};
            end
            for u = 1:length(uni)
                figure;
                sp(1) = subplot(4,1,1);
                obj.PlotRastersByPos(window, alignment, uni{u}, 'corr');
                axis tight;
                sp(2) = subplot(4,1,2:4);
                obj.PlotGaussFireByPos(window, alignment, uni{u}, 'corr');
                axis tight;
                linkaxes(sp, 'x');
                drawnow;
            end
        end
        %% Plot Rasters By Position
        function PlotRastersByPos(obj, window, alignment, uni, trlType)
            uniCol = strcmp(obj.ensembleMatrixColIDs, uni);
            eventSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix(:,uniCol), window, alignment);
            [evtSpks, trlIDs] = obj.ExtractTrialSpikes(eventSpikes, trlType);
            trlIDspkDta = sortrows([trlIDs(2,:); trlIDs(1,:); squeeze(evtSpks)]');
            sortedTrlIDs = trlIDspkDta(:,1);
            sortedSpkDta = trlIDspkDta(:,3:end);
            rowIDs = repmat((1:size(sortedSpkDta,1))', [1,size(sortedSpkDta,2)]);
            colIDs = repmat(window(1):window(2), [size(sortedSpkDta,1),1]);
            for pos = 1:max(sortedTrlIDs)
                posRowLog = sortedTrlIDs==pos;
                posSpksLog = sortedSpkDta(posRowLog,:)~=0 & ~isnan(sortedSpkDta(posRowLog,:));
                posRows = rowIDs(posRowLog,:);
                posCols = colIDs(posRowLog,:);
                scatter(posCols(posSpksLog), posRows(posSpksLog),5,...
                    'markeredgecolor', 'none', 'markerfacecolor', obj.PositionColors(pos,:))
                hold on;
            end
            set(gca, 'ydir', 'reverse');            
        end        
        %% Plot Gaussian Firing Rate By Position
        function PlotGaussFireByPos(obj, window, alignment, uni, trlType)
            uniCol = strcmp(obj.ensembleMatrixColIDs, uni);
            eventSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix(:,uniCol), [window(1)-obj.gaussWinDur/2, window(2)+obj.gaussWinDur/2], alignment);
            [evtSpks, trlIDs] = ExtractTrialSpikes(obj, eventSpikes, trlType);
            
            instFRgauss = gausswin(obj.gaussWinDur);
            instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(obj.tsVect)));
            
            instFR = nan(size(squeeze(evtSpks)));
            for trl = 1:size(trlIDs,2)
                instFR(:,trl) = conv(evtSpks(:,:,trl), instFRgauss, 'same');
            end
            instFR = instFR(obj.gaussWinDur/2+1:size(instFR,1)-(obj.gaussWinDur/2),:);
            for pos = 1:max(trlIDs(2,:))
                curPosTrls = instFR(:,trlIDs(2,:)==pos);
                obj.PlotMeanVarLine(window(1):window(2),curPosTrls,2,0.05,obj.PositionColors(pos,:));
            end
        end
    end
end