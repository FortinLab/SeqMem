classdef SingleUnit_SM < SeqMem
    % Single unit analyses from the Sequence Memory project
    properties % General Variables
        sessionSpikes
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
            [iscPRE, trlIDs] = ExtractTrialSpikes(obj, squeeze(preTrialSpikes), 'isc');
            rlyTrialSpikes = obj.BinTrialEventSpikes([0 500], 'PokeIn');
            iscRLY = ExtractTrialSpikes(obj, squeeze(rlyTrialSpikes), 'isc');
            latTrialSpikes = obj.BinTrialEventSpikes([-500 0], 'PokeOut');
            iscLAT = ExtractTrialSpikes(obj, squeeze(latTrialSpikes), 'isc');
            pstTrialSpikes = obj.BinTrialEventSpikes([0 500], 'PokeOut');
            iscPST = ExtractTrialSpikes(obj, squeeze(pstTrialSpikes), 'isc');

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
            rwdSpikes = obj.BinTrialEventSpikes([-250 250], 'RewardSignal');
            [iscRWD, rwdTrls] = ExtractTrialSpikes(obj, squeeze(rwdSpikes), 'corr');
            errSpikes = obj.BinTrialEventSpikes([-250 250], 'ErrorSignal');
            [iscERR, errTrls] = ExtractTrialSpikes(obj, squeeze(errSpikes), 'incorr');
            
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
            [preERR, trialIDs] = ExtractTrialSpikes(obj, squeeze(preErrSigSpikes), 'incorr');
            pstErrSigSpikes = obj.BinTrialEventSpikes([0 250], 'ErrorSignal');
            pstERR = ExtractTrialSpikes(obj, squeeze(pstErrSigSpikes), 'incorr');
            
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
        %% Evaluate Reward Responses
        function [anovaTable, anovaStats, anovaData, groupIDs] = RewardSpiking(obj)
            preRwdSigSpikes = obj.BinTrialEventSpikes([-250 0], 'RewardSignal');
            [preRWD, trialIDs] = ExtractTrialSpikes(obj, squeeze(preRwdSigSpikes), 'corr');
            pstRwdSigSpikes = obj.BinTrialEventSpikes([0 250], 'RewardSignal');
            pstRWD = ExtractTrialSpikes(obj, squeeze(pstRwdSigSpikes), 'corr');
            
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
end