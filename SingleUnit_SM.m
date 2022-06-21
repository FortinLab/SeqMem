classdef SingleUnit_SM < SeqMem
    % Single unit analyses from the Sequence Memory project
    properties % General Variables
        gaussWinDur = 200
        gaussWinWidth = 2.5
        numChancePerms = 100
        binSize = 200
        dsRate = 5
        savePlots = false
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
            elseif strcmp(trlType,'iscSANSA')
                iscLog = isLog & perfLog;
                posLog = posVect~=1;
                if trlDim == 1
                    trlSpikes = spiking(iscLog & posLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,iscLog & posLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,iscLog & posLog);
                end
                trlIDs = ssnLogs(:,iscLog & posLog);
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
            elseif strcmp(trlType, 'corrSANSA')
                posLog = posVect~=1;
                if trlDim == 1
                    trlSpikes = spiking(perfLog & posLog,:,:);
                elseif trlDim == 2
                    trlSpikes = spiking(:,perfLog & posLog,:);
                elseif trlDim == 3
                    trlSpikes = spiking(:,:,perfLog & posLog);
                end
                trlIDs = ssnLogs(:,perfLog & posLog);
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
        %% Bin trial event spiking
        function binnedSpikes = BinTrialEventSpikes(obj,window,alignment)
            trialSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix, window, alignment);
            binnedSpikes = sum(trialSpikes,1);
        end
        %% Sliding bin trial event spiking
        function [binnedSpikes, tsVect] = SlidingBinTrialEventSpikes(obj,window,alignment,binSize,dsRate)
            trialSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix,[window(1)-binSize/2 window(2)+binSize/2], alignment);
            binSpikes = nan(size(trialSpikes));
            for t = 1:size(trialSpikes,3)
                for u = 1:size(trialSpikes,2)
                    binSpikes(:,u,t) = conv(trialSpikes(:,u,t),ones(1,binSize)./(binSize/obj.sampleRate), 'same');
                end
            end
            binnedSpikes = downsample(binSpikes(binSize/2+1:size(binSpikes,1)-(binSize/2),:,:),dsRate);
            tsVect = window(1):dsRate:window(2);
        end        
        %% Sliding window F-Test
        function [fVals] = SlidingWindowFtest(~,spiking,idVect)
            % Written assuming spiking is NxMxO N=time; M=unit; O=trial
            % I'm feeling a little too lazy to bother with matching vect sizes and coding contingencies right now
            fVals = nan(size(spiking,1),size(spiking,2));
            for t = 1:size(spiking,1)
                for u = 1:size(spiking,2)
                    [~,tab] = anova1(squeeze(spiking(t,u,:)), idVect, 'off');
                    if ~isempty(tab{2,5})
                        fVals(t,u) = tab{2,5};
                    else
                        fVals(t,u) = 1;
                    end
                end
            end
        end
         %% Quantify Modulation with F-vals
        function [rawF, tsVect, evtSpks, trlIDs] = SlidingWindowModIC(obj, window, alignment, critVar, trlType)
            [binnedSpikes, tsVect] = obj.SlidingBinTrialEventSpikes(window,alignment,obj.binSize,obj.dsRate);
            [evtSpks, trlIDs] = obj.ExtractTrialSpikes(binnedSpikes, trlType);
            if strcmp(critVar, 'pos')
                varID = trlIDs(2,:);
            elseif strcmp(critVar, 'odr')
                varID = trlIDs(3,:);
            end
            [rawF] = obj.SlidingWindowFtest(evtSpks,varID);
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
        %% Quantify Position Information
        function [posInfo, tsVect, binnedSpikes, trlIDs] = QuantPosInfo(obj,window,alignment,critVar,trlType)
            [posInfo, tsVect, binnedSpikes, trlIDs] = obj.SlidingWindowModIC(window,alignment,critVar,trlType);
        end
        %% Quantify Position Info & Spike Rates
        function [posSpkCorr, posInfo, varSpks, tsVect] = QuantPosSpkCorr(obj,window,alignment,critVar,trlType)
            [posInfo,tsVect] = obj.SlidingWindowModIC(window,alignment,critVar,trlType);
            [binnedSpikes, ~] = obj.SlidingBinTrialEventSpikes(window,alignment,obj.binSize,obj.dsRate);
            [evtSpks, trlIDs] = obj.ExtractTrialSpikes(binnedSpikes, trlType);
            if strcmp(critVar, 'pos')
                varID = trlIDs(2,:);
            elseif strcmp(critVar, 'odr')
                varID = trlIDs(3,:);
            end
            vars = unique(varID);
            if strcmp(trlType,'corrSANSA')
                vars = vars(2:end);
            end
            varSpks = nan(size(evtSpks,1), size(evtSpks,2), length(vars)+1);
            posSpkCorr = nan(size(evtSpks,2), length(vars)+1);
            for v = 1:length(vars)
                varSpks(:,:,v) = mean(evtSpks(:,:,varID==vars(v)),3);
                for u = 1:size(evtSpks,2)
%                     posSpkCorr(u,v) = pdist([posInfo(:,u), varSpks(:,u,v)]', 'cosine');
                    posSpkCorr(u,v) = pdist([posInfo(:,u), varSpks(:,u,v)]', 'correlation');
                end
            end
            varSpks(:,:,length(vars)+1) = mean(evtSpks,3);
            for u = 1:size(evtSpks,2)
%                 posSpkCorr(u,length(vars)+1) = pdist([posInfo(:,u), varSpks(:,u,length(vars)+1)]', 'cosine');
                posSpkCorr(u,length(vars)+1) = pdist([posInfo(:,u), varSpks(:,u,length(vars)+1)]', 'correlation');
            end
        end
        %% Examine Firing Rate Changes Over Time During the Trial
        function [meanSpkDiff,tsVect] = SpikeDiffAcrossTime(obj,window,alignment,trlType,winStepSize,horizon)
            [binnedSpikes, tsVect] = obj.SlidingBinTrialEventSpikes(window,alignment,winStepSize,winStepSize);
            tempTrlSpks = obj.ExtractTrialSpikes(binnedSpikes,trlType);
            meanSpkDiff = nan(length(tsVect)-horizon,horizon, size(tempTrlSpks,2));
            for u = 1:size(tempTrlSpks,2)
                tempSpkDiff = nan(length(tsVect)-horizon,horizon, size(tempTrlSpks,3));
                for trl = 1:size(tempTrlSpks,3)
                    for t = 1:size(tempTrlSpks,1)-horizon
                        for tt = 1:horizon
                            tempSpkDiff(t,tt,trl) = tempTrlSpks(t,u,trl)-tempTrlSpks(t+tt,u,trl);
                        end
                    end
                end
                tempSpkDiff(isinf(tempSpkDiff)) = nan;
                meanSpkDiff(:,:,u) = mean(tempSpkDiff,3,'omitnan');
%                 tempMeanSpks = mean(squeeze(tempTrlSpks(:,u,:)),2);
%                 for t = 1:length(tempMeanSpks)-horizon
%                     for tt = 1:horizon
%                         meanSpkDiff(t,tt,u) = (tempMeanSpks(t)-tempMeanSpks(t+tt))/tempMeanSpks(t);
%                     end
%                 end
            end
        end
        %% Examine Ensemble Similarity Over Time During the Trial
        function [meanStateSim, tsVect] = EnsembleStateDiffAcrossTime(obj,window,alignment,trlType,winStepSize,horizon)
            [binnedSpikes, tsVect] = obj.SlidingBinTrialEventSpikes(window,alignment,winStepSize,winStepSize);
            tempTrlSpks = obj.ExtractTrialSpikes(binnedSpikes,trlType);
            tempStateDiff = nan(length(tsVect)-horizon,horizon, size(tempTrlSpks,3));
            for trl = 1:size(tempTrlSpks,3)
                for t = 1:size(tempTrlSpks,1)-horizon
                    for tt = 1:horizon
                        tempStateDiff(t,tt,trl) = pdist([tempTrlSpks(t,:,trl);tempTrlSpks(t+tt,:,trl)],'euclidean');
                    end
                end
            end
            meanStateSim = mean(tempStateDiff,3,'omitnan');
        end
        %% Examine Firing Rate Changes Relative to Peak Rate
        function [changeVals,ndxs,maxs,windowTSvect] = SpikeDiffWindowedOnMaxRate(obj,window,alignment,trlType,winStepSize,horizon)
            binnedSpikes= obj.SlidingBinTrialEventSpikes(window,alignment,winStepSize,winStepSize);
            tempTrlSpks = obj.ExtractTrialSpikes(binnedSpikes,trlType);
            changeVals = nan(size(tempTrlSpks,2),horizon*2+1);
            ndxs = nan(size(tempTrlSpks,2),1);
            maxs = nan(size(tempTrlSpks,2),1);
            for u = 1:size(tempTrlSpks,2)
                tempChangeVals = nan(size(tempTrlSpks,3),horizon*2+1);
                tempMaxNdxs = nan(size(tempTrlSpks,3),1);
                tempMaxSpk = nan(size(tempTrlSpks,3),1);
                for trl = 1:size(tempTrlSpks,3)
                    tempChangeVals(trl,horizon+1) = 0;
                    curMaxNdx = find(tempTrlSpks(:,u,trl)==max(tempTrlSpks(:,u,trl)),1,'first');
                    for t = 1:horizon
                        if curMaxNdx+t<size(tempTrlSpks,1)
                            tempChangeVals(trl,horizon+t+1) = tempTrlSpks(curMaxNdx+t,u,trl)-tempTrlSpks(curMaxNdx,u,trl);
                        end
                        if curMaxNdx-t>0
                            tempChangeVals(trl,horizon-t+1) = tempTrlSpks(curMaxNdx-t,u,trl)-tempTrlSpks(curMaxNdx,u,trl);
                        end
                    end
                    tempMaxNdxs(trl) = curMaxNdx;
                    tempMaxSpk(trl) = tempTrlSpks(curMaxNdx,u,trl);
                end
                changeVals(u,:) = mean(tempChangeVals,1,'omitnan');
                ndxs(u) = mean(tempMaxNdxs);
                maxs(u) = mean(tempMaxSpk);
            end
            windowTSvect = winStepSize*-horizon:winStepSize:winStepSize*horizon;
        end
        %% Examine Structure of Ordinal Information at Peak Info
        function [unitPosInfoPeaksInfo] = InterrogatePeakInfo(obj,window,alignment,critVar,trlType,critF)
            [posInfo, tsVect, evtSpks, trlIDs] = obj.SlidingWindowModIC(window,alignment,critVar,trlType);
            if strcmp(critVar, 'pos')
                varID = trlIDs(2,:);
            elseif strcmp(critVar, 'odr')
                varID = trlIDs(3,:);
            end
            unqIDs = unique(varID);
            unitPosInfoPeaksInfo = struct('Unit', obj.ensembleMatrixColIDs,...
                'Peaks', cell(1,size(posInfo,2)),...
                'Locations', cell(1,size(posInfo,2)),...
                'Widths', cell(1,size(posInfo,2)),...
                'Prominance', cell(1,size(posInfo,2)),...
                'Means', cell(1,size(posInfo,2)),...
                'LinearFit', cell(1,size(posInfo,2)),...
                'PolyRatio', cell(1,size(posInfo,2)));
            for u = 1:size(posInfo,2)
%                 [pks,locs,widths,proms] = findpeaks(posInfo(:,u), 'MinPeakProminence', critF, 'MinPeakDistance', obj.binSize/4, 'MinPeakHeight', critF);
                [pks,locs,widths,proms] = findpeaks(smooth(posInfo(:,u)), 'MinPeakDistance', obj.binSize/4, 'MinPeakHeight', critF);
%                 [pks,locs,widths,proms] = findpeaks(posInfo(:,u), 'MinPeakProminence', critF, 'MinPeakHeight', critF);                
%                 [pks,locs,widths,proms] = findpeaks(posInfo(:,u), 'MinPeakHeight', critF);
                peakRvals = nan(size(pks));
                tempMeans = nan(length(pks), length(unqIDs));
                polyErrorRat = nan(size(pks));
                for p = 1:length(pks)
                    tempFR = squeeze(evtSpks(locs(p),u,:));
                    peakRvals(p) = corr(varID', tempFR, 'rows', 'pairwise');
                    p1 = polyfit(varID',tempFR,1);
                    x1 = polyval(p1,unqIDs);
                    p2 = polyfit(varID',tempFR,2);
                    x2 = polyval(p2,unqIDs);
                    sse1 = cell(length(unqIDs),1);
                    sse2 = cell(length(unqIDs),1);
                    for ui = 1:length(unqIDs)
                        idFR = tempFR(varID==unqIDs(ui));
                        tempMeans(p,ui) = mean(idFR,'omitnan');
                        sse1{ui} = (idFR-x1(ui)).^2;
                        sse2{ui} = (idFR-x2(ui)).^2;
                    end
                    polyErrorRat(p) = sum(cell2mat(sse1))/sum(cell2mat(sse2));
                end
                unitPosInfoPeaksInfo(u).Peaks = pks;
                unitPosInfoPeaksInfo(u).Locations = tsVect(locs);
                unitPosInfoPeaksInfo(u).Widths = widths;
                unitPosInfoPeaksInfo(u).Prominance = proms;
                unitPosInfoPeaksInfo(u).Means = tempMeans;
                unitPosInfoPeaksInfo(u).LinearFit = peakRvals;
                unitPosInfoPeaksInfo(u).PolyRatio = polyErrorRat;
            end
        end
    end
    %% Visualizations
    methods
        %% Plot Unit Summary
        function PlotUnitSummary(obj,window,alignment,uni)
            if nargin == 1
                uni = obj.ensembleMatrixColIDs;
                window = [-1200 2000];
                alignment = 'PokeIn';
            else
                uni = {uni};
            end
            aniSsn = strsplit(obj.pathDir(regexp(obj.pathDir, '([a-z]*|[A-Z]*[0-9]*)_Session[0-9*]'):end), '_');
            [posInfo, tsVect] = obj.QuantPosInfo(window,alignment,'pos','corr');
            posInfoSANSA = obj.QuantPosInfo(window,alignment,'pos','corrSANSA');
            for u = 1:length(uni)
                curUniInfo = obj.unitInfo(strcmp({obj.unitInfo.UnitName}, uni{u}));
                wfVect = 1:length(curUniInfo.TemplateMean{1});
                figure; 
                annotation(gcf,'textbox', 'Units', 'normalized', 'Position', [0.05 0.9 0.5 0.1],...
                    'String', sprintf('%s %s %s',uni{u}, aniSsn{1}, aniSsn{2}),...
                    'FontSize',15, 'FontWeight', 'bold', 'EdgeColor', 'none',...
                    'horizontalalignment', 'left', 'interpreter', 'none');
                wfSP(1) = axes(gcf, 'Units', 'Normalized', 'Position', [0.05 0.8 0.08 0.15]);
                plot(wfVect, curUniInfo.TemplateMean{1}, 'k');
                hold on;
                patch('XData', [wfVect(:); flipud(wfVect(:))],...
                    'YData', [(curUniInfo.TemplateMean{1}+curUniInfo.TemplateStDev{1})'; flipud((curUniInfo.TemplateMean{1}-curUniInfo.TemplateStDev{1})')],...
                    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
                axis tight;
                wfSP(2) = axes(gcf, 'Units', 'Normalized', 'Position', [0.15 0.8 0.08 0.15]);
                plot(wfVect, curUniInfo.TemplateMean{2}, 'k');
                hold on;
                patch('XData', [wfVect(:); flipud(wfVect(:))],...
                    'YData', [(curUniInfo.TemplateMean{2}+curUniInfo.TemplateStDev{2})'; flipud((curUniInfo.TemplateMean{2}-curUniInfo.TemplateStDev{2})')],...
                    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
                axis tight;
                wfSP(3) = axes(gcf, 'Units', 'Normalized', 'Position', [0.25 0.8 0.08 0.15]);
                plot(wfVect, curUniInfo.TemplateMean{3}, 'k');
                hold on;
                patch('XData', [wfVect(:); flipud(wfVect(:))],...
                    'YData', [(curUniInfo.TemplateMean{3}+curUniInfo.TemplateStDev{3})'; flipud((curUniInfo.TemplateMean{3}-curUniInfo.TemplateStDev{3})')],...
                    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
                axis tight;
                wfSP(4) = axes(gcf, 'Units', 'Normalized', 'Position', [0.35 0.8 0.08 0.15]);
                plot(wfVect, curUniInfo.TemplateMean{4}, 'k');
                hold on;
                patch('XData', [wfVect(:); flipud(wfVect(:))],...
                    'YData', [(curUniInfo.TemplateMean{4}+curUniInfo.TemplateStDev{4})'; flipud((curUniInfo.TemplateMean{4}-curUniInfo.TemplateStDev{4})')],...
                    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
                axis tight;
                linkaxes(wfSP,'xy');
                
                autocorrSP = axes(gcf, 'Units', 'Normalized', 'Position', [0.06 0.5 0.15 0.25]);
                isi = diff(find(obj.ensembleMatrix(:,strcmp(obj.ensembleMatrixColIDs, uni{u}))==1));
                histogram([isi;isi*-1],-200:1:200);
                xlabel('ITI (ms)');
                
                %%%%%%%%%%%%%%%%%% Spike phase plots would go in here whenever they get coded up...

                sp(1) = axes(gcf, 'Units', 'Normalized', 'Position', [0.45 0.65 0.5 0.3]);
                obj.PlotRastersByPos(window,alignment,uni{u},'isc');
                axis tight;
                axis off;
                sp(2) = axes(gcf, 'Units', 'Normalized', 'Position', [0.45 0.5 0.5 0.125]);
                uniLog = strcmp(obj.ensembleMatrixColIDs, uni{u});
                yyaxis(sp(2),'left');
                plot(tsVect, posInfo(:,uniLog), 'k');
                yyaxis(sp(2),'right');
                plot(tsVect, posInfoSANSA(:,uniLog), 'r');
                set(gca, 'xticklabel', []);
                axis tight;
                sp(3) = axes(gcf, 'Units', 'Normalized', 'Position', [0.45 0.05 0.5 0.425]);
                obj.PlotGaussFireByPos(window,alignment,uni{u},'isc');
                axis tight;
                linkaxes(sp, 'x');
                drawnow;
                if obj.savePlots                    
                    orient(gcf, 'tall');
                    orient(gcf, 'landscape');
                    print('-painters', gcf, '-dpdf', sprintf('%s_%s_%s',aniSsn{1}, aniSsn{2}, uni{u}));
%                     saveas(gcf, sprintf('%s_%s_%s.png',aniSsn{1}, aniSsn{2}, uni{u}), 'png');
                end
            end
        end
        %% Plot Trial Summary
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
            plot([0 0], get(gca, 'ylim'), '-k');
            scatter([obj.trialInfo(trlIDspkDta(:,2)).PokeDuration]*1000, 1:size(trlIDspkDta,1), 20, '|k', 'linewidth', 1);
            set(gca, 'ydir', 'reverse');            
        end        
        %% Plot Gaussian Firing Rate By Position
        function PlotGaussFireByPos(obj, window, alignment, uni, trlType)
            uniCol = strcmp(obj.ensembleMatrixColIDs, uni);
            eventSpikes = obj.ExtractTrialMatrix(obj.ensembleMatrix(:,uniCol), [window(1)-obj.gaussWinDur/2, window(2)+obj.gaussWinDur/2], alignment);
            [evtSpks, trlIDs] = obj.ExtractTrialSpikes(eventSpikes, trlType);
            
            instFRgauss = gausswin(obj.gaussWinDur,obj.gaussWinWidth);
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
            plot([0 0], get(gca, 'ylim'), '-k');
            plot(repmat(median([obj.trialInfo(trlIDs(1,:)).PokeDuration]), [1,2])*1000, get(gca, 'ylim'), '-k');
            plot(repmat(median([obj.trialInfo(trlIDs(1,:)).RewardIndex]-[obj.trialInfo(trlIDs(1,:)).PokeInIndex]), [1,2]), get(gca, 'ylim'), '--k');
        end
    end
end