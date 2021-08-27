classdef PFC_TrialPhaseTime_MLB_SM < MLB_SM
    properties
        alignment = 'PokeIn'
        trialWindow = [-500 1200]
        freqWin = [4 12]
        radPrSamp
        phaseBins = deg2rad(0:10:359)
        binSizeDeg = 40
        dsRateDeg = 1
        cycs2pad = 2
        numCycs
        numCycsPre
    end
    properties % Trial Organized Data
        lfpCyc
        lfpPowerTrialTimeMatrix
        spikeTrialTimeMatrix        
    end    
    properties
        unwrapIdealPhaseBins        
        unwrapIdealPhaseLog
        unwrapOdorLog
    end
    properties
        fiscLikes
        fiscPosts
        fiscOdorDecode
        fiscCycleDecode
        fiscPhaseDecode
    end
    %% Methods
    methods
        %% Object Creation
        %#ok<*EXIST>
        function obj = PFC_TrialPhaseTime_MLB_SM(fileDir, dsRate, binSize, trialWindow, freqWin, phaseBins, alignment)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@MLB_SM(fileDir);
            if nargin >= 2
                if exist('dsRate')==1
                    obj.dsRate = dsRate;
                end
                if exist('binSize')==1
                    obj.binSize = binSize;
                end
                if exist('trialWindow')==1
                    obj.trialWindow = trialWindow;
                end
                if exist('freqWin')==1
                    obj.freqWin = freqWin;
                end                
                if exist('phaseBins')==1
                    obj.phaseBins = phaseBins;
                end
                if exist('alignment')==1
                    obj.alignment = alignment;
                end
                obj.RunAnalysis;
            else
            end
            
        end
    end
    
    methods
        %% Run the full analysis
        function RunAnalysis(obj)
            fprintf('Analyzing\n');
            obj.PreProcessLFPandSpikes;
            obj.CompileFISlikes;
            obj.DecodeFIS_L1O;
        end
        %% Plot FIS Odor Decode
        function PlotFIS_OdorDecode(obj)
            figure;
            for o = 1:4
                plot(smooth(nanmean(obj.fiscOdorDecode==o,2))); 
                hold on; 
            end
            boundaries = find(diff(obj.unwrapOdorLog));
            for b = 1:length(boundaries)
                plot([boundaries(b) boundaries(b)], get(gca, 'ylim'), '--k');
            end
            drawnow;
            annotation(gcf,'textbox', [0 0.95 1 0.05],'String', obj.behavMatFile,...
                'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
        end

        %%
        function DecodeFIS_L1O(obj)
            obj.fiscPosts = nan(length(obj.unwrapIdealPhaseLog), length(obj.unwrapIdealPhaseLog), size(obj.fiscTrials,2));
            for trl = 1:size(obj.fiscTrials,2)
                tempLikes = obj.fiscLikes;
                tempObsv = obj.fiscLikes(:,:,trl);
                tempLikes(:,:,trl) = [];
                obj.fiscPosts(:,:,trl) = obj.CalcStaticBayesPost(mean(tempLikes,3), tempObsv);
            end
            [obj.fiscOdorDecode, ~] = obj.DecodeBayesPost(obj.fiscPosts, obj.unwrapOdorLog);
        end
        %%
        function CompileFISlikes(obj)
            if isempty(obj.spikeTrialTimeMatrix)
                obj.PreProcessLFPandSpikes;
            end
            obj.fiscLikes = nan(length(obj.unwrapIdealPhaseLog),size(obj.spikeTrialTimeMatrix,2), size(obj.fiscTrials,2));
            for seq = 1:size(obj.fiscTrials,2)
                tempMat = [];
                for trl = 1:size(obj.fiscTrials,1)
                    tempMat = [tempMat; obj.spikeTrialTimeMatrix(:,:,obj.fiscTrials(trl,seq))];
                end
                obj.fiscLikes(:,:,seq) = tempMat;
            end
        end
        %%
        function PreProcessLFPandSpikes(obj)
            %#ok<*AGROW>
            %% First extract the trial data
            if isempty(obj.lfpMatrix)
                obj.CompileLFPmatrix;
            end
            tetWMU = find(obj.numUniPerTet==max(obj.numUniPerTet), 1, 'first');
            % To remove the pre/post session recordings from merged files that were cut to have unit IDs across recordings (for LFP analysis only)
            sessionWindow = [obj.trialInfo(1).PokeInIndex-1000; obj.trialInfo(end).PokeOutIndex+2000];
           
            [~, lfpPhase, lfpPwr] = obj.SimpleFilter(obj.lfpMatrix(sessionWindow(1):sessionWindow(2),tetWMU), obj.freqWin);
            obj.radPrSamp = mean(diff(unwrap(lfpPhase)));
            % Unless a binSize is defined a priori, calculate binsize based on defined degrees to make the gaussian based in degree
            % (e.g. 45degrees by default)
            if isempty(obj.binSize)
                obj.binSize = floor(deg2rad(obj.binSizeDeg)/mean(diff(obj.phaseBins)));
            end
            % Unless a setpsize is defined a priori, calculate dsRate based on defined degrees to make the steps based in degrees
            % (default 1 degree since we're not hurting for space)
            if isempty(obj.dsRate)
                obj.dsRate = ceil(deg2rad(obj.dsRateDeg)/mean(diff(obj.phaseBins)));
            end
            lfpPhase = [zeros(sessionWindow(1)-1,1); lfpPhase + pi; zeros(size(obj.lfpMatrix,1)-sessionWindow(2),1)]; % Adding pi to lfpPhase here shifts the cycles from -pi:pi to 0:2*pi and makes the trough = 0. Makes it a bit easier to handle with the unwrapping since nothing is now negative
            lfpPwr = [zeros(sessionWindow(1)-1,1); lfpPwr; zeros(size(obj.lfpMatrix,1)-sessionWindow(2),1)];

            % Pad the window for extraction by the equivalent of ~2 cycles for the sake of smoothing/binning everything
            paddedWindow = round([obj.trialWindow(1)-(1000/(mean(obj.freqWin)))*obj.cycs2pad obj.trialWindow(2)+(1000/(mean(obj.freqWin)))*obj.cycs2pad]);
            trialTime = obj.ExtractTrialMatrix(obj.tsVect, paddedWindow, obj.alignment);
            zroTime = obj.ExtractTrialMatrix(obj.tsVect, [0 0], obj.alignment);
            trialTimeLog = trialTime(:,:,1)-zroTime(:,:,1);
            tempLFPphaseTrialTimeMatrix = obj.ExtractTrialMatrix(lfpPhase, paddedWindow, obj.alignment);
            tempLFPpowerTrialTimeMatrix = obj.ExtractTrialMatrix(lfpPwr, paddedWindow, obj.alignment);
            tempSpikeTrialTimeMatrix = obj.ExtractTrialMatrix(obj.ensembleMatrix, paddedWindow, obj.alignment);
            
            if isempty(obj.numCycs) % if this value isn't pre-defined when the analysis is initiated
                % Use InSeq correct trials to identify the number of cycles
                % that at least 75% of inseq correct trials have
                iscLog = ([obj.trialInfo.Odor] == [obj.trialInfo.Position]) &...
                    ([obj.trialInfo.Performance] == 1);
                iscPhase = tempLFPphaseTrialTimeMatrix(:,:,iscLog);
                tempNumCycs = nan(1,sum(iscLog));
                for trl = 1:sum(iscLog)
                    [~, phaseBounds] = findpeaks(diff(iscPhase(:,:,trl))*-1, 'MinPeakProminence', 3);
                    tempNumCycs(trl) = length(phaseBounds);
                end
                [counts, bins] = histcounts(tempNumCycs, mean(obj.freqWin)-20:mean(obj.freqWin)+20);
                csHC = cumsum(counts)./sum(counts);
                obj.numCycs = min(bins(csHC>=0.25));
            end
            if isempty(obj.numCycsPre)
                obj.numCycsPre = round(abs((obj.trialWindow(1)-0)/(1000/mean(obj.freqWin)))) + obj.cycs2pad;
            end
            % Create idealized phase vectors for binning and analysis purposes
            obj.unwrapIdealPhaseBins = unwrap(repmat(obj.phaseBins', [obj.numCycs,1])) - (2*pi*obj.numCycsPre);
            tempTrialCycPower = nan(obj.numCycs, length(obj.trialInfo));
            tempPhaseSpikeBins = nan(length(obj.unwrapIdealPhaseBins)-1, size(obj.ensembleMatrix,2), length(obj.trialInfo));
            % Now bin the spiking data on each trial and determine average power per cycle
            for trl = 1:length(obj.trialInfo)
                %%%% Could potentially align phases relative to trial start
                %%%% here
                tempTrlPhase = tempLFPphaseTrialTimeMatrix(:,:,trl);
                [~, phaseBounds] = findpeaks(diff(tempTrlPhase)*-1, 'MinPeakProminence', 3);
                phaseBoundsLog = false(size(tempTrlPhase,1),1);
                phaseBoundsLog(phaseBounds) = true;
                startFirstCyclePostPokeIn = find(phaseBoundsLog & trialTimeLog>=0, 1, 'first');
                tempLFPuwPhase = unwrap(tempTrlPhase);
                tempLFPuwPhase = tempLFPuwPhase - tempLFPuwPhase(startFirstCyclePostPokeIn);
                for u = 1:size(obj.ensembleMatrix,2)
                    tempCounts = histcounts(tempLFPuwPhase(tempSpikeTrialTimeMatrix(:,u,trl)>=1), obj.unwrapIdealPhaseBins);
%                     binnedMtx(:,u,t) = conv(tempMtx(:,u,t),ones(1,obj.binSize)./(obj.binSize/obj.sampleRate), 'same');
                    tempPhaseSpikeBins(:,u,trl) = conv(tempCounts, ones(1,obj.binSize)./(obj.binSize*obj.radPrSamp), 'same');
                end
                tempCycPwr = nan(length(phaseBounds),1);
                tempPhaseBounds = [1; phaseBounds];
                for cyc = 2:length(tempPhaseBounds)
                    tempCycPwr(cyc-1) = mean(tempLFPpowerTrialTimeMatrix(tempPhaseBounds(cyc-1):tempPhaseBounds(cyc),:,trl));
                end
                firstCycNum = find(phaseBounds==startFirstCyclePostPokeIn)+1;
                if firstCycNum - obj.numCycsPre <= 0
                    tempAlignedCycs = [nan(firstCycNum-obj.numCycsPre+1,1);...
                    tempCycPwr(1:firstCycNum-1+((size(tempTrialCycPower,1) - obj.numCycsPre) - ((size(tempTrialCycPower,1) - obj.numCycsPre) - (length(tempCycPwr)-firstCycNum))))];
                else
                    tempAlignedCycs = tempCycPwr(firstCycNum-obj.numCycsPre:firstCycNum-1+((size(tempTrialCycPower,1) - obj.numCycsPre) - ((size(tempTrialCycPower,1) - obj.numCycsPre) - (length(tempCycPwr)-firstCycNum))));
                end
                if length(tempAlignedCycs) > size(tempTrialCycPower,1)
                    tempAlignedCycs = tempAlignedCycs(1:size(tempTrialCycPower,1));
                end
                tempTrialCycPower(1:length(tempAlignedCycs),trl) = tempAlignedCycs;
            end
            radsUnpadLog = [false(length(obj.phaseBins)*obj.cycs2pad,1);...
                true(length(obj.phaseBins)*(obj.numCycs-(obj.cycs2pad*2)),1);...
                false((length(obj.phaseBins)*obj.cycs2pad)-1,1)];
            tempUnwrapPhaseLog = obj.unwrapIdealPhaseBins(1:end-1)-mean(diff(obj.unwrapIdealPhaseBins))/2;
            obj.unwrapIdealPhaseLog = repmat(tempUnwrapPhaseLog(radsUnpadLog), [4 1]);
            obj.unwrapOdorLog = [ones(sum(radsUnpadLog),1); ones(sum(radsUnpadLog),1)*2; ones(sum(radsUnpadLog),1)*3; ones(sum(radsUnpadLog),1)*4];
            obj.spikeTrialTimeMatrix = tempPhaseSpikeBins(radsUnpadLog,:,:);
            cycUnpadLog = [false(obj.cycs2pad,1); true(obj.numCycs-(obj.cycs2pad*2),1); false(obj.cycs2pad,1)];
            obj.lfpPowerTrialTimeMatrix = tempTrialCycPower(cycUnpadLog,:);
        end
    end
    methods % MLB Analysis Methods
        %%%%%% NEED TO UPDATE THIS TO RUN ON SPK/RAD RATHER THAN SPK/SEC
        %% Calculate MLB
        function post = CalcStaticBayesPost(obj,likely, obsv)
            post = nan(size(obsv,1), size(likely,1), size(obsv,3));
            for trl = 1:size(obsv,3)
                for t = 1:size(obsv,1)
                    p = nan(size(likely));
                    curPopVect = round(obsv(t,:,trl)*(obj.binSize*obj.radPrSamp)); % use round here to get rid of floating point bs
                    curPopFact = factorial(curPopVect);
                    for u = 1:size(likely,2)
                        curAvgUniFR = likely(:,u);
                        p(:,u) = (((obj.binSize*obj.radPrSamp).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
                    end
                    pp = prod(p,2);
                    ee = exp(-((obj.binSize*obj.radPrSamp)*sum(likely,2)));
                    tempPost = pp.*ee;
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end
        end
        %% Decode MLB
        function [decode, maxPost] = DecodeBayesPost(~, post, id)
            % Assumes post is in the structure of ObservTime X LikelyTime X Trial
            decode = nan(size(post,1),size(post,3));
            maxPost = nan(size(post,1),size(post,3));
            for o = 1:size(post,3)
                for d = 1:size(post,1)
                    if ~isnan(post(d,1,o))
                        maxPost(d,o) = max(post(d,:,o));
                        decode(d,o) = id(find(post(d,:,o)==maxPost(d,o),1,'first'));
                    end
                end
            end
        end
    end
end