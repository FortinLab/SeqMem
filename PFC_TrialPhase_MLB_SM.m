classdef PFC_TrialPhase_MLB_SM < MLB_SM
    properties % Parameters
        beginTrialAlignment = 'PokeIn'
        endTrialAlignment = 'PokeOut'
        beginTrialWindow = [-500 500]
        endTrialWindow = [-500 500]
        phaseBins = [-pi/4, pi/4;...
            pi/4, pi*3/4;...
            pi*3/4 -pi*3/4;...
            -pi*3/4 -pi/4]
        phaseBinNames = [{'Peak'}, {'Fall'}, {'Trough'}, {'Rise'}]
        oscBandLims = [4 12]
        coordLFPchanID
    end
    properties % Data Structures
        beginTrialSpkPhaseMtx
        beginTrialLFPsig
        beginTrialLFPphase
        beginTrialLFPphaseID
        beginTrialLFPpower
        
        endTrialSpkPhaseMtx
        endTrialLFPsig
        endTrialLFPphase
        endTrialLFPphaseID
        endTrialLFPpower
        
        fisSeqSpikePhaseMtx
        fisSeqLFPphaseIDmtx
    end
    properties % Data Organization Vectors
        beginTrialTime
        endTrialTime
    end
    properties % Posterior Decoding Logicals
        fisSeqSpikeTimeLog
        fisSeqSpikeOdorLog
        fisSeqSpikePhaseLog
        fisSeqTrialPeriodLog
    end
    properties % Trial Logicals
        beginTrialPhaseLog
        endTrialPhaseLog
        beginTrialPeriodLog
        endTrialPeriodLog
    end
    properties % Posterior Distributions
        fisSeqPosts
    end
    properties % Decodings
        fisL1OdecodeOdrTime
        fisL1OdecodeTimeTime
        fisL1OdecodePhaseTime
        
        fisL1OphaseDecodePhase_TrlPrd
        fisL1OphaseDecodePhaseOdor_TrlPrd
        
        fisL1OdecodeOdrTime_PhaseDecode
    end
    %%
    methods
        %% Object Creation
        %#ok<*EXIST>
        function obj = PFC_TrialPhase_MLB_SM(fileDir, binSize, dsRate, oscBandLims, beginWindow, endWindow, phaseBins)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@MLB_SM(fileDir);
            if nargin >= 3
%                 if exist('binSize')==1
                    obj.binSize = binSize;
%                 end
%                 if exist('dsRate')==1
                    obj.dsRate = dsRate;
%                 else
%                     obj.dsRate = 1;
%                 end
                if exist('oscBandLims')==1
                    obj.oscBandLims = oscBandLims;
                end                    
                if exist('beginWindow')==1
                    obj.beginTrialWindow = beginWindow;
                end
                if exist('endWindow')==1
                    obj.endTrialWindow = endWindow;
                end
                if exist('phaseBins')==1
                    obj.phaseBins = phaseBins;
                end
                obj.RunAnalysis;
            else
            end
            
        end
        %% Run the full analysis
        function RunAnalysis(obj)
            fprintf('Analysis Beginning\n');
            obj.IdentifyCoordinatingLFP;
            obj.CompileMLBmtx(obj.beginTrialAlignment);
            obj.CompileMLBmtx(obj.endTrialAlignment);
            obj.CompileFISlikes;
            obj.DecodeFIS_L1O;
        end
        %% Decode FIS trials via leave 1 out
        function DecodeFIS_L1O(obj)
            if isempty(obj.fisSeqSpikePhaseMtx)
                obj.CompileFISlikes;
            end
            fprintf('Calculating MLB\n');
            obj.fisSeqPosts = nan(size(obj.fisSeqSpikePhaseMtx,1), size(obj.fisSeqSpikePhaseMtx,1), size(obj.fisSeqSpikePhaseMtx,3));
            tic;
            for s = 1:size(obj.fisSeqSpikePhaseMtx,3)
                tempISS = obj.fisSeqSpikePhaseMtx;
                tempISS(:,:,s) = [];
                obj.fisSeqPosts(:,:,s) = obj.CalcStaticBayesPost(mean(tempISS,3), obj.fisSeqSpikePhaseMtx(:,:,s));
            end
            toc;
%             for n = 1:4; for m = 1:4;figure; imagesc(mean(obj.fisSeqPosts(obj.fisSeqSpikePhaseLog==n,obj.fisSeqSpikePhaseLog==m,:),3), [0 0.001]); title(sprintf('phase1 = %s; phase2 = %s', obj.phaseBinNames{n}, obj.phaseBinNames{m})); end; end
            %   Decode Odor across trial time
            [tempDecodeOdrTm, ~] = obj.DecodeBayesPost(obj.fisSeqPosts, obj.fisSeqSpikeOdorLog);
            obj.fisL1OdecodeOdrTime = nan(size(tempDecodeOdrTm,1)/(size(obj.phaseBins,1)),4,size(obj.phaseBins,1));
            for phase = 1:size(obj.phaseBins,1)
                tempPhaseDecodeOdor = tempDecodeOdrTm(obj.fisSeqSpikePhaseLog==phase,:);
                for o = 1:4
                    obj.fisL1OdecodeOdrTime(:,o,phase) = sum(tempPhaseDecodeOdor==o,2)./sum(~isnan(tempPhaseDecodeOdor),2);
                end
            end
            obj.fisL1OdecodeOdrPhaseTrlPrd = repmat({nan(4,4,size(obj.phaseBins,1))}, [1,4]);
            %   Decode Time across trial time
            [tempDecodeTime, ~] = obj.DecodeBayesPost(obj.fisSeqPosts, obj.fisSeqSpikeTimeLog);
            for s = 1:size(tempDecodeTime,2)
                tempDecodeTime(:,s) = tempDecodeTime(:,s) - obj.fisSeqSpikeTimeLog;
            end
            obj.fisL1OdecodeTimeTime = nan(size(tempDecodeOdrTm,1)/size(obj.phaseBins,1),size(tempDecodeTime,2),size(obj.phaseBins,1));
            for phase = 1:size(obj.phaseBins,1)
                obj.fisL1OdecodeTimeTime(:,:,phase) = tempDecodeTime(obj.fisSeqSpikePhaseLog==phase,:);
            end
            %   Decode Phase across trial time
            [tempDecodePhase, ~] = obj.DecodeBayesPost(obj.fisSeqPosts, obj.fisSeqSpikePhaseLog);
            obj.fisL1OdecodePhaseTime = nan(size(tempDecodePhase,1)/size(obj.phaseBins,1),size(obj.phaseBins,1),size(obj.phaseBins,1));
            for phase = 1:size(obj.phaseBins,1)
                tempPhaseDecodePhase = tempDecodePhase(obj.fisSeqSpikePhaseLog==phase,:);
                for p = 1:size(obj.phaseBins,1)
                    obj.fisL1OdecodePhaseTime(:,p,phase) = sum(tempPhaseDecodePhase==p,2)./sum(~isnan(tempPhaseDecodePhase),2);
                end
            end    
            
            obj.fisL1OphaseDecodePhase_TrlPrd = cell(1,4);
            obj.fisL1OphaseDecodePhaseOdor_TrlPrd = cell(1,4);
            for prd = 1:4
                tempPrdLog = obj.fisSeqTrialPeriodLog==prd;
                tempOdrLog = repmat(obj.fisSeqSpikeOdorLog(tempPrdLog), [1, size(tempDecodeOdrTm,2)]);
                tempOdrDecode = tempDecodeOdrTm(tempPrdLog,:);
                tempPhaseDecode = tempDecodePhase(tempPrdLog,:);
                tempTruePhase = reshape(obj.fisSeqLFPphaseIDmtx(tempPrdLog,:,:), [sum(tempPrdLog), size(obj.fisSeqLFPphaseIDmtx,3)]);
                
                tempPhaseDecodePhase = nan(size(obj.phaseBins,1));
                tempPhaseDecodePhaseOdor = repmat({nan(4)}, size(obj.phaseBins,1));
                for phase1 = 1:size(obj.phaseBins,1)
                    truePhaseLog = tempTruePhase==phase1;
                    for phase2 = 1:size(obj.phaseBins,1)
                        decodePhaseLog = tempPhaseDecode==phase2;
                        trueDecodePhaseLog = truePhaseLog & decodePhaseLog;
                        tempPhaseDecodePhase(phase2,phase1) = sum(trueDecodePhaseLog(:))/sum(truePhaseLog(:));
                        for odr1 = 1:4
                            trueOdorLog = tempOdrLog == odr1;
                            for odr2 = 1:4
                                decodeOdorLog = tempOdrDecode == odr2;
                                trueDecodeOdorLog = trueOdorLog & decodeOdorLog;
                                tempPhaseDecodePhaseOdor{phase2,phase1}(odr2,odr1) = sum(trueDecodePhaseLog(:) & trueDecodeOdorLog(:))/sum(trueDecodePhaseLog(:) & trueOdorLog(:));
                            end
                        end
                    end
                end
                obj.fisL1OphaseDecodePhase_TrlPrd{prd} = tempPhaseDeocdePhase;
                obj.fisL1OphaseDecodePhaseOdor_TrlPrd{prd} = tempPhaseDecodePhaseOdor;
            end
        end
        %% Compile Fully InSeq Likelihoods
        function CompileFISlikes(obj)
            if isempty(obj.beginTrialSpkPhaseMtx)
                obj.CompileMLBmtx(obj.beginTrialAlignment);
            elseif isempty(obj.endTrialSpkPhaseMtx)
                obj.CompileMLBmtx(obj.endTrialAlignment);
            end
            obj.fisSeqLFPphaseIDmtx = nan((size(obj.beginTrialLFPphaseID,1) + size(obj.endTrialLFPphaseID,1))*4, 1, size(obj.fiscTrials,2));
            obj.fisSeqSpikePhaseMtx = nan((size(obj.beginTrialSpkPhaseMtx,1) + size(obj.endTrialSpkPhaseMtx,1))*4, length(obj.ensembleMatrixColIDs), size(obj.fiscTrials,2));
            for seq = 1:size(obj.fiscTrials,2)
                obj.fisSeqSpikePhaseMtx(:,:,seq) = [obj.beginTrialSpkPhaseMtx(:,:,obj.fiscTrials(1,seq)); obj.endTrialSpkPhaseMtx(:,:,obj.fiscTrials(1,seq));...
                    obj.beginTrialSpkPhaseMtx(:,:,obj.fiscTrials(2,seq)); obj.endTrialSpkPhaseMtx(:,:,obj.fiscTrials(2,seq));...
                    obj.beginTrialSpkPhaseMtx(:,:,obj.fiscTrials(3,seq)); obj.endTrialSpkPhaseMtx(:,:,obj.fiscTrials(3,seq));...
                    obj.beginTrialSpkPhaseMtx(:,:,obj.fiscTrials(4,seq)); obj.endTrialSpkPhaseMtx(:,:,obj.fiscTrials(4,seq))];
                obj.fisSeqLFPphaseIDmtx(:,:,seq) = [obj.beginTrialLFPphaseID(:,:,obj.fiscTrials(1,seq)); obj.endTrialLFPphaseID(:,:,obj.fiscTrials(1,seq));...
                    obj.beginTrialLFPphaseID(:,:,obj.fiscTrials(2,seq)); obj.endTrialLFPphaseID(:,:,obj.fiscTrials(2,seq));...
                    obj.beginTrialLFPphaseID(:,:,obj.fiscTrials(3,seq)); obj.endTrialLFPphaseID(:,:,obj.fiscTrials(3,seq));...
                    obj.beginTrialLFPphaseID(:,:,obj.fiscTrials(4,seq)); obj.endTrialLFPphaseID(:,:,obj.fiscTrials(4,seq))];
            end
            obj.fisSeqSpikeOdorLog = [ones(size(obj.beginTrialSpkPhaseMtx,1)+size(obj.endTrialSpkPhaseMtx,1),1)*1;...
                ones(size(obj.beginTrialSpkPhaseMtx,1)+size(obj.endTrialSpkPhaseMtx,1),1)*2;...
                ones(size(obj.beginTrialSpkPhaseMtx,1)+size(obj.endTrialSpkPhaseMtx,1),1)*3;...
                ones(size(obj.beginTrialSpkPhaseMtx,1)+size(obj.endTrialSpkPhaseMtx,1),1)*4];
            obj.fisSeqSpikeTimeLog = repmat([obj.beginTrialTime; obj.endTrialTime+1+obj.dsRate/obj.sampleRate], [4,1]);
            obj.fisSeqSpikePhaseLog = repmat([obj.beginTrialPhaseLog; obj.endTrialPhaseLog], [4,1]);
            obj.fisSeqTrialPeriodLog = repmat([obj.beginTrialPeriodLog; obj.endTrialPeriodLog], [4, 1]);
        end
        %% Extract Matrix
        function CompileMLBmtx(obj, alignment)
            %% First check the data to make sure everything's here
            if isempty(obj.binSize)
                error('Define binSize');
            end
            if isempty(obj.coordLFPchanID)
                obj.IdentifyCoordinatingLFP;
            end
            if isempty(obj.lfpMatrix)
                obj.CompileLFPmatrix;
            end
            %% Now figure out which window is being used
            switch alignment
                case obj.beginTrialAlignment
                    window = obj.beginTrialWindow;
                case obj.endTrialAlignment
                    window = obj.endTrialWindow;
            end            
            paddedWindow = [window(1)-(obj.binSize/2) window(2)+(obj.binSize/2)];            
            
            %% Determine the time vector for the trial period in question            
            trialTime = obj.ExtractTrialMatrix(obj.tsVect, paddedWindow, alignment);
            firstValid = 1;
            switch alignment
                case 'Odor'
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).OdorIndex);
                case 'PokeIn'
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).PokeInIndex);
                case 'PokeOut'
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).PokeOutIndex);
                case 'FrontReward'
                    firstValid = find(~isnan([obj.trialInfo.RewardIndex]),1, 'first');
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).RewardIndex);
                case 'RewardSignal'
                    firstValid = find(~isnan([obj.trialInfo.RewardSignalIndex]),1, 'first');
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).RewardSignalIndex);
                case 'RearReward'
                    firstValid = find(~isnan([obj.trialInfo.RearRewardIndex]),1, 'first');
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).RearRewardIndex);
                case 'ErrorSignal'
                    firstValid = find(~isnan([obj.trialInfo.ErrorIndex]),1, 'first');
                    firstAlignNdx = obj.tsVect(obj.trialInfo(firstValid).ErrorIndex);
            end
            tempTS = trialTime(:,1,firstValid) - firstAlignNdx;
            
            %% Extract LFP data
            tempLFP = obj.lfpMatrix(:,obj.coordLFPchanID);
            % Filter the LFP
            [tempSig, tempPhase, tempPower] = obj.SimpleFilter(tempLFP, obj.oscBandLims);
            % Extract LFP and Spiking into time X uni X trial matrices
            tempSigMtx = obj.ExtractTrialMatrix(tempSig, paddedWindow, alignment);
            tempPhaseMtx = obj.ExtractTrialMatrix(tempPhase, paddedWindow, alignment);
            tempPowerMtx = obj.ExtractTrialMatrix(tempPower, paddedWindow, alignment);
            tempSpkMtx = obj.ExtractTrialMatrix(obj.ensembleMatrix, paddedWindow, alignment);
            % Remove any units being ignored from the population
            if ~isempty(obj.popVectIncludeLog)
                tempSpkMtx(:,~obj.popVectIncludeLog,:) = [];
            end
            % Now go through and bin the spiking data based on phase
            binnedMtxs = repmat({nan(size(tempSpkMtx))}, [1, size(obj.phaseBins,1)]);
            tempPhaseIdMtx = nan(size(tempSpkMtx,1),1,size(tempSpkMtx,3));
            for trl = 1:size(tempSpkMtx,3)
                tempTrlPhase = tempPhaseMtx(:,:,trl);
                tempTrlPhaseBinIDs = nan(size(tempTrlPhase));
                for t = 1:size(tempSpkMtx,1)
                    if tempTrlPhase(t)>=max(obj.phaseBins(:,1)) || tempTrlPhase(t)<min(obj.phaseBins(:,2))
                        tempTrlPhaseBinIDs(t) = find(obj.phaseBins(:,1)==max(obj.phaseBins(:,1)));
                    else
                        tempTrlPhaseBinIDs(t) = find((tempTrlPhase(t)>=obj.phaseBins(:,1))...
                            & (tempTrlPhase(t)<obj.phaseBins(:,2)));
                    end                    
                end
                tempPhaseIdMtx(:,:,trl) = tempTrlPhaseBinIDs;
                for u = 1:size(tempSpkMtx,2)
                    for prd = 1:size(obj.phaseBins,1)
                        tempSpkVect = tempSpkMtx(:,u,trl);
                        tempSpkVect(tempTrlPhaseBinIDs~=prd) = 0;
                        binnedMtxs{prd}(:,u,trl) = conv(tempSpkVect,ones(1,obj.binSize)./(obj.binSize/obj.sampleRate), 'same');
                    end
                end
            end
            % Remove the data pad
            unpaddedBinnedMtxs = cell(size(binnedMtxs));
            for bin = 1:size(obj.phaseBins,1)
                unpaddedBinnedMtxs{bin} = binnedMtxs{bin}((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            end
            unpaddedTS = tempTS((obj.binSize/2)+1:end-(obj.binSize/2));
            unpaddedSigMtx = tempSigMtx((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            unpaddedPhaseMtx = tempPhaseMtx((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            unpaddedPowerMtx = tempPowerMtx((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            unpaddedPhaseIdMtx = tempPhaseIdMtx((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            
            % Downsample as necessary
            dsVect = downsample(1:length(unpaddedTS),obj.dsRate);
            dataMtx = cell(size(unpaddedBinnedMtxs));
            for bin = 1:size(obj.phaseBins,1)
                dataMtx{bin} = unpaddedBinnedMtxs{bin}(dsVect,:,:);
            end                
            % Create the phase log vector
            tempPhaseLog = cell2mat(cellfun(@(a,b)ones(size(a,1),1)*b, dataMtx, num2cell(1:4), 'uniformoutput',0));
            tempTrlPrdLog = sum([unpaddedTS<0, (unpaddedTS>0)*2],2);
            
            % Store the data!
            switch alignment
                case obj.beginTrialAlignment
                    obj.beginTrialLFPsig = unpaddedSigMtx(dsVect,:,:);
                    obj.beginTrialLFPphase = unpaddedPhaseMtx(dsVect,:,:);
                    obj.beginTrialLFPpower = unpaddedPowerMtx(dsVect,:,:);
                    obj.beginTrialPhaseLog = tempPhaseLog(:);
                    obj.beginTrialSpkPhaseMtx = cell2mat(dataMtx');
                    obj.beginTrialTime = repmat(unpaddedTS(dsVect), [size(obj.phaseBins,1),1]);
                    obj.beginTrialPeriodLog = repmat(tempTrlPrdLog(dsVect), [size(obj.phaseBins,1),1]);
                    obj.beginTrialLFPphaseID = repmat(unpaddedPhaseIdMtx(dsVect,:,:), [size(obj.phaseBins,1),1]);
                case obj.endTrialAlignment
                    obj.endTrialLFPsig = unpaddedSigMtx(dsVect,:,:);
                    obj.endTrialLFPphase = unpaddedPhaseMtx(dsVect,:,:);
                    obj.endTrialLFPpower = unpaddedPowerMtx(dsVect,:,:);
                    obj.endTrialPhaseLog = tempPhaseLog(:);
                    obj.endTrialSpkPhaseMtx = cell2mat(dataMtx');
                    obj.endTrialTime = repmat(unpaddedTS(dsVect), [size(obj.phaseBins,1),1]);
                    tempTrlPrdLog(tempTrlPrdLog==1) = 3;
                    tempTrlPrdLog(tempTrlPrdLog==2) = 4;
                    obj.endTrialPeriodLog = repmat(tempTrlPrdLog(dsVect), [size(obj.phaseBins,1),1]);
                    obj.endTrialLFPphaseID = repmat(unpaddedPhaseIdMtx(dsVect,:,:), [size(obj.phaseBins,1),1]);
            end
        end
        %% Identify LFP
        function IdentifyCoordinatingLFP(obj)
            if isempty(obj.lfpMatrix)
                obj.CompileLFPmatrix;
            end
            obj.coordLFPchanID = find(obj.numUniPerTet==max(obj.numUniPerTet),1,'first');
        end

    end
end