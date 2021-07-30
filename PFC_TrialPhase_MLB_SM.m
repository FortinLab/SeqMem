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
        oscBandLims = [16 32] % Beta by default
        coordLFPchanID
    end
    properties % Data Structures
        beginTrialSpkMtx
        beginTrialLFPsig
        beginTrialLFPphase
        beginTrialLFPpower
        
        endTrialSpkMtx
        endTrialLFPsig
        endTrialLFPphase
        endTrialLFPpower
    end
    properties % Data Organization Vectors
        beginTrialTime
        endTrialTime
    end
    properties
    end
    %%
    methods
        %% Object Creation
        %#ok<*EXIST>
        function obj = PFC_TrialPhase_MLB_SM(fileDir, binSize, dsRate, beginWindow, endWindow, phaseBins)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@MLB_SM(fileDir);
            if nargin >= 2
                if exist('binSize')==1
                    obj.binSize = binSize;
                end
                if exist('dsRate')==1
                    obj.dsRate = dsRate;
                else
                    obj.dsRate = 1;
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
            fprintf('Analyzing\n');
            obj.IdentifyCoordinatingLFP;
            obj.CompileMLBmtx(obj.beginTrialAlignment);
            obj.CompileMLBmtx(obj.endTrialAlignment);
            obj.CompileTrialObsvs;
            obj.CompileFISlikes;
            obj.DecodeFIS_L1O;
            obj.DecodeOS;
            obj.DecodeTAO;
        end
        %% Identify LFP
        function IdentifyCoordinatingLFP(obj)
            if isempty(obj.lfpMatrix)
                obj.CompileLFPmatrix;
            end
            obj.coordLFPchanID = find(obj.numUniPerTet==max(obj.numUniPerTet),1,'first');
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
            
            % Downsample as necessary
            dsVect = downsample(1:length(unpaddedTS),obj.dsRate);
            dataMtx = cell(size(unpaddedBinnedMtxs));
            for bin = 1:size(obj.phaseBins,1)
                dataMtx{bin} = unpaddedBinnedMtxs{bin}(dsVect,:,:);
            end                
            % Store the data!
            switch alignment
                case obj.beginTrialAlignment
                    obj.beginTrialLFPsig = unpaddedSigMtx(dsVect,:,:);
                    obj.beginTrialLFPphase = unpaddedPhaseMtx(dsVect,:,:);
                    obj.beginTrialLFPpower = unpaddedPowerMtx(dsVect,:,:);
                    obj.beginTrialSpkMtx = dataMtx;
                    obj.beginTrialTime = unpaddedTS(dsVect);
                case obj.endTrialAlignment
                    obj.endTrialLFPsig = unpaddedSigMtx(dsVect,:,:);
                    obj.endTrialLFPphase = unpaddedPhaseMtx(dsVect,:,:);
                    obj.endTrialLFPpower = unpaddedPowerMtx(dsVect,:,:);
                    obj.endTrialSpkMtx = dataMtx;
                    obj.endTrialTime = unpaddedTS(dsVect);
            end
        end
    end
end