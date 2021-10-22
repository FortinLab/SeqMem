classdef MLB_SM < SeqMem
    properties % Analysis Variables
        binSize
        dsRate
        binType = 'box'
    end
    properties % Data Structures
        fiscTrials
    end
    methods
        function obj = MLB_SM(fileDir)
            if nargin == 0 
                fileDir = uigetdir;
            end
            obj@SeqMem(fileDir);
            obj.PP_IdentifyFISCseqs
        end
    end
    methods % Data Pre-Processing Methods
        %% Pre-Process Trial Spiking Data
        function [dataMtx, timeVect] = PP_TrialMatrix_Spiking(obj, window, alignment)
            if isempty(obj.binSize)
                error('Define binSize');
            elseif isempty(obj.dsRate)
                error('Define dsRate');
            end
            paddedWindow = [window(1)-(obj.binSize/2) window(2)+(obj.binSize/2)];
            trialTime = obj.ExtractTrialMatrix(obj.tsVect, paddedWindow, alignment);
            tempMtx = obj.ExtractTrialMatrix(obj.ensembleMatrix, paddedWindow, alignment);
            if ~isempty(obj.popVectIncludeLog)
                tempMtx(:,~obj.popVectIncludeLog,:) = [];
            end
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
            binnedMtx = nan(size(tempMtx));
            for t = 1:size(tempMtx,3)
                for u = 1:size(tempMtx,2)
                    if strcmp(obj.binType, 'box')
                        binnedMtx(:,u,t) = conv(tempMtx(:,u,t),ones(1,obj.binSize)./(obj.binSize/obj.sampleRate), 'same');
                    elseif strcmp(obj.binType, 'gauss')
                        binnedMtx(:,u,t) = conv(tempMtx(:,u,t),gausswin(obj.binSize)./(obj.binSize/obj.sampleRate), 'same');
                    end
                end
            end
            unpaddedBinnedMtx = binnedMtx((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            unpaddedTS = tempTS((obj.binSize/2)+1:end-(obj.binSize/2));
            dsVect = downsample(1:size(unpaddedBinnedMtx,1),obj.dsRate);
            dataMtx = unpaddedBinnedMtx(dsVect,:,:);
            timeVect = unpaddedTS(dsVect);
        end
        %% Pre-Process Trial LFP Data
        function [phaseMtx, powerMtx] = PP_TrialMatrix_LFP(obj, freqWin, window, alignment)
            if isempty(obj.binSize)
                error('Define binSize');
            end
            if isempty(obj.dsRate)
                error('Define dsRate');
            end
            if isempty(obj.lfpMatrix)
                obj.CompileLFPmatrix;
            end
            paddedWindow = [window(1)-(obj.binSize/2) window(2)+(obj.binSize/2)];
            if isempty(obj.lfpRefTet)
                obj.lfpRefTet = find(obj.numUniPerTet==max(obj.numUniPerTet), 1, 'first');
            end
            [~, phase, power] = obj.SimpleFilter(obj.lfpMatrix(:,obj.lfpRefTet), freqWin);

            tempPhase = obj.ExtractTrialMatrix(phase, paddedWindow, alignment);
            tempPower = obj.ExtractTrialMatrix(power, paddedWindow, alignment);
            unpaddedPhase = tempPhase((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            unpaddedPower = tempPower((obj.binSize/2)+1:end-(obj.binSize/2),:,:);
            dsVect = downsample(1:size(unpaddedPhase,1),obj.dsRate);
            phaseMtx = unpaddedPhase(dsVect,:,:);
            powerMtx = unpaddedPower(dsVect,:,:);
        end
        %% Identify Fully InSeq Correct Sequences
        function PP_IdentifyFISCseqs(obj)
            odrVect = [obj.trialInfo.Odor];
            posVect = [obj.trialInfo.Position];
            
            pos1 = find(posVect==1);
            potentialSeqs = repmat(pos1,[obj.seqLength,1])+ repmat((0:obj.seqLength-1)', [1,length(pos1)]);
            fiscLog = false(1,size(potentialSeqs,2));
            for seq = 1:size(potentialSeqs,2)
                if potentialSeqs(end,seq) <= length(obj.trialInfo)
                    if sum(odrVect(potentialSeqs(:,seq)) == posVect(potentialSeqs(:,seq))) == obj.seqLength ...
                        && sum([obj.trialInfo(potentialSeqs(:,seq)).Performance]) == obj.seqLength
                        fiscLog(seq) = true;
                    elseif sum(odrVect(potentialSeqs(:,seq))-10 == posVect(potentialSeqs(:,seq))) == obj.seqLength ...
                        && sum([obj.trialInfo(potentialSeqs(:,seq)).Performance]) == obj.seqLength
                        fiscLog(seq) = true;
                    end
                end
            end
            obj.fiscTrials = potentialSeqs(:,fiscLog);     
        end
    end
    methods % MLB Analysis Methods
        %% Calculate MLB
        function post = CalcStaticBayesPost(obj,likely, obsv)
            post = nan(size(obsv,1), size(likely,1), size(obsv,3));
            for trl = 1:size(obsv,3)
                for t = 1:size(obsv,1)
                    p = nan(size(likely));
                    curPopVect = floor(obsv(t,:,trl)*(obj.binSize/obj.sampleRate));
                    curPopFact = factorial(curPopVect);
                    for u = 1:size(likely,2)
                        curAvgUniFR = likely(:,u);
                        p(:,u) = (((obj.binSize/obj.sampleRate).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
                    end
                    p(p==0) = 0.0000000000000000000001;
                    pp = prod(p,2, 'omitnan');
                    ee = exp(-((obj.binSize/obj.sampleRate)*sum(likely,2)));
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
                        tempDecode = find(post(d,:,o)==maxPost(d,o));
                        select = rand(1,length(tempDecode));
                        decode(d,o) = id(tempDecode(select==max(select)));
                    end
                end
            end
        end
        %% Tabluate MLB
        function decode = TabulateBayesPost(~, post, id)
            idS = unique(id);
            decode = nan(size(post,1), size(post,3), length(idS));
            for trl = 1:size(post,3)
                for t = 1:size(post,1)
                    for iD = 1:length(idS)
                        decode(t,trl,iD) = sum(post(t,id==idS(iD),trl));
                    end
                end
            end
        end
    end
end