classdef MLB_SM < SeqMem
    % Implementation of a naive bayes classifier for spiking data
    properties % Analysis Variables
        binSize
        dsRate
        binType = 'box'
        bayesType %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts
        windows
        alignments
    end
    properties % Data Structures
        fiscTrials        
        trialIDs = [{'Time'}, {'Window'}, {'Position'}, {'Odor'}];
        likeTrlSpikes
        likeTrlIDs
        likeTimeVect
        obsvTrlSpikes
        obsvTrlIDs
        obsvTimeVect
        post
        postTrlIDs
        decode
        decodeIDvects
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
        %% Pre-Process Trial Spiking Data as Binary
        function [binMtx, timeVect] = PP_TrialMatrix_SpikeBin(obj, window, alignment)
            [dataMtx, timeVect] = obj.PP_TrialMatrix_Spiking(window, alignment);
            binMtx = dataMtx.*(obj.binSize/obj.sampleRate);
            binMtx(binMtx>=1) = 1;
        end
        %% Pre-Process Trial Spiking Data as Timebin z-scored
        function [timeZMtx, timeVect] = PP_TrialMatrix_SpikeTimeZ(obj,window,alignment)
            [dataMtx, timeVect] = obj.PP_TrialMatrix_Spiking(window, alignment);
            iscLog = [obj.trialInfo.TranspositionDistance]==0 & [obj.trialInfo.Performance]==1;
            cellTimeMeans = mean(dataMtx(:,:,iscLog),3);
            cellTimeVar = std(dataMtx(:,:,iscLog),0,3);
            for trl = 1:size(dataMtx,3)
                dataMtx(:,:,trl) = (dataMtx(:,:,trl)-cellTimeMeans)./cellTimeVar;
            end
            timeZMtx = dataMtx;
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
            fiscSeqs = potentialSeqs(:,fiscLog);
            if size(obj.odrSeqs,1)>=2
                tempSeqs = cell(size(obj.odrSeqs,1),1);
                for list = 1:size(obj.odrSeqs,1)
                    tempSeqs{list} = fiscSeqs(:,odrVect(fiscSeqs(1,:))==obj.odrSeqs(list,1));
                end
                numSeqs = cellfun(@(a)size(a,2),tempSeqs);
                padSize = max(numSeqs)-numSeqs;
                for list = 1:length(numSeqs)
                    tempSeqs{list} = [tempSeqs{list}, nan(obj.seqLength,padSize(list))];
                end
                obj.fiscTrials = cell2mat(tempSeqs);
            else
                obj.fiscTrials = fiscSeqs;
            end
        end
        %% Group Trial Data
        function [ssnSpikes, ssnID] = PP_ConcatTrialData(obj)
            if isempty(obj.bayesType)
                error('Specify bayesian type to be used before processing data');
            end
            trlDta = cell(size(obj.windows,1),1);
            trlTimeVects = cell(size(obj.windows,1),1);
            for win = 1:size(obj.windows,1)
                if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                    [trlDta{win}, trlTimeVects{win}] = obj.PP_TrialMatrix_Spiking(obj.windows{win}, obj.alignments{win});
                elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                    [trlDta{win}, trlTimeVects{win}] = obj.PP_TrialMatrix_SpikeBin(obj.windows{win}, obj.alignments{win});
                elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                    [trlDta{win}, trlTimeVects{win}] = obj.PP_TrialMatrix_SpikeTimeZ(obj.windows{win}, obj.alignments{win});
                end                
            end
            if ~isempty(obj.popVectIncludeLog)
                ssnSpikes = nan(size(cell2mat(trlTimeVects),1), sum(obj.popVectIncludeLog), length(obj.trialInfo));
            else
                ssnSpikes = nan(size(cell2mat(trlTimeVects),1), size(obj.ensembleMatrix,2), length(obj.trialInfo));
            end
            ssnID = nan(size(cell2mat(trlTimeVects),1), 5, length(obj.trialInfo));
            for trl = 1:size(trlDta{win},3)
                tempTrialData = cell(size(obj.windows,1),1);
                tempTrialTimeVect = cell(size(obj.windows,1),1);
                tempTrialWindowVect = cell(size(obj.windows,1),1);
                for win = 1:size(obj.windows,1)
                    tempTrialData{win} = trlDta{win}(:,:,trl);
                    if win==1
                        tempTrialTimeVect{win} = trlTimeVects{win};
                    else
                        tempTrialTimeVect{win} = cumsum([trlTimeVects{win-1}(end) + mode(diff(trlTimeVects{win})); diff(trlTimeVects{win})]);
                    end
                    tempTrialWindowVect{win} = ones(size(trlTimeVects{win},1),1).*win;
                end
                ssnSpikes(:,:,trl) = cell2mat(tempTrialData);
                ssnID(:,1,trl) = round(cell2mat(tempTrialTimeVect)*1000)/1000;
                ssnID(:,2,trl) = cell2mat(tempTrialWindowVect);
                ssnID(:,3,trl) = ones(size(ssnID,1),1).*obj.trialInfo(trl).Position;
                ssnID(:,4,trl) = ones(size(ssnID,1),1).*obj.trialInfo(trl).Odor;
                ssnID(:,5,trl) = ones(size(ssnID,1),1).*obj.trialInfo(trl).TrialNum;
            end
        end
    end
    methods % Configuration Methods
        %% Set Likelihoods as Fully InSeq Correct
        function SetLikes_FISC(obj)
            % Set likelihoods & observations using FISC trials where FISC are likelihoods and all other trials are observations
            [ssnSpikes, ssnID] = obj.PP_ConcatTrialData;
            obj.likeTrlSpikes = cell(size(obj.fiscTrials,1),1,size(obj.fiscTrials,2));
            obj.likeTrlIDs = cell(size(obj.fiscTrials,1),1,size(obj.fiscTrials,2));
            for pos = 1:size(obj.fiscTrials,1)
                for seq = 1:size(obj.fiscTrials,2)
                    if ~isnan(obj.fiscTrials(pos,seq))
                        obj.likeTrlSpikes{pos,seq} = ssnSpikes(:,:,obj.fiscTrials(pos,seq));
                        obj.likeTrlIDs{pos,seq} = ssnID(:,:,obj.fiscTrials(pos,seq));
                    else
                        obj.likeTrlSpikes{pos,seq} = nan(size(ssnSpikes,1), size(ssnSpikes,2));
                        obj.likeTrlIDs{pos,seq} = nan(size(ssnID,1), size(ssnID,2));
                    end
                end
            end
            obj.likeTimeVect = cell2mat(cellfun(@(a){a(:,1)}, obj.likeTrlIDs(:,:,1)));
            obj.decodeIDvects = cell2mat(cellfun(@(a){a(:,1:end-1)}, obj.likeTrlIDs(:,:,1)));
            
            fiscTrls = unique(obj.fiscTrials(~isnan(obj.fiscTrials)));
            nonFiscLog = true(1,size(ssnSpikes,3));
            nonFiscLog(fiscTrls) = false;
            obj.obsvTrlSpikes = ssnSpikes(:,:,nonFiscLog);
            obj.obsvTrlIDs = ssnID(:,:,nonFiscLog);
            obj.obsvTimeVect = ssnID(:,1,1);
        end
        %% Set Likelihoods as All InSeq Correct
        function SetLikes_ISC(obj)
            % Set likelihoods & observations using All ISC trials where ISC are likelihoods and all other trials are Observations
            [ssnSpikes, ssnID] = obj.PP_ConcatTrialData;
            posTrlIDs = [obj.trialInfo.Position];
            odrTrlIDs = [obj.trialInfo.Odor];
            perfTrlIDs = [obj.trialInfo.Performance];
            iscTrls = cell(fliplr(size(obj.odrSeqs)));
            for s = 1:size(obj.odrSeqs,1)
                for p = 1:size(obj.odrSeqs,2)
                    iscTrls{p,s} = find(posTrlIDs==p & odrTrlIDs==obj.odrSeqs(s,p) & perfTrlIDs==1);
                end
            end
            iscTrls = iscTrls(:);
            trlCounts = cellfun(@(a)length(a),iscTrls);
            trlPad = max(trlCounts)-trlCounts;
            for o = 1:length(iscTrls)
                iscTrls{o} = [iscTrls{o}, nan(1,trlPad(o))];
            end
            iscTrls = cell2mat(iscTrls);
            
            obj.likeTrlSpikes = cell(size(iscTrls,1),1,size(iscTrls,2));
            obj.likeTrlIDs = cell(size(iscTrls,1),1,size(iscTrls,2));
            for pos = 1:size(iscTrls,1)
                for seq = 1:size(iscTrls,2)
                    if ~isnan(iscTrls(pos,seq))
                        obj.likeTrlSpikes{pos,seq} = ssnSpikes(:,:,iscTrls(pos,seq));
                        obj.likeTrlIDs{pos,seq} = ssnID(:,:,iscTrls(pos,seq));
                    else
                        obj.likeTrlSpikes{pos,seq} = nan(size(ssnSpikes,1), size(ssnSpikes,2));
                        obj.likeTrlIDs{pos,seq} = nan(size(ssnID,1), size(ssnID,2));
                    end
                end
            end
            obj.likeTimeVect = cell2mat(cellfun(@(a){a(:,1)}, obj.likeTrlIDs(:,:,1)));
            obj.decodeIDvects = cell2mat(cellfun(@(a){a(:,1:end-1)}, obj.likeTrlIDs(:,:,1)));
            
            iscTrls = unique(iscTrls(~isnan(iscTrls)));
            nonIscLog = true(1,size(ssnSpikes,3));
            nonIscLog(iscTrls) = false;
            obj.obsvTrlSpikes = ssnSpikes(:,:,nonIscLog);
            obj.obsvTrlIDs = ssnID(:,:,nonIscLog);
            obj.obsvTimeVect = ssnID(:,1,1);
        end
    end
    methods % MLB Processing Methods
        %% Process via Leave-1-Out leaving out individual trials during calculation
        function Process_LikelyL1O(obj)
            tempPost = cell(size(obj.likeTrlSpikes));
            obj.postTrlIDs = permute(cellfun(@(a)a(1,end),obj.likeTrlIDs), [1,3,2]);
            for odr = 1:size(obj.postTrlIDs,1)
                for seq = 1:size(obj.postTrlIDs,2)
                    if ~isnan(obj.postTrlIDs(odr,seq))
                        tempObsv = obj.likeTrlSpikes{odr,seq};
                        tempLike = obj.likeTrlSpikes;
                        tempLike{odr,:,seq} = nan(size(tempLike{1}));
                        tempLike = cell2mat(tempLike);
                        tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));                        
                        if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Poisson(mean(tempLike,3, 'omitnan'), tempObsv, tempProb);
                        elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Bernoulli(mean(tempLike,3, 'omitnan'), tempObsv, tempProb);
                        elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Gaussian(mean(tempLike,3, 'omitnan'), std(tempLike,0,3), tempObsv, tempProb);
                        end
                    else
                        tempPost{odr,seq} = nan(size(tempObsv,1), size(tempLike,1));
                    end
                end
            end
            obj.post = cell(size(obj.postTrlIDs,1),1);
            for odr = 1:size(obj.postTrlIDs,1)
                obj.post{odr} = cell2mat(tempPost(odr,:,:));
            end
        end
        %% Process via Leave-1-Out Iteratively
        function Process_IterativeLikelyL1O(obj)
            obj.post = cell(size(obj.likeTrlSpikes));
            obj.postTrlIDs = permute(cellfun(@(a)a(1,end),obj.likeTrlIDs), [1,3,2]);
            for pos = 1:size(obj.likeTrlSpikes,1)
                for seq = 1:size(obj.likeTrlSpikes,3)
                    if ~isnan(obj.postTrlIDs(pos,seq))
                        tempObsv = obj.likeTrlSpikes{pos,1,seq};
                        tempLike = obj.likeTrlSpikes;
                        tempLike{pos,:,seq} = nan(size(tempLike{1}));
                        tempLike = cell2mat(tempLike);
                        tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));
                        if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                            obj.post{pos,seq} = obj.CalcIterativeBayesPost_Poisson(mean(tempLike,3, 'omitnan'), tempObsv, obj.decodeIDvects(:,1), obj.decodeIDvects(:,4), tempProb);
                        elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                            error('Not Implemented Yet');
                        elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                            error('Not Implemented Yet');
                        end
                    end
                end
            end
        end
        %% Process all Observations
        function Process_Observes(obj)
            if iscell(obj.likeTrlSpikes)
                tempLike = cell2mat(obj.likeTrlSpikes);
            else
                tempLike = obj.likeTrlSpikes;
            end
            tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));
            if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                obj.post = obj.CalcStaticBayesPost_Poisson(mean(tempLike,3, 'omitnan'), obj.obsvTrlSpikes, tempProb);
            elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                obj.post = obj.CalcStaticBayesPost_Bernoulli(mean(tempLike,3, 'omitnan'), obj.obsvTrlSpikes, tempProb);
            elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                obj.post = obj.CalcStaticBayesPost_Gaussian(mean(tempLike,3, 'omitnan'), std(obj.likeTrlSpikes,0,3,'omitnan'), obj.obsvTrlSpikes, tempProb);
            end
            obj.postTrlIDs = permute(obj.obsvTrlIDs(1,end,:), [1,3,2]);
        end
        %% Process all Observations Iteratively (cross-temporal decoding)
        function Process_IterativeObserves(obj)    
            if iscell(obj.likeTrlSpikes)
                tempLike = cell2mat(obj.likeTrlSpikes);
            else
                tempLike = obj.likeTrlSpikes;
            end
            tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));
            if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                obj.post = obj.CalcIterativeBayesPost_Poisson(mean(tempLike,3, 'omitnan'), obj.obsvTrlSpikes, obj.decodeIDvects(:,1), obj.decodeIDvects(:,4), tempProb);
            elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                error('Not implemented yet');
%                     obj.post{perm} = obj.CalcIterativeBayesPost_Bernoulli(mean(obj.likeTrlSpikes{perm},3, 'omitnan'), obj.obsvTrlSpikes{perm});
            elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                error('Not implemented yet');
%                     obj.post{perm} = obj.CalcIterativeBayesPost_Gaussian(mean(obj.likeTrlSpikes{perm},3, 'omitnan'), std(obj.likeTrlSpikes{perm},0,3), obj.obsvTrlSpikes{perm});
            end
            obj.postTrlIDs = permute(obj.obsvTrlIDs(1,end,:), [1,3,2]);
        end
    end
    methods % MLB Algorithms
        %% Calculate Static MLB Poisson
        function post = CalcStaticBayesPost_Poisson(obj,likely, obsv, prob)
%             tic;
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
                    p(p==0) = 0.00000000000000001;
                    pp = prod(p,2, 'omitnan');
                    ee = exp(-((obj.binSize/obj.sampleRate)*sum(likely,2)));
                    tempPost = prob.*pp.*ee;
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end
%             toc
        end
        %% Calculate Static MLB Bernoulli
        function post = CalcStaticBayesPost_Bernoulli(~,likely, obsv, prob)
            % Developed by M.Saraf... 
            % MS code: (organized neuron X time)
            %   for i = 1: size(expected_prob,2)
            %       indiv_neurons = [];
            %           for n = 1:size(expected_prob,1)
            %               indiv_neurons(n,:) = (expected_prob(n,:).^test_data(n,i)).*((1-expected_prob(n,:)).^(1-test_data(n,i)));
            %           end
            %       posteriors(i,:) = prod(indiv_neurons,1);
            %   end (edited)
            post = nan(size(obsv,1), size(likely,1), size(obsv,3));
            for trl = 1:size(obsv,3)
                for t = 1:size(obsv,1)
                    indiv_neurons = nan(size(likely));
                    for u = 1:size(likely,2)
                        indiv_neurons(:,u) = (likely(:,u).^obsv(t,u,trl)).*((1-likely(:,u)).^(1-obsv(t,u,trl)));
                    end
                    tempPost = prob.*prod(indiv_neurons,2);
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end
        end
        %% Calculate Static MLB Gaussian
        function post = CalcStaticBayesPost_Gaussian(~,meanLikely, varLikely, obsv, prob)
            % Developed by M.Saraf...
            % MS Code: (organized neuron X time)
            % %get the first term:
            %   first_term = prod(1./(sigma.*(2*pi)^0.5),1);
            % %get the second term :
            %   second_term = [];
            %   for i = 1: size(mu,2)
            %       indiv_neurons = [];
            %       for n = 1:size(mu,1)
            %           a = test_data(n,i) - (mu(n,:));
            %           b = sigma(n,:);
            %           indiv_neurons(n,:) = (-0.5).*(a./b).^2;
            %       end
            %       second_term(:,i) = exp(sum(indiv_neurons,1));
            %   end
            %   posteriors = first_term.*second_term;
            varLikely(varLikely==0) = min(varLikely(varLikely~=0));
%             varLikely(varLikely==0) = nan;
%             meanLikely(varLikely==0) = nan;
            post = nan(size(obsv,1), size(meanLikely,1), size(obsv,3));
            for trl = 1:size(obsv,3)
                for t = 1:size(obsv,1)
                    first_term = prod(1./(varLikely*sqrt(2*pi)),2, 'omitnan');
                    second_term = nan(size(meanLikely));
                    for u = 1:size(meanLikely,2)
                        second_term(:,u) = -0.5.*(((obsv(t,u,trl)-meanLikely(:,u))./varLikely(:,u)).^2);
                    end
                    second_term = exp(sum(second_term,2, 'omitnan'));
                    tempPost = prob.*first_term.*second_term;
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end                        
        end
        %% Calculate Iterative MLB Poisson
        % **** TO CREATE ****
        function post = CalcIterativeBayesPost_Poisson(obj, likely, obsv, depVar, grpVar, prob)
            % Calculate posterior probability using a dependant variable and a grouping variable
            %   The dependant variable (mainly time) is the variable that is iterated across
            %   The grouping variable (mainly odor or position) is the variable used as the likelihoods at each level of the dependant variable
            % For example, typical use case here is to look for temporally invariant coding where
            %   depVar = time, i.e. the algorithm will step through different time points
            %   grpVar = odor/position, i.e. the algorithm will then choose likelihoods from different levels of odor or position
            % Note: The input "likely" is a concatenated vector containing depVar & grpVar arranged data. 
            %   Original use case had the likely arranged with spiking bins from trial types arranged by time and concatenated together into a single
            %   vector such that depVar (time) and grpVar (odor/position) are both 1xN vectors with N = (trial time bins X #positions).
            rateScalar = obj.binSize/obj.sampleRate;
            lvlsDepVar = unique(depVar);
            lvlsGrpVar = unique(grpVar);
            post = nan(size(obsv,1), size(likely,1)/length(lvlsGrpVar), length(lvlsGrpVar), size(obsv,3));
            for dv = 1:length(lvlsDepVar)
                tempLikely = nan(length(lvlsGrpVar), size(likely,2));
                tempProb = nan(length(lvlsGrpVar),size(prob,2));
                for gv = 1:length(lvlsGrpVar)
                    tempLikely(gv,:) = likely(depVar==lvlsDepVar(dv) & grpVar==lvlsGrpVar(gv),:);
                    tempProb(gv,:) = prob(depVar==lvlsDepVar(dv) & grpVar==lvlsGrpVar(gv),:);
                end
                for trl = 1:size(obsv,3)
                    %                     tic;
                    for t = 1:size(obsv,1)
                        p = nan(size(tempLikely));
                        curPopVect = floor(obsv(t,:,trl)*rateScalar);
                        curPopFact = factorial(curPopVect);
                        for u = 1:size(tempLikely,2)
                            curAvgUniFR = tempLikely(:,u);
                            p(:,u) = ((rateScalar.*curAvgUniFR).^curPopVect(u))./curPopFact(u);
                        end
                        p(p==0) = 0.00000000000000001;
                        pp = prod(p,2, 'omitnan');
                        ee = exp(-(rateScalar*sum(tempLikely,2)));
                        tempTempPost = tempProb.*pp.*ee;
                        post(t,dv,:,trl) = tempTempPost ./ sum(tempTempPost);
                    end
                    %                     toc
                end
            end
            %             toc
        end
        %% Calculate Iterative MLB Bernoulli
        % **** TO CREATE ****
        function post = CalcIterativeBayesPost_Bernoulli(obj, likely, obsv, depVar, grpVar)
        end
        %% Calculate Iterative MLB Gaussian
        % **** TO CREATE ****
        function post = CalcIterativeBayesPost_Gaussian(obj, likely, obsv, depVar, grpVar)
        end
    end
    methods % Decoding Methods
        %% Decode MLB
        function [decode, maxPost] = DecodeBayesPost(obj, post, id)
            if ~iscell(post)
                if ndims(post) < 4
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
                else
                    % Assumes post is in the structure of ObservTime X LikelyTime X GrpVar X Trial
                    decode = nan(size(post,1),size(post,1),size(post,4));
                    maxPost = nan(size(post,1),size(post,1),size(post,4));
                    for trl = 1:size(post,4)
                        curTrl = post(:,:,:,trl);
                        for obsvTime = 1:size(curTrl,1)
                            for likeTime = 1:size(curTrl,2)
                                maxPost(obsvTime,likeTime,trl) = max(curTrl(obsvTime,likeTime,:));
                                if ~isnan(maxPost(obsvTime,likeTime,trl))
                                    tempDecode = find(curTrl(obsvTime,likeTime,:)==maxPost(obsvTime,likeTime,trl));
                                    select = rand(1,length(tempDecode));
                                    decode(obsvTime,likeTime,trl) = tempDecode(select==max(select));
                                end
                            end
                        end
                    end
                end
            else
                if ndims(post{1}) < 4 % Output from StaticBayes versions are 3d
                    if size(post{1},3) == size(obj.postTrlIDs,2)
                        decode = cell(size(post,1),1);
                        maxPost = cell(size(post,1),1);
                        for o = 1:size(post,1)
                            tempPost = post{o};
                            tempDecodes = nan(size(tempPost,1), size(tempPost,3));
                            tempMax = nan(size(tempPost,1), size(tempPost,3));
                            for trl = 1:size(tempPost,3)
                                for t = 1:size(tempPost,1)
                                    tempMax(t,trl) = max(tempPost(t,:,trl));
                                    tempDecode = find(tempPost(t,:,trl)==tempMax(t,trl));
                                    select = rand(1,length(tempDecode));
                                    tempDecodes(t,trl) = id(tempDecode(select==max(select)));
                                end
                            end
                            decode{o} = tempDecodes;
                            maxPost{o} = tempMax;
                        end
                    elseif size(post,3) == size(obj.postTrlIDs,2)
                        decode = cell(size(post,1),1,size(post,3));
                        maxPost = cell(size(post,1), size(post,3));
                        for o = 1:size(post,1)
                            for seq = 1:size(post,3)
                                curPost = post{o,1,seq};
                                curDecode = nan(size(post{1},1), size(post{1},2));
                                curMax = nan(size(post{1},1), size(post{1},2));
                                if ~isempty(curPost)
                                    for tO = 1:size(curPost,1)
                                        for tL = 1:size(curPost,1)
                                            curMax(tO,tL) = max(curPost(tO,tL,:));
                                            tempDecode = find(curPost(tO,tL,:)==curMax(tO,tL));
                                            select = rand(1,length(tempDecode));
                                            curDecode(tO,tL) = tempDecode(find(select==max(select)));
                                        end
                                    end
                                end
                                decode{o,1,seq} = curDecode;
                                maxPost{o,1,seq} = curMax;
                            end
                        end
                    end         
                elseif ndims(post) == 4 % Output from IterativeBayes versions are 4d
%                     % Assumes post is in the structure of ObservTime X LikelyTime X GrpVar X Trial
%                     decode = cell(size(post,1),1);
%                     maxPost = cell(size(post,1),1);
% %                     for 
%                     decode = nan(size(post,1),size(post,1),size(post,4));
%                     maxPost = nan(size(post,1),size(post,1),size(post,4));
%                     for trl = 1:size(post,4)
%                         curTrl = post(:,:,:,trl);
%                         for obsvTime = 1:size(curTrl,1)
%                             for likeTime = 1:size(curTrl,2)
%                                 maxPost(obsvTime,likeTime,trl) = max(curTrl(obsvTime,likeTime,:));
%                                 if ~isnan(maxPost(obsvTime,likeTime,trl))
%                                     tempDecode = find(curTrl(obsvTime,likeTime,:)==maxPost(obsvTime,likeTime,trl));
%                                     select = rand(1,length(tempDecode));
%                                     decode(obsvTime,likeTime,trl) = tempDecode(select==max(select));
%                                 end
%                             end
%                         end
%                     end
                end
            end
        end
        %% Tabluate MLB
        function decode = TabulateBayesPost(~, post, id)
            idS = unique(id);
            if ~iscell(post)
                decode = nan(size(post,1), size(post,3), length(idS));
                for trl = 1:size(post,3)
                    for t = 1:size(post,1)
                        for iD = 1:length(idS)
                            decode(t,trl,iD) = sum(post(t,id==idS(iD),trl));
                        end
                    end
                end
            else
                decode = cell(length(idS),1);
                for p = 1:size(post,1)
                    tempPost = post{p};
                    tempDecode = nan(size(tempPost,1), size(tempPost,3), length(idS));
                    for trl = 1:size(tempPost,3)
                        for t = 1:size(tempPost,1)
                            for iD = 1:length(idS)
                                tempDecode(t,trl,iD) = sum(tempPost(t,id==idS(iD),trl));
                            end
                        end
                    end
                    decode{p} = tempDecode;
                end
            end
        end
    end
    methods % "Analyses"
        %%
        function dOvrT = CalcDprmVectFromDecode(obj, decodeMtx, trlIDvect)
            % decodeMtx here is a matrix NxM where N=time and M=trials
            % trlIDvect here is a vector 1xM containing the trial IDs
%             trlIDs = sort(obj.odrSeqs(:));
            trlIDs = 1:obj.seqLength;
            dOvrT = nan(size(decodeMtx,1),length(trlIDs));
            for trl = 1:length(trlIDs)
                tempTrlLogVect = trlIDvect==trlIDs(trl);
                for t = 1:size(decodeMtx,1)
                    decodeCounts(1,1) = sum(decodeMtx(t,tempTrlLogVect)==trlIDs(trl));
                    decodeCounts(1,2) = sum(decodeMtx(t,tempTrlLogVect)~=trlIDs(trl));
                    decodeCounts(2,1) = sum(decodeMtx(t,~tempTrlLogVect)==trlIDs(trl));
                    decodeCounts(2,2) = sum(decodeMtx(t,~tempTrlLogVect)~=trlIDs(trl));
                    dOvrT(t,trl) = obj.CalculateDprime(decodeCounts);
                end
            end
        end
        %%
        function [dOvrTraw, dOvrTchance] = CalcDprmMtxFromDecode(obj, decodeMtx, trlIDvect)
            % decodeMtx here is a NxMxP matrix (tensor?) where N = observation time, M = likelihood time and P = trials
            % trlIDvect here is a 1xP matrix containing the trial IDs
%             trlIDs = sort(obj.odrSeqs(:));
            trlIDs = 1:obj.seqLength;
            rndDecodeMtx = nan(size(decodeMtx));
            for trl = 1:size(decodeMtx,3)
                tempDecode = decodeMtx(:,:,trl);
                shuffDecode = sortrows([randperm(numel(tempDecode))',tempDecode(:)]);
                rndDecodeMtx(:,:,trl) = reshape(shuffDecode(:,2), [size(decodeMtx,1),size(decodeMtx,2)]);
            end
            dOvrTraw = nan(size(decodeMtx,1),size(decodeMtx,2), length(trlIDs));
            dOvrTchance = nan(size(decodeMtx,1),size(decodeMtx,2), length(trlIDs));
            for trl = 1:length(trlIDs)
                tempTrlLogVect = trlIDvect==trlIDs(trl);
                for to = 1:size(decodeMtx,1)
                    for tl = 1:size(decodeMtx,2)
                        decodeCounts(1,1) = sum(decodeMtx(to,tl,tempTrlLogVect)==trlIDs(trl));
                        decodeCounts(1,2) = sum(decodeMtx(to,tl,tempTrlLogVect)~=trlIDs(trl));
                        decodeCounts(2,1) = sum(decodeMtx(to,tl,~tempTrlLogVect)==trlIDs(trl));
                        decodeCounts(2,2) = sum(decodeMtx(to,tl,~tempTrlLogVect)~=trlIDs(trl));
                        dOvrTraw(to,tl,trl) = obj.CalculateDprime(decodeCounts);
                        randCounts(1,1) = sum(rndDecodeMtx(to,tl,tempTrlLogVect)==trlIDs(trl));
                        randCounts(1,2) = sum(rndDecodeMtx(to,tl,tempTrlLogVect)~=trlIDs(trl));
                        randCounts(2,1) = sum(rndDecodeMtx(to,tl,~tempTrlLogVect)==trlIDs(trl));
                        randCounts(2,2) = sum(rndDecodeMtx(to,tl,~tempTrlLogVect)~=trlIDs(trl));
                        dOvrTchance(to,tl,trl) = obj.CalculateDprime(randCounts);
                    end
                end
            end
        end
        %% Integrate Anti-diagonal
        function int = IntegrateAntiDiagonal(~, mtx)
            int = nan(1,size(mtx,1)*2);
            for t = 1:2:size(mtx,1)*2
                int(t) = trapz(diag(flip(mtx), t-size(mtx,1)));
            end
            int(isnan(int)) = [];
        end
        %% Quantify Persistance
        function [avg, up, down, left, right] = QuantPersist(obj, mtx, threshold)
            % Designed to run with output of decoding of iterative bayes
            % Assumes matrix of decodings organized by Observation Time X Likelihood (training) Time
            up = nan(1,size(mtx,1));
            down = nan(1,size(mtx,1));
            left = nan(1,size(mtx,1));
            right = nan(1,size(mtx,1));
            for t = 1:size(mtx,1)
                fwdLike = mtx(t:end,t);
                frstFwdLikeDiff = find(abs(fwdLike-fwdLike(1))>=threshold,1,'first')-1;
                if ~isempty(frstFwdLikeDiff)
                    up(t) = frstFwdLikeDiff;
                end
                revLike = flipud(mtx(1:t,t));
                frstRevLikeDiff = find(abs(revLike-revLike(1))>=threshold,1,'first')-1;
                if ~isempty(frstRevLikeDiff)
                    down(t) = frstRevLikeDiff;
                end
                fwdObsv = mtx(t,t:end);
                frstFwdObsvDiff = find(abs(fwdObsv-fwdObsv(1))>=threshold,1,'first')-1;
                if ~isempty(frstFwdObsvDiff)
                    right(t) = frstFwdObsvDiff;
                end
                revObsv = fliplr(mtx(t,1:t));
                frstRevObsvDiff = find(abs(revObsv-revObsv(1))>=threshold,1,'first')-1;
                if ~isempty(frstRevObsvDiff)
                    left(t) = frstRevObsvDiff;
                end                    
            end
            avg = max([up(:),down(:),right(:),left(:)],[],2);
%             avg = mean([up(:),down(:),right(:),left(:)], 2, 'omitnan');
        end
        %% Visualize Trial Data
        function PlotTrial(obj, spkMtx, lfpMtx, trlInfo)
            for t = 1:length(trlInfo)
                rwdLat = (trlInfo(t).RewardIndex-trlInfo(t).PokeInIndex)/obj.sampleRate;
                figure;
                subplot(4,1,1);
                plot(obj.obsvTimeVect, lfpMtx(:,:,t), '-k');
                set(gca, 'ylim', [min(lfpMtx(:)) max(lfpMtx(:))]);
                hold on;
                box off;
                plot([0 0], get(gca, 'ylim'), '--k', 'linewidth', 2);
                plot(repmat(trlInfo(t).TargetDuration, [1,2]), get(gca, 'ylim'), '-k', 'linewidth', 2);
                plot(repmat(trlInfo(t).PokeDuration, [1,2]), get(gca, 'ylim'), '--k', 'linewidth', 2);
                plot(repmat(rwdLat, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
                title(sprintf('Trial#=%i; Odor=%i; Position=%i; Performance=%i', trlInfo(t).TrialNum, trlInfo(t).Odor, trlInfo(t).Position, trlInfo(t).Performance));
                                
                subplot(4,1,2:4);
                imagesc(obj.obsvTimeVect, 1:size(spkMtx,2), spkMtx(:,:,t)', [0 10*(1/(obj.binSize/obj.sampleRate))]);
%                 set(gca, 'ydir', 'normal');
                hold on;
                plot([0 0], [1, size(spkMtx,2)], '--k', 'linewidth', 2);
                plot(repmat(trlInfo(t).TargetDuration, [1,2]), [1,size(spkMtx,2)], '-k', 'linewidth', 2);
                plot(repmat(trlInfo(t).PokeDuration, [1,2]), [1,size(spkMtx,2)], '--k', 'linewidth', 2);
                plot(repmat(rwdLat, [1,2]), [1,size(spkMtx,2)], ':k', 'linewidth', 2);
                drawnow;
            end
        end
        %% Visualize PCA
        % NOTE: This will move to it's own class eventually
        function PFC_PCA(obj)
            curData = obj.obsvTrlSpikes{1};
            curTrlInfo = obj.trialInfo(obj.obsvTrlIDs{1});
            
%             curData = curData(:,:,obj.fiscTrials(:));
%             curTrlInfo = curTrlInfo(obj.fiscTrials(:));            
            
            perfLog = [curTrlInfo.Performance]==1;
            isLog = [curTrlInfo.TranspositionDistance]==0;
            posVect = [curTrlInfo.Position];
            dims = [1 2 3];
            
            trlPV = cell(size(curData,3),1);
            for t = 1:size(curData,3)
                trlPV{t} = reshape(curData(:,:,t), [1,numel(curData(:,:,t))]);
                %     trlPV{t} = (mean(mlb.obsvTrlSpikes{1}(:,:,t),2))';
            end
            
            [~,score] = pca(cell2mat(trlPV(perfLog,:)));
%             [~,score] = pca(cell2mat(trlPV));
            figure;
            h = gca;                
            hold on;
%             scatter3(score(~perfLog & isLog,dims(1)), score(~perfLog & isLog,dims(2)), score(~perfLog & isLog,dims(3)), 'marker', 'x', 'markeredgecolor', 'k');
%             scatter3(score(~perfLog & ~isLog,dims(1)), score(~perfLog & ~isLog,dims(2)), score(~perfLog & ~isLog,dims(3)),'marker', 'x', 'markeredgecolor', 'r');
            isLog = isLog(perfLog);
            posVect = posVect(perfLog);
            perfLog = perfLog(perfLog);
            for pos = 1:4
                curPosLog = perfLog & isLog & posVect==pos;
                scatter3(score(curPosLog,dims(1)), score(curPosLog,dims(2)), score(curPosLog,dims(3)), 20, 'filled', 'marker', 'o', 'markerfacecolor', obj.PositionColors(pos,:));
                scatter3(mean(score(curPosLog,dims(1))), mean(score(curPosLog,dims(2))), mean(score(curPosLog,dims(3))), 100, 'filled', 'marker', 'o', 'markerfacecolor', obj.PositionColors(pos,:));
            end
            h.XRuler.FirstCrossoverValue = 0;
            h.XRuler.SecondCrossoverValue = 0;
            h.YRuler.FirstCrossoverValue = 0;
            h.YRuler.SecondCrossoverValue = 0;
            h.ZRuler.FirstCrossoverValue = 0;
            h.ZRuler.SecondCrossoverValue = 0;
            drawnow;
        end
        %% Validate PCA
        % NOTE: This will move to it's own class eventually
        function distMtx = PFC_validPCA(obj)
            distMtx = nan(obj.seqLength, obj.seqLength, obj.numPerms);
            for perm = 1:obj.numPerms
                trainData = obj.likeTrlSpikes{perm};
                %             trainPV = nan(size(trainData,3), size(trainData,1)*size(trainData,2));
                trainPV = cell(size(trainData,3), 1);
                posNdx = [find(obj.likeTimeVect==min(obj.likeTimeVect)); length(obj.likeTimeVect)+1];
                for t = 1:size(trainData,3)
                    %                 trainPV(t,:) = reshape(trainData(:,:,t), [1,numel(trainData(:,:,t))]);
                    tempSeqPV = trainData(:,:,t);
                    tempPV = nan(obj.seqLength, length(obj.obsvTimeVect)*size(tempSeqPV,2));
                    for p = 1:obj.seqLength
                        tempPV(p,:) = reshape(tempSeqPV(posNdx(p):posNdx(p+1)-1,:), [1,length(obj.obsvTimeVect)*size(tempSeqPV,2)]);
                    end
                    trainPV{t} = tempPV;
                end
                
                testData = obj.obsvTrlSpikes{perm};
                testPV = nan(size(testData,3), size(testData,1)*size(testData,2));
                for t = 1:size(testData,3)
                    testPV(t,:) = reshape(testData(:,:,t), [1,numel(testData(:,:,t))]);
                end
                
                % Train Likelihoods - Test Observations
                [coeff,~,~,~,explained,mu] = pca(cell2mat(trainPV));
                dim = find(cumsum(explained)>95,1, 'first');                
                testScores = (testPV-mu)*coeff;
                curTrlInfo = obj.trialInfo(obj.obsvTrlIDs{perm});
                % Train Observations - Test Likelihoods
%                 [coeff,~,~,~,explained,mu] = pca(testPV);
%                 dim = find(cumsum(explained)>95,1, 'first');                
%                 testScores = (cell2mat(trainPV)-mu)*coeff;
%                 curTrlInfo = obj.trialInfo(obj.likeTrlIDs{perm});
                                
                perfLog = [curTrlInfo.Performance]==1;
                isLog = [curTrlInfo.TranspositionDistance]==0;
                posVect = [curTrlInfo.Position];
                dims = [1 2 3];
                figure;
                h = gca;
                hold on;
                %             scatter3(score(~perfLog & isLog,dims(1)), score(~perfLog & isLog,dims(2)), score(~perfLog & isLog,dims(3)), 'marker', 'x', 'markeredgecolor', 'k');
                %             scatter3(score(~perfLog & ~isLog,dims(1)), score(~perfLog & ~isLog,dims(2)), score(~perfLog & ~isLog,dims(3)),'marker', 'x', 'markeredgecolor', 'r');
                isLog = isLog(perfLog);
                posVect = posVect(perfLog);
                perfLog = perfLog(perfLog);
                simMtx = squareform(pdist(testScores(:,1:dim), 'cosine'));
                for pos = 1:obj.seqLength
                    curPosLog = perfLog & isLog & posVect==pos;
                    scatter3(testScores(curPosLog,dims(1)), testScores(curPosLog,dims(2)), testScores(curPosLog,dims(3)), 20, 'filled', 'marker', 'o', 'markerfacecolor', obj.PositionColors(pos,:));
                    scatter3(mean(testScores(curPosLog,dims(1))), mean(testScores(curPosLog,dims(2))), mean(testScores(curPosLog,dims(3))), 100, 'filled', 'marker', 'o', 'markerfacecolor', obj.PositionColors(pos,:));
                    for odr = 1:obj.seqLength
                        tempDist = simMtx(curPosLog,(perfLog & isLog & posVect==odr));
                        if pos==odr
                            tempDist = tempDist(logical(triu(ones(size(tempDist)),1)));
                        end
                        distMtx(pos,odr,perm) = mean(tempDist(:));
                    end
                end
                h.XRuler.FirstCrossoverValue = 0;
                h.XRuler.SecondCrossoverValue = 0;
                h.YRuler.FirstCrossoverValue = 0;
                h.YRuler.SecondCrossoverValue = 0;
                h.ZRuler.FirstCrossoverValue = 0;
                h.ZRuler.SecondCrossoverValue = 0;
                title(dim);
                drawnow;
            end
            figure; 
            imagesc(1:obj.seqLength, 1:obj.seqLength, mean(distMtx,3));
        end
    end
end