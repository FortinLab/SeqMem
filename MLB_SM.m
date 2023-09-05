classdef MLB_SM < SeqMem
    % Implementation of a naive bayes classifier for spiking data
    properties % Analysis Variables
        binSize
        dsRate
        binType = 'box'
        bayesType %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts
        windows
        alignments
        alignmentOffset = 0
    end
    properties % Data Structures
        fiscTrials        
        trialIDs = [{'Time'}, {'Window'}, {'Position'}, {'Odor'}];
        likeTrlSpikes
        likeIDvects
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
    properties % Trial Event Masks
        mask_Trial
        mask_Int
        mask_IntMid
        mask_PreInt
        mask_PreIntCON
        mask_PstInt
        mask_PstIntCON
    end
    methods
        function obj = MLB_SM(fileDir)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@SeqMem(fileDir);
            fprintf('Completed\n');
            obj.PP_IdentifyFISCseqs
        end
    end
    %% Data Pre-Processing Methods
    methods
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
    %% MLB Configuration Methods
    methods
        %% Set Likelihoods as Fully InSeq Correct
        function SetLikes_FISC(obj)
            % Set likelihoods & observations using FISC trials where FISC are likelihoods and all other trials are observations
            [ssnSpikes, ssnID] = obj.PP_ConcatTrialData;
            obj.likeTrlSpikes = cell(size(obj.fiscTrials,1),1,size(obj.fiscTrials,2));
            obj.likeIDvects = cell(size(obj.fiscTrials,1),1,size(obj.fiscTrials,2));
            for pos = 1:size(obj.fiscTrials,1)
                for seq = 1:size(obj.fiscTrials,2)
                    if ~isnan(obj.fiscTrials(pos,seq))
                        obj.likeTrlSpikes{pos,seq} = ssnSpikes(:,:,obj.fiscTrials(pos,seq));
                        obj.likeIDvects{pos,seq} = ssnID(:,:,obj.fiscTrials(pos,seq));
                    else
                        obj.likeTrlSpikes{pos,seq} = nan(size(ssnSpikes,1), size(ssnSpikes,2));
                        obj.likeIDvects{pos,seq} = nan(size(ssnID,1), size(ssnID,2));
                    end
                end
            end
            obj.likeTrlIDs = nan(size(obj.likeIDvects));
            for ind = 1:numel(obj.likeIDvects)
                if ~isempty(obj.likeIDvects{ind})
                    obj.likeTrlIDs(ind) = obj.likeIDvects{ind}(1,end);
                end
            end
            obj.likeTimeVect = cell2mat(cellfun(@(a){a(:,1)}, obj.likeIDvects(:,:,1)));
            obj.decodeIDvects = cell2mat(cellfun(@(a){a(:,1:end-1)}, obj.likeIDvects(:,:,1)));
            
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
            obj.likeIDvects = cell(size(iscTrls,1),1,size(iscTrls,2));
            for pos = 1:size(iscTrls,1)
                for seq = 1:size(iscTrls,2)
                    if ~isnan(iscTrls(pos,seq))
                        obj.likeTrlSpikes{pos,seq} = ssnSpikes(:,:,iscTrls(pos,seq));
                        obj.likeIDvects{pos,seq} = ssnID(:,:,iscTrls(pos,seq));
                    else
                        obj.likeTrlSpikes{pos,seq} = nan(size(ssnSpikes,1), size(ssnSpikes,2));
                        obj.likeIDvects{pos,seq} = nan(size(ssnID,1), size(ssnID,2));
                    end
                end
            end
            obj.likeTrlIDs = nan(size(obj.likeIDvects));
            for ind = 1:numel(obj.likeIDvects)
                if ~isempty(obj.likeIDvects{ind})
                    obj.likeTrlIDs(ind) = obj.likeIDvects{ind}(1,end);
                end
            end
            obj.likeTimeVect = cell2mat(cellfun(@(a){a(:,1)}, obj.likeIDvects(:,:,1)));
            obj.decodeIDvects = cell2mat(cellfun(@(a){a(:,1:end-1)}, obj.likeIDvects(:,:,1)));
            
            iscTrls = unique(iscTrls(~isnan(iscTrls)));
            nonIscLog = true(1,size(ssnSpikes,3));
            nonIscLog(iscTrls) = false;
            obj.obsvTrlSpikes = ssnSpikes(:,:,nonIscLog);
            obj.obsvTrlIDs = ssnID(:,:,nonIscLog);
            obj.obsvTimeVect = ssnID(:,1,1);
        end
        %% Set Observations as OutSeq Correct
        function SetObsvs_OSC(obj)
            % I need to be made!!!
        end
    end
    %% MLB Processing Methods
    methods 
        %% Permute the likelihood spiking data
        function [permLikes, permLikesTrlIDs] = RandPermLikes(obj, shuffType)
            tempLikes = obj.likeTrlSpikes;
            % Permute Trial IDs
            tempLikesPos = tempLikes;
            realTrlIDsNdx = find(~isnan(obj.likeTrlIDs));
            permTrlIDsNdx = realTrlIDsNdx(randperm(length(realTrlIDsNdx))');
            tempLikesPos(realTrlIDsNdx) = tempLikes(permTrlIDsNdx);
            % Permute Time IDs
            tempLikesTime = tempLikes;
            tempLikesFull = tempLikesPos;
            for trl = 1:numel(tempLikes)
                chancePerm = randperm(length(obj.obsvTimeVect));
                tempLikesTime{trl} = tempLikesTime{trl}(chancePerm,:);
                tempLikesFull{trl} = tempLikesFull{trl}(chancePerm,:);
            end
            permLikesTrlIDs = obj.likeTrlIDs;
            if strcmp(shuffType, 'Trial')
                permLikes = tempLikesPos;
                permLikesTrlIDs(realTrlIDsNdx) = permLikesTrlIDs(permTrlIDsNdx);
            elseif strcmp(shuffType, 'Time')
                permLikes = tempLikesTime;
            elseif strcmp(shuffType, 'Full')
                permLikes = tempLikesFull;
                permLikesTrlIDs(realTrlIDsNdx) = permLikesTrlIDs(permTrlIDsNdx);
            end
        end
        %% Process via Leave-1-Out leaving out individual trials during calculation
        function Process_LikelyL1O(obj, shuffType)
            if nargin == 1
                shuffYN = false;
            elseif nargin == 2
                shuffYN = true;
            end
            tempPost = cell(size(obj.likeTrlSpikes));
            obj.postTrlIDs = permute(cellfun(@(a)a(1,end),obj.likeIDvects), [1,3,2]);
            tempObsvs = obj.likeTrlSpikes;
            if shuffYN
                [tempLikes, tempTrlIDs] = obj.RandPermLikes(shuffType);
            else
                tempLikes = obj.likeTrlSpikes;
                tempTrlIDs = obj.postTrlIDs;
            end
            for odr = 1:size(obj.postTrlIDs,1)
                for seq = 1:size(obj.postTrlIDs,2)
                    if ~isnan(obj.postTrlIDs(odr,seq))
                        currObsv = tempObsvs{odr,seq};
                        curLikes = tempLikes;
                        curLikes{find(tempTrlIDs==obj.postTrlIDs(odr,seq))} = nan(size(currObsv)); %#ok<FNDSB> 
                        curLikes = cell2mat(curLikes);
                        tempProb = sum(~isnan(curLikes(:,1,:)),3)./sum(sum(~isnan(curLikes(:,1,:))));
                        if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Poisson(mean(curLikes,3, 'omitnan'), currObsv, tempProb);
                        elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Bernoulli(mean(curLikes,3, 'omitnan'), currObsv, tempProb);
                        elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                            tempPost{odr,seq} = obj.CalcStaticBayesPost_Gaussian(mean(curLikes,3, 'omitnan'), std(curLikes,0,3), currObsv, tempProb);
                        end
                    else
                        tempPost{odr,seq} = nan(size(currObsv,1), size(curLikes,1));
                    end
                end
            end
            obj.post = cell(size(obj.postTrlIDs,1),1);
            for odr = 1:size(obj.postTrlIDs,1)
                obj.post{odr} = cell2mat(tempPost(odr,:,:));
            end
        end
        %% Process via Leave-1-Out Iteratively
        function Process_IterativeLikelyL1O(obj, shuffType)
            if nargin == 1
                shuffYN = false;
            elseif nargin ==2
                shuffYN = true;
            end
            tempPost = cell(size(obj.likeTrlSpikes));
            obj.postTrlIDs = permute(cellfun(@(a)a(1,end),obj.likeIDvects), [1,3,2]);
            tempObsvs = obj.likeTrlSpikes;
            if shuffYN 
                [tempLikes, tempTrlIDs] = obj.RandPermLikes(shuffType);
            else
                tempLikes = obj.likeTrlSpikes;
                tempTrlIDs = obj.postTrlIDs;
            end
            for pos = 1:size(tempLikes,1)
                for seq = 1:size(tempLikes,3)
                    if ~isnan(obj.postTrlIDs(pos,seq))
                        cur_TempObsv = tempObsvs{pos,1,seq};
                        cur_TempLike = tempLikes;
                        cur_TempLike{find(tempTrlIDs==obj.postTrlIDs(pos,seq))} = nan(size(cur_TempObsv)); %#ok<FNDSB> 
                        cur_TempLike = cell2mat(cur_TempLike);
                        tempProb = sum(~isnan(cur_TempLike(:,1,:)),3)./sum(sum(~isnan(cur_TempLike(:,1,:))));
                        if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                            tempPost{pos,seq} = obj.CalcIterativeBayesPost_Poisson(mean(cur_TempLike,3, 'omitnan'), cur_TempObsv, obj.decodeIDvects(:,1), obj.decodeIDvects(:,4), tempProb);
                        elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                            error('Not Implemented Yet');
                        elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                            error('Not Implemented Yet');
                        end
                    end
                end
            end
            obj.post = tempPost;
        end
        %% Process all Observations
        function Process_Observes(obj, shuffType)
            if nargin == 1
                shuffYN = false;
            elseif nargin == 2
                shuffYN = true;
            end
            if shuffYN
                [tempLike, ~] = obj.RandPermLikes(shuffType);
            else
                tempLike = obj.likeTrlSpikes;
            end
            if iscell(obj.likeTrlSpikes)
                tempLike = cell2mat(tempLike);
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
            % Pull out Trial IDs
            trlIDs = squeeze(obj.obsvTrlIDs(1,end,:));
            trlPosIDs = [obj.trialInfo(trlIDs).Position];
            % Identify maximum number of sequence trials
            maxNumSeqs = sum(trlPosIDs==1);
            obj.post = cell(obj.seqLength, 1, maxNumSeqs);
            obj.postTrlIDs = nan(obj.seqLength,maxNumSeqs);
            for pos = 1:obj.seqLength
                posIDs = trlIDs(trlPosIDs==pos);
                obj.postTrlIDs(pos,1:length(posIDs)) = posIDs;
                for trl = 1:length(posIDs)
                    tempObsv = obj.obsvTrlSpikes(:,:,trlIDs==posIDs(trl));
                    tempProb = sum(~isnan(tempLike(:,1,:)),3)./sum(sum(~isnan(tempLike(:,1,:))));
                    if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                        obj.post{pos,1,trl} = obj.CalcIterativeBayesPost_Poisson(mean(tempLike,3, 'omitnan'), tempObsv, obj.decodeIDvects(:,1), obj.decodeIDvects(:,4), tempProb);
                    elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                        error('Not implemented yet');
                        %                     obj.post{perm} = obj.CalcIterativeBayesPost_Bernoulli(mean(obj.likeTrlSpikes{perm},3, 'omitnan'), obj.obsvTrlSpikes{perm});
                    elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                        error('Not implemented yet');
                        %                     obj.post{perm} = obj.CalcIterativeBayesPost_Gaussian(mean(obj.likeTrlSpikes{perm},3, 'omitnan'), std(obj.likeTrlSpikes{perm},0,3), obj.obsvTrlSpikes{perm});
                    end
                end
            end
        end
    end
    %% MLB Algorithms
    methods
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
%                                         tic;
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
%                                         toc
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
    %% Posterior Processing Methods
    methods
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
        %% New Decode MLB
        function [decode, maxPost] = DecodeBayesPostNew(~,post,decodeDim)
            if decodeDim==1
                decode = nan(size(post,2), size(post,3));
                maxPost = nan(size(post,2), size(post,3));
                for r = 1:size(post,2)
                    for c = 1:size(post,3)
                        decode(r,c) = find(post(:,r,c)==max(post(:,r,c)),1,'first');
                        maxPost(r,c) = max(post(:,r,c));
                    end
                end
            elseif decodeDim==2
                decode = nan(size(pos,1), size(post,3));
                maxPost = nan(size(post,2), size(post,3));
                for r = 1:size(post,1)
                    for c = 1:size(post,3)
                        decode(r,c) = find(post(r,:,c)==max(post(r,:,c)),1,'first');
                        maxPost(r,c) = max(post(r,:,c));
                    end
                end
            elseif decodeDim==3
                decode = nan(size(post,1), size(post,2));
                maxPost = nan(size(post,1), size(post,2));
                for r = 1:size(post,1)
                    for c = 1:size(post,2)
                        if sum(~isnan(post(r,c,:)),3)~=0
                            decode(r,c) = find(post(r,c,:)==max(post(r,c,:)),1,'first');
                            maxPost(r,c) = max(post(r,c,:));
                        end
                    end
                end
            end
                
        end
        %% Tabluate MLB
        function decode = TabulateBayesPost(~, post, id)
            idS = unique(id);
            if ~iscell(post)
                decode = nan(size(post,1), length(idS), size(post,3));
                for trl = 1:size(post,3)
                    for t = 1:size(post,1)
                        for iD = 1:length(idS)
                            decode(t,iD,trl) = sum(post(t,id==idS(iD),trl));
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
        %% Calculate d' vector from decodings
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
        %% Calculate d' matrix from decodings
        function [dOvrTraw, dOvrTchance] = CalcDprmMtxFromDecode(obj, decodeMtx, trlIDvect)
            % decodeMtx here is a NxMxP array where N = observation time, M = likelihood time and P = trials
            % trlIDvect here is a 1xP matrix containing the trial IDs
            %             trlIDs = sort(obj.odrSeqs(:));
            trlIDs = unique(trlIDvect);
            rndDecodeMtx = nan(size(decodeMtx));
            for trl = 1:size(decodeMtx,3)
                tempDecode = decodeMtx(:,:,trl);
                shuffDecode = sortrows([randperm(numel(tempDecode))',tempDecode(:)]);
                rndDecodeMtx(:,:,trl) = reshape(shuffDecode(:,2), [size(decodeMtx,1),size(decodeMtx,2)]);
            end
            dOvrTraw = repmat({nan(size(decodeMtx,1),size(decodeMtx,2))},length(trlIDs),length(trlIDs));
            dOvrTchance = repmat({nan(size(decodeMtx,1),size(decodeMtx,2))},length(trlIDs),length(trlIDs));
            for trl1 = 1:length(trlIDs)
                tempTrlLogVect = trlIDvect==trlIDs(trl1);
                for trl2 = 1:length(trlIDs)
                    for to = 1:size(decodeMtx,1)
                        for tl = 1:size(decodeMtx,2)
                            decodeCounts(1,1) = sum(decodeMtx(to,tl,tempTrlLogVect)==trlIDs(trl2));
                            decodeCounts(1,2) = sum(decodeMtx(to,tl,tempTrlLogVect)~=trlIDs(trl2));
                            decodeCounts(2,1) = sum(decodeMtx(to,tl,~tempTrlLogVect)==trlIDs(trl2));
                            decodeCounts(2,2) = sum(decodeMtx(to,tl,~tempTrlLogVect)~=trlIDs(trl2));
                            dOvrTraw{trl1,trl2}(to,tl) = obj.CalculateDprime(decodeCounts);
                            randCounts(1,1) = sum(rndDecodeMtx(to,tl,tempTrlLogVect)==trlIDs(trl2));
                            randCounts(1,2) = sum(rndDecodeMtx(to,tl,tempTrlLogVect)~=trlIDs(trl2));
                            randCounts(2,1) = sum(rndDecodeMtx(to,tl,~tempTrlLogVect)==trlIDs(trl2));
                            randCounts(2,2) = sum(rndDecodeMtx(to,tl,~tempTrlLogVect)~=trlIDs(trl2));
                            dOvrTchance{trl1,trl2}(to,tl) = obj.CalculateDprime(randCounts);
                        end
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
        %% Calculate d', HR, FAR from posteriors
        function [tMat_HR, tMat_FAR, tMat_D] = CalcDecodabilityTransMatFromPost(obj)
            tMat_HR = cell(obj.seqLength, obj.seqLength);
            tMat_FAR = cell(obj.seqLength, obj.seqLength);
            tMat_D = cell(obj.seqLength, obj.seqLength);
            % To avoid nans and infs any values ==1 or ==0 need to be replaced with the next nearest value
            allVals = unique([obj.post{:}]);
            if min(allVals)==0
                minVal = allVals(2);
            else
                minVal = allVals(1);
            end
            if max(allVals)==1
                maxVal = allVals(end-1);
            else
                maxVal = allVals(end);
            end
            clear allVals;
            for trlPos = 1:obj.seqLength
                for decPos = 1:obj.seqLength
                    tempHRs = obj.post(trlPos,1,:);
                    tempHRs(cellfun(@(a)isempty(a),tempHRs)) = [];
                    temp_tempHRs = cell2mat(cellfun(@(a){a(:,:,decPos)}, tempHRs));
                    temp_tempHRs(temp_tempHRs==0) = minVal;
                    temp_tempHRs(temp_tempHRs==1) = maxVal;
                    tMat_HR{trlPos, decPos} = temp_tempHRs;

                    faLog = (1:obj.seqLength ~= trlPos);
                    tempFARs = reshape(squeeze(obj.post(faLog,1,:)), [1,1,numel(obj.post)-(size(obj.post,3)*sum(~faLog))]);
                    tempFARs(cellfun(@(a)isempty(a),tempFARs)) = [];
                    temp_tempFARs = repmat(mean(cell2mat(cellfun(@(a){a(:,:,decPos)},tempFARs)),3,'omitnan'), [1,1,size(temp_tempHRs,3)]);
                    temp_tempFARs(temp_tempFARs==0) = minVal;
                    temp_tempFARs(temp_tempFARs==1) = maxVal;
                    tMat_FAR{trlPos,decPos} = temp_tempFARs;

                    tMat_D{trlPos,decPos} = norminv(temp_tempHRs) - norminv(temp_tempFARs);
                end
            end
        end
        %% Organize Posteriors into 3D TransMat
        function [tm_HR, tm_D, tm_TrialInfo] = OrganizeDecodabilityTrialHistoryTransMat(obj,hr,d)
            if nargin == 1
                [hr, ~, d] = obj.CalcDecodabilityTransMatFromPost;
            elseif nargin == 2
                d = hr;
            end
            tm_HR = cell(obj.seqLength,obj.seqLength,obj.seqLength,obj.seqLength);
            tm_D = cell(obj.seqLength,obj.seqLength,obj.seqLength,obj.seqLength);
            tm_TrialInfo = cell(obj.seqLength,obj.seqLength,obj.seqLength);
            seqStarts = find([obj.trialInfo.Position]==1);
            for trlOP = 1:obj.seqLength
                curTrls = obj.postTrlIDs(trlOP,:);
                curTrls(isnan(curTrls)) = [];
                for trl = 1:length(curTrls)
                    cur_TrialInfo = obj.trialInfo(curTrls(trl));
                    if curTrls(trl) ~= 1
                        cur_TrialInfo.PrevTrlEnd = obj.trialInfo(curTrls(trl)-1).PokeOutIndex;
                    else
                        cur_TrialInfo.PrevTrlEnd = nan;
                    end
                    if curTrls(trl) ~= length(obj.trialInfo)
                        cur_TrialInfo.NextTrlStart = obj.trialInfo(curTrls(trl)+1).PokeInIndex;
                    else
                        cur_TrialInfo.NextTrlStart = nan;
                    end
                    curSeqStart = seqStarts(find(seqStarts<=cur_TrialInfo.TrialNum,1,'last'));
                    curSeqVect = [obj.trialInfo(curSeqStart:cur_TrialInfo.TrialNum).Odor];
                    osLog = mod(curSeqVect,10)~=1:trlOP;
                    if sum(osLog)==0
                        osPos = trlOP;
                        osOdr = trlOP;
                    else
                        osPos = find(osLog);
                        osOdr = mod(curSeqVect(osPos),10);
                    end
                    tm_TrialInfo{osOdr,osPos,trlOP} = cat(2,tm_TrialInfo{osOdr,osPos,trlOP}, cur_TrialInfo);
                    for decPos = 1:obj.seqLength
                        tm_HR{osOdr,osPos,trlOP,decPos} = cat(3,tm_HR{osOdr,osPos,trlOP,decPos}, hr{trlOP,decPos}(:,:,trl));
                        tm_D{osOdr,osPos,trlOP,decPos} = cat(3,tm_D{osOdr,osPos,trlOP,decPos}, d{trlOP,decPos}(:,:,trl));
                    end
                end
            end
% 
        end
    end
    %% Mask Creation Methods
    methods
        %% Create TrialInfo Mask
        function Create_TrialInfoMasks(obj)
            % Find trial start anchors
            if strcmp(obj.alignments{1}, 'PokeIn')
                tStart = zeros(size(obj.trialInfo));
                startNdx = repmat(find(obj.obsvTimeVect==0), size(obj.trialInfo));
                tEnd = [obj.trialInfo.PokeDuration];
                endNdx = arrayfun(@(a){find(a<obj.obsvTimeVect,1,'first')},tEnd);
                messedUpDur = find(cellfun(@(a)isempty(a),endNdx));
                for mud = 1:length(messedUpDur)
                    if obj.trialInfo(messedUpDur(mud)).Performance == 0
                        obj.trialInfo(messedUpDur(mud)).PokeOutIndex = obj.trialInfo(messedUpDur(mud)).ErrorIndex;
                        newDur = (obj.trialInfo(messedUpDur(mud)).ErrorIndex - obj.trialInfo(messedUpDur(mud)).PokeInIndex)/obj.sampleRate;
                    elseif obj.trialInfo(messedUpDur(mud)).Performance == 1
                        obj.trialInfo(messedUpDur(mud)).PokeOutIndex = obj.trialInfo(messedUpDur(mud)).RewardSignalIndex;
                        newDur = (obj.trialInfo(messedUpDur(mud)).RewardSignalIndex - obj.trialInfo(messedUpDur(mud)).PokeInIndex)/obj.sampleRate;
                    end
                    tEnd(messedUpDur(mud)) = newDur;
                    endNdx{messedUpDur(mud)} = find(newDur<obj.obsvTimeVect,1,'first');
                end
                endNdx = cell2mat(endNdx);
            elseif strcmp(obj.alignments{1}, 'PokeOut')
                tStart = [obj.trialInfo.PokeDuration]*-1;
                startNdx = arrayfun(@(a){find(a<obj.obsvTimeVect,1,'first')}, tStart);
                messedUpDur = find(cellfun(@(a)isempty(a),startNdx));
                for mud = 1:length(messedUpDur)
                    if obj.trialInfo(messedUpDur(mud)).Performance == 0
                        obj.trialInfo(messedUpDur(mud)).PokeOutIndex = obj.trialInfo(messedUpDur(mud)).ErrorIndex;
                        newDur = (obj.trialInfo(messedUpDur(mud)).ErrorIndex - obj.trialInfo(messedUpDur(mud)).PokeInIndex)/obj.sampleRate;
                    elseif obj.trialInfo(messedUpDur(mud)).Performance == 1
                        obj.trialInfo(messedUpDur(mud)).PokeOutIndex = obj.trialInfo(messedUpDur(mud)).RewardSignalIndex;
                        newDur = (obj.trialInfo(messedUpDur(mud)).RewardSignalIndex - obj.trialInfo(messedUpDur(mud)).RewardSignalIndex)/obj.sampleRate;
                    end
                    tStart(messedUpDur(mud)) = newDur*-1;
                    startNdx{messedUpDur(mud)} = find(newDur<obj.obsvTimeVect,1,'first');
                end
                startNdx = cell2mat(startNdx);
                tEnd = zeros(size(obj.trialInfo));
                endNdx = repmat(find(obj.obsvTimeVect==0), size(obj.trialInfo));
            end
            trlStarts = [obj.trialInfo.PokeInIndex];
            trlEnds = [obj.trialInfo.PokeOutIndex];
            pos1log = [obj.trialInfo.Position]==1;
            pos1ndx = find(pos1log);
            preTrlIntDur = [nan, trlStarts(2:end)-trlEnds(1:end-1)]/obj.sampleRate;            
            pstTrlIntDur = [trlStarts(2:end) - trlEnds(1:end-1), nan]/obj.sampleRate;

            meanNonP1Int = mean(preTrlIntDur(~pos1log));
            preTrlIntDur(pos1log) = meanNonP1Int;
            prevTrlStartNdx = arrayfun(@(a)find(a<obj.obsvTimeVect, 1, 'first'), preTrlIntDur*-1);
            preTrlIntDurMids = tStart-(preTrlIntDur./2);
            preTrlIntDurMidNdx = arrayfun(@(a)find(a<obj.obsvTimeVect,1,'first'), preTrlIntDurMids);
            pstTrlIntDur([pos1ndx(2:end)-1,end]) = meanNonP1Int;
            nextTrlStartNdx = arrayfun(@(a)find(a>obj.obsvTimeVect, 1, 'last'), pstTrlIntDur);
            pstTrlIntDurMids = tEnd+(pstTrlIntDur./2);
            pstTrlIntDurMidNdx = arrayfun(@(a)find(a>obj.obsvTimeVect, 1, 'last'), pstTrlIntDurMids);

            obj.mask_Trial = obj.BuildVectMaskMtx(startNdx,[startNdx', endNdx']-startNdx',length(obj.obsvTimeVect));
            if strcmp(obj.alignments{1}, 'PokeIn')
                obj.mask_Int = obj.BuildVectMaskMtx(startNdx, [prevTrlStartNdx', startNdx']-startNdx', length(obj.obsvTimeVect));
                obj.mask_IntMid = obj.BuildVectMaskMtx(startNdx, repmat(preTrlIntDurMidNdx', [1,2]), length(obj.obsvTimeVect));
            elseif strcmp(obj.alignments{1}, 'PokeOut')
                obj.mask_Int = obj.BuildVectMaskMtx(startNdx, [startNdx', nextTrlStartNdx]-startNdx', length(obj.obsvTimeVect));
                obj.mask_IntMid = obj.BuildVectMaskMtx(startNdx, repmat(pstTrlIntDurMidNdx', [1,2]), length(obj.obsvTimeVect));
            else
                error('What, poke in and poke out alignments aren''t good enough for you you? Get outta here and code if yourself ya lazy bum!' )
            end
            obj.mask_PreInt = obj.BuildVectMaskMtx(startNdx,[preTrlIntDurMidNdx', startNdx']-startNdx',length(obj.obsvTimeVect));
            obj.mask_PreIntCON = obj.BuildVectMaskMtx(startNdx,[preTrlIntDurMidNdx', startNdx']-startNdx'-(floor((startNdx-preTrlIntDurMidNdx)/2))',length(obj.obsvTimeVect));
            obj.mask_PstInt = obj.BuildVectMaskMtx(startNdx,[endNdx' pstTrlIntDurMidNdx']-startNdx',length(obj.obsvTimeVect));
            obj.mask_PstIntCON = obj.BuildVectMaskMtx(startNdx,[endNdx' pstTrlIntDurMidNdx']-startNdx'+(floor((pstTrlIntDurMidNdx-endNdx)/2))',length(obj.obsvTimeVect));
        end
    end
    %% Analyses
    methods
        %% Calc Model Persistence Fits
        function [trlFits, intFits] = CalcModelPersistenceFit_XTD(obj,data,trialInfo, varargin)
             if nargin==1 || isempty(data)
                [data, ~, trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
             elseif nargin==2
                 [~,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
             end
             if isempty(trialInfo)
                 [~,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
             end
             if ismatrix(data)
                 [data,~,~] = obj.OrganizeDecodabilityTrialHistoryTransMat(data);
             end
             trlFits = cell(obj.seqLength, obj.seqLength, obj.seqLength);
             intFits = cell(obj.seqLength, obj.seqLength, obj.seqLength);
             for odrPos = 1:obj.seqLength
                 for histOSpos = 1:obj.seqLength
                     for trlPos = 1:obj.seqLength
                         tempData = data{odrPos,histOSpos,trlPos,trlPos};
                         tempTrialInfo = trialInfo{odrPos,histOSpos,trlPos};

                         temp_trlFits = nan(size(tempData,3), length(obj.obsvTimeVect));
                         temp_intFits = nan(size(tempData,3), length(obj.obsvTimeVect));
                         if ~isempty(tempData)
                             for trl = 1:size(tempData,3)
                                 cur_tempData = tempData(:,:,trl);
                                 
                                 curTrialMask = obj.mask_Trial(:,tempTrialInfo(trl).TrialNum);
                                 trialDynCap = sum(curTrialMask)-1;
                                 trialPeriod = zscore(cur_tempData(curTrialMask,curTrialMask),0,'all');
                                 for d = 1:trialDynCap-1
                                     tempDyn = zscore(triu(true(sum(curTrialMask)), d*-1) & tril(true(sum(curTrialMask)), d), 0, 'all');
                                     temp_trlFits(trl,d) = pdist([trialPeriod(:)';tempDyn(:)'],'cosine');
                                 end

                                 curIntMask = obj.mask_Int(:,tempTrialInfo(trl).TrialNum);
                                 intDynCap = sum(curIntMask)-1;
                                 intPeriod = zscore(cur_tempData(curIntMask,curIntMask),0,'all');
                                 for d = 1:intDynCap-1
                                     tempDyn = zscore(triu(true(sum(curIntMask)), d*-1) & tril(true(sum(curIntMask)), d), 0, 'all');
                                     temp_intFits(trl,d) = pdist([intPeriod(:)';tempDyn(:)'],'cosine');
                                 end
                             end
                         end
                         trlFits{odrPos,histOSpos,trlPos} = temp_trlFits;
                         intFits{odrPos,histOSpos,trlPos} = temp_intFits;
                     end
                 end
             end
             if strcmp(varargin, 'latency')>=1
                 tempTrialFits = cell(size(trlFits));
                 tempIntFits = cell(size(intFits));
                 for p = 1:numel(tempTrialFits)
                     if ~isempty(trlFits{p})
                         curTrlFits = trlFits{p};
                         [~,tempTrialFitLats]=find(curTrlFits==repmat(min(curTrlFits,[],2),[1,size(curTrlFits,2)]));
                         tempTrialFits{p} = obj.obsvTimeVect(tempTrialFitLats);
                         curIntFIts = intFits{p};                         
                         [~,tempIntFitLats]=find(curIntFIts==repmat(min(curIntFIts,[],2),[1,size(curIntFIts,2)]));
                         tempIntFits{p} = obj.obsvTimeVect(tempIntFitLats);
                     end
                 end
                 trlFits = tempTrialFits;
                 intFits = tempIntFits;
             end
        end
        %% Calc Decoding Peaks
%         function [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodabilityPeaks_XTD(obj,data,dataType)
%             if nargin<2
%                 [data, ~, ~] = obj.OrganizeDecodabilityTrialHistoryTransMat;
%             elseif nargin==2 && ismatrix(data)
%                 [data,~,~] = obj.OrganizeDecodabilityTrialHistoryTransMat(data);
%             elseif nargin==3 
%                 if isempty(data) && ~isempty(dataType)
%                     [hr,d,~] = obj.OrganizeDecodabilityTrialHistoryTransMat;
%                     if strcmp(dataType, 'd') || strcmp(dataType, 'D')
%                         data = d;
%                     else
%                         data = hr;
%                     end
%                 else
%                     [data, ~, ~] = obj.OrganizeDecodabilityTrialHistoryTransMat;
%                 end
%             end
%             tm_PeakNdx = cell(obj.seqLength, obj.seqLength, obj.seqLength);
%             tm_PeakVal = cell(obj.seqLength, obj.seqLength, obj.seqLength);
%             tm_PeakWid = cell(obj.seqLength, obj.seqLength, obj.seqLength);
%             for odrPos = 1:obj.seqLength
%                 for histOSpos = 1:obj.seqLength
%                     for trlPos = 1:obj.seqLength
%                         temp_Data = squeeze(data(odrPos,histOSpos,trlPos,:));
%                         if ~isempty(temp_Data)
%                             temp_PeakNdx = nan(obj.seqLength,size(temp_Data{1},2),size(temp_Data{1},3));
%                             temp_PeakVal = nan(obj.seqLength,size(temp_Data{1},2),size(temp_Data{1},3));
%                             temp_PeakWid = nan(obj.seqLength,size(temp_Data{1},2),size(temp_Data{1},3));
%                             for pos = 1:size(temp_Data,1)
%                                 cur_temp_Data = temp_Data{pos};
%                                 for t = 1:size(cur_temp_Data,2)
%                                     for trl = 1:size(cur_temp_Data,3)
%                                         [pks,loc,wid,prom] = findpeaks(cur_temp_Data(:,t),'minpeakdistance', obj.binSize/obj.dsRate);
%                                         if ~isempty(pks)
%                                             featWeightPeaks = pks.*wid.*prom;
%                                             temp_PeakNdx(pos,t,trl) = obj.obsvTimeVect(loc(featWeightPeaks==max(featWeightPeaks)));
%                                             temp_PeakVal(pos,t,trl) = pks(featWeightPeaks==max(featWeightPeaks));
%                                             temp_PeakWid(pos,t,trl) = wid(featWeightPeaks==max(featWeightPeaks));
%                                         end
%                                     end
%                                 end
%                             end
%                             tm_PeakNdx{odrPos,histOSpos,trlPos} = temp_PeakNdx;
%                             tm_PeakVal{odrPos,histOSpos,trlPos} = temp_PeakVal;
%                             tm_PeakWid{odrPos,histOSpos,trlPos} = temp_PeakWid;
%                         end
%                     end
%                 end                
%             end
%         end
        %% Calc Windowed Decoding Peaks
        function [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodablityPeaks_XTD_Windowed(obj,data,window)
            if nargin==1 || isempty(data)
                [data, ~, ~] = obj.OrganizeDecodabilityTrialHistoryTransMat;
                window = [min(obj.obsvTimeVect), max(obj.obsvTimeVect)];
            end
            if max(abs(window))>obj.sampleRate
                window = window./obj.sampleRate;
            end
            if ismatrix(data)
                [data,~,~] = obj.OrganizeDecodabilityTrialHistoryTransMat(data);
            end
            windowLog = obj.obsvTimeVect>=window(1) & obj.obsvTimeVect<=window(2);
            tempData = cell(size(data));
            for d = 1:numel(data)
                if ~isempty(data{d})
                    tempData{d} = obj.MaskArrayExtract(data{d},windowLog, windowLog);
                end
            end
            [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodabilityPeaks_XTD(tempData);
            
        end
        %% Calc Masked Decoding Peaks
        function [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodabilityPeaks_XTD_Mask(obj,data,trialInfo,mask1,mask2)
            if nargin==1 || isempty(data)
                [data,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
                mask1 = obj.mask_Trial;
                mask2 = obj.mask_Trial;
            end
            if isempty(trialInfo)
                [~,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
            end
            if nargin==3
                mask1 = obj.mask_Trial;
                mask2 = obj.mask_Trial;
            elseif nargin==4
                mask2 = mask1;
            end
            tempData = cell(size(data));
            
            for odrPos = 1:obj.seqLength
                for histOSpos = 1:obj.seqLength
                    for trlPos = 1:obj.seqLength
                        temp_tempData = squeeze(data(odrPos,histOSpos,trlPos,:));
                        if ~isempty(temp_tempData{1})
                            temp_trialIDs = [trialInfo{odrPos,histOSpos,trlPos}.TrialNum];
                            for p = 1:length(temp_tempData)
                                tempData{odrPos,histOSpos,trlPos,p} = obj.MaskArrayExtract(temp_tempData{p},mask1(:,temp_trialIDs),mask2(:,temp_trialIDs));
                            end
                        end
                    end
                end
            end
            [tm_PeakNdx, tm_PeakVal, tm_PeakWid] = CalcDecodabilityPeaks_XTD(obj,tempData);
        end
        %% Extract Exemplar XTD Models
        function [intDecMod, preDecMod, pstDecMod, offsets] = ExtractExemplarXTDmodels(obj,data,dataType,trialInfo,varargin)
            if nargin==1 || isempty(data)
                [data, ~, trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
            elseif nargin==2 || ismatrix(data)
                [data,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat(data);
            elseif nargin==3 && isempty(data)
                [hr,d,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
                if ~isempty(dataType) && (strcmp(dataType, 'd') || strcmp(dataType, 'D'))
                    data = d;
                else
                    data = hr;
                end            
            end
            if sum(strcmp(varargin, 'aniOffset'))>=1
                oscLog = [obj.trialInfo.Performance] == 1 & [obj.trialInfo.TranspositionDistance]~=0;
                decLat = mean([obj.trialInfo(oscLog).PokeDuration]);
                decLatOffset = floor(decLat*obj.sampleRate/obj.dsRate);
                preDecOffset = floor(decLatOffset/2);
                pstDecOffset = floor(decLatOffset*1.5);
            else
                preDecOffset = 9;
                pstDecOffset = 17;
            end
            if sum(strcmp(varargin,'windowed'))>=1
                tempWindow = varargin{find(strcmp(varargin,'windowed'))+1};
                if max(abs(tempWindow))>obj.sampleRate
                    tempWindow = tempWindow./obj.sampleRate;
                end
                trlLog = obj.obsvTimeVect>tempWindow(1) & obj.obsvTimeVect<tempWindow(2);
            else
                trlLog = true(size(obj.obsvTimeVect));
            end
            offsets = [preDecOffset, pstDecOffset];

            intDecMod = cell(obj.seqLength,obj.seqLength,obj.seqLength);
            preDecMod = cell(obj.seqLength,obj.seqLength,obj.seqLength);
            pstDecMod = cell(obj.seqLength,obj.seqLength,obj.seqLength);

            for odrPos = 1:obj.seqLength
                for histOSpos = 1:obj.seqLength
                    for trlPos = 1:obj.seqLength
                        if ~isempty(data{odrPos,histOSpos,trlPos,1})
                            cur_tempData = squeeze(data(odrPos,histOSpos,trlPos,:));
                            temp_intDecMod = nan(size(obj.mask_IntMid,1),length(cur_tempData),size(cur_tempData{1},3));
                            temp_preDecMod = nan(size(obj.mask_Trial,1),length(cur_tempData),size(cur_tempData{1},3));
                            temp_pstDecMod = nan(size(obj.mask_Trial,1),length(cur_tempData),size(cur_tempData{1},3));
                            temp_trialIDs = [trialInfo{odrPos,histOSpos,trlPos}.TrialNum];
                            for p = 1:length(cur_tempData)
                                temp_intDecMod(trlLog,p,:) = obj.MaskArrayExtract(cur_tempData{p},trlLog,obj.mask_IntMid(:,temp_trialIDs),'extract');
                                temp_preDecMod(trlLog,p,:) = obj.MaskArrayExtract(cur_tempData{p},trlLog,[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',preDecOffset,preDecOffset)}], 'extract');
                                temp_pstDecMod(trlLog,p,:) = obj.MaskArrayExtract(cur_tempData{p},trlLog,[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',pstDecOffset,pstDecOffset)}], 'extract');
                            end
                            intDecMod{odrPos,histOSpos,trlPos} = temp_intDecMod;
                            preDecMod{odrPos,histOSpos,trlPos} = temp_preDecMod;
                            pstDecMod{odrPos,histOSpos,trlPos} = temp_pstDecMod;
                        end
                    end
                end
            end
        end
        %% Extract Exemplar XTD Model Peaks
        function [intNdx, intVals, intWids, preDecNdx, preDecVals, preDecWids,pstDecNdx, pstDecVals, pstDecWids] = CalcDecodabilityPeaks_ExemplarXTDmodels(obj,data,dataType,trialInfo,window,varargin)
             if nargin==1 || isempty(data)
                [data, ~, trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
            elseif nargin==2 || ismatrix(data)
                [data,~,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat(data);
            elseif nargin==3 && isempty(data)
                [hr,d,trialInfo] = obj.OrganizeDecodabilityTrialHistoryTransMat;
                if ~isempty(dataType) && (strcmp(dataType, 'd') || strcmp(dataType, 'D'))
                    data = d;
                else
                    data = hr;
                end            
             end
             if sum(strcmp(varargin, 'aniOffset'))>=1
                oscLog = [obj.trialInfo.Performance] == 1 & [obj.trialInfo.TranspositionDistance]~=0;
                decLat = mean([obj.trialInfo(oscLog).PokeDuration]);
                decLatOffset = floor(decLat*obj.sampleRate/obj.dsRate);
                preDecOffset = floor(decLatOffset/2);
                pstDecOffset = floor(decLatOffset*1.5);
            else
                preDecOffset = 9;
                pstDecOffset = 17;
            end
            [tempPeakNdx, tempPeakVals, tempPeakWids] = obj.CalcDecodablityPeaks_XTD_Windowed(data,window);
            intNdx = cell(size(tempPeakNdx));
            intVals = cell(size(tempPeakVals));
            intWids = cell(size(tempPeakWids));
            preDecNdx = cell(size(tempPeakNdx));
            preDecVals = cell(size(tempPeakVals));
            preDecWids = cell(size(tempPeakWids));
            pstDecNdx = cell(size(tempPeakNdx));
            pstDecVals = cell(size(tempPeakVals));
            pstDecWids = cell(size(tempPeakWids));
            for p = 1:numel(tempPeakNdx)
                if ~isempty(tempPeakNdx{p})
                    temp_trialIDs = [trialInfo{p}.TrialNum];
                    intNdx{p} = obj.MaskVectExtract(tempPeakNdx{p},obj.mask_IntMid(:,temp_trialIDs),2,3);
                    intVals{p} = obj.MaskVectExtract(tempPeakVals{p},obj.mask_IntMid(:,temp_trialIDs),2,3);
                    intWids{p} = obj.MaskVectExtract(tempPeakWids{p},obj.mask_IntMid(:,temp_trialIDs),2,3);
                    preDecNdx{p} = obj.MaskVectExtract(tempPeakNdx{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',preDecOffset,preDecOffset)}],2,3);
                    preDecVals{p} = obj.MaskVectExtract(tempPeakVals{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',preDecOffset,preDecOffset)}],2,3);
                    preDecWids{p} = obj.MaskVectExtract(tempPeakWids{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',preDecOffset,preDecOffset)}],2,3);
                    pstDecNdx{p} = obj.MaskVectExtract(tempPeakNdx{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',pstDecOffset,pstDecOffset)}],2,3);
                    pstDecVals{p} = obj.MaskVectExtract(tempPeakVals{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',pstDecOffset,pstDecOffset)}],2,3);
                    pstDecWids{p} = obj.MaskVectExtract(tempPeakWids{p},[{obj.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',pstDecOffset,pstDecOffset)}],2,3);
                end
            end
        end              
    end
    %% Visualizations
    methods
        %% Trial-wise probability density matrix
        function PlotTrialPDF(obj, realTrialMtx, varargin)
            if sum(strcmp(varargin, 'x') | strcmp(varargin, 'X'))>=1
                x = varargin{find(strcmp(varargin, 'x') | strcmp(varargin, 'X'))+1};
            else
                x = 1:size(realTrialMtx,2);
            end
            if sum(strcmp(varargin, 'y') | strcmp(varargin, 'Y'))>=1
                y = varargin{find(strcmp(varargin, 'y') | strcmp(varargin, 'Y'))+1};
            else
                y = 1:size(realTrialMtx,1);
            end
%             realMean = median(realTrialMtx,3,'omitnan');         
            realMean = mean(realTrialMtx,3,'omitnan');
%             realMean = median(realTrialMtx,3,'omitnan')./std(realTrialMtx,0,3);
%             realMean = mean(realTrialMtx,3,'omitnan')./std(realTrialMtx,0,3);
            if sum(strcmp(varargin, 'rotate') | strcmp(varargin, 'Rotate'))>=1
                imagesc(x, y, realMean');
            else
                imagesc(x, y, realMean);
            end
            if sum(strcmp(varargin, 'xlabel') | strcmp(varargin, 'xLabel') | strcmp(varargin, 'Xlabel') | strcmp(varargin, 'XLabel'))>=1
                xlabel(varargin{find(strcmp(varargin, 'xlabel') | strcmp(varargin, 'xLabel') | strcmp(varargin, 'Xlabel') | strcmp(varargin, 'XLabel'))+1});
            end
            if sum(strcmp(varargin, 'ylabel') | strcmp(varargin, 'yLabel') | strcmp(varargin, 'Ylabel') | strcmp(varargin, 'YLabel'))>=1
                ylabel(varargin{find(strcmp(varargin, 'ylabel') | strcmp(varargin, 'yLabel') | strcmp(varargin, 'Ylabel') | strcmp(varargin, 'YLabel'))+1});
            end
            colorbar
            set(gca, 'ydir', 'normal');
            hold on;
            if sum(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))
                if ~ischar(varargin{find(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))+1})
                    set(gca, 'clim', varargin{find(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))+1});
                elseif strcmp(varargin{find(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))+1}, 'Min')
                    if min(get(gca, 'clim'))<0
                        set(gca, 'clim', [min(get(gca, 'clim')), min(get(gca, 'clim'))*-1]);
                    else
                        set(gca, 'clim', [max(abs(get(gca, 'clim')))*-1, max(abs(get(gca, 'clim')))]);
                    end                        
                elseif strcmp(varargin{find(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))+1}, 'Max')
                    if max(get(gca, 'clim'))>0
                        set(gca, 'clim', [max(get(gca, 'clim'))*-1, max(get(gca, 'clim'))]);
                    else                        
                        set(gca, 'clim', [max(abs(get(gca, 'clim')))*-1, max(abs(get(gca, 'clim')))]);
                    end
                elseif strcmp(varargin{find(strcmp(varargin, 'clim') | strcmp(varargin, 'cLim'))+1}, 'Abs')
                    set(gca, 'clim', [max(abs(get(gca, 'clim')))*-1, max(abs(get(gca, 'clim')))]);
                end
            end
                    
            if sum(strcmp(varargin, 'thresh') | strcmp(varargin, 'Thresh'))>=1
                thresh = varargin{find(strcmp(varargin, 'thresh') | strcmp(varargin, 'Thresh'))+1};
            else
                thresh = 0.975;
            end
            if sum(strcmp(varargin, 'chance') | strcmp(varargin, 'Chance'))>=1
                chanceTrialMtx = varargin{find(strcmp(varargin, 'chance') | strcmp(varargin, 'Chance'))+1};
                realSEM = obj.SEMcalc(realTrialMtx,0,3);
                realCI = tinv(thresh, size(realTrialMtx,3)-1).*realSEM;
                tempDthresh = chanceTrialMtx(:,:,1)+(tinv(thresh, numChancePerms-1).*(chanceTrialMtx(:,:,2)./sqrt(numChancePerms-1)));
                abvThresh = (realMean-realCI)>tempDthresh;
                bounds = bwBoundaries(abvThresh);
                for b = 1:length(bounds)
                    if numel(bounds{b})>4
                        plot(obj.obsvTimeVect(bounds{b}(:,1)), mlb.obsvTimeVect(bounds{b}(:,2)), 'k', 'linewidth', 2);
                    end
                end                        
            end
        end            
        %% Sequence-wise probability density matrix
        function PlotSequencePDF(obj)
        end
        %% Line Plots w/Chance
        function PlotLinesWithChance(obj)
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
    end
    %% Analyses
    methods
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
    end
    methods % PCA stuff... prototyped.... may not currently work...
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
%                 curTrlInfo = obj.trialInfo(obj.likeIDvects{perm});
                                
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