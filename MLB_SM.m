classdef MLB_SM < SeqMem
    % Implementation of a naive bayes classifier assuming a poisson distribution for spiking data
    properties % Analysis Variables
        binSize
        dsRate
        binType = 'box'
        numPerms = 10
        ssProportion = 0.5
        ssType % 0 = use all ISC for decoding; 1 = use subsampled ISC types
        bayesType %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts
        windows
        alignments
    end
    properties % Data Structures
        fiscTrials        
        trialIDs = [{'Time'}, {'Window'}, {'Position'}, {'Odor'}, {'TrialNum'}];
        likeTrlSpikes
        likeIDvects
        obsvTrlSpikes
        obsvIDvects
        post
        postIDvects
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
            obj.fiscTrials = potentialSeqs(:,fiscLog);
        end
        %% Group Trial Data
        function [ssnSpikes, ssnID] = PP_ConcatTrialData(obj)
            trlDta = cell(size(obj.windows,1),1);
            trlTimeVects = cell(size(obj.windows,1),1);
            for win = 1:size(obj.windows,1)
                [trlDta{win}, trlTimeVects{win}] = obj.PP_TrialMatrix_Spiking(obj.windows{win}, obj.alignments{win});
            end
            ssnSpikes = nan(size(cell2mat(trlTimeVects),1), size(obj.ensembleMatrix,2), length(obj.trialInfo));
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
                ssnID(:,1,trl) = cell2mat(tempTrialTimeVect);
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
            % Set likelihoods & observations using FISC trials where FISC are likelihoods and all other ISC are observations
            if size(obj.odrSeqs,1)>=2
                warning('FISC likes should only be used with single list sessions (use sub-sampling for dual list). Ignore if ONLY examining ordinal coding, or if task conditions changed and code was updated but warning statement wasn''t (check github/readme).');
            end
            [ssnSpikes, ssnID] = obj.PP_ConcatTrialData;
            obj.likeTrlSpikes = {nan(size(ssnSpikes,1)*size(obj.fiscTrials,1), size(ssnSpikes,2), size(obj.fiscTrials,2))};
            obj.likeIDvects = {nan(size(ssnID,1)*size(obj.fiscTrials,1), size(ssnID,2), size(obj.fiscTrials,2))};
            for seq = 1:size(obj.fiscTrials,2)
                tempSeqData = cell(size(obj.fiscTrials,1),1);
                tempSeqIDlog = cell(size(obj.fiscTrials,1),1);
                for pos = 1:size(obj.fiscTrials,1)
                    tempSeqData{pos} = ssnSpikes(:,:,obj.fiscTrials(pos,seq));
                    tempSeqIDlog{pos} = ssnID(:,:,obj.fiscTrials(pos,seq));
                end
                obj.likeTrlSpikes{1}(:,:,seq) = cell2mat(tempSeqData);
                obj.likeIDvects{1}(:,:,seq) = cell2mat(tempSeqIDlog);
            end
            iscLog = ([obj.trialInfo.TranspositionDistance]==0 & [obj.trialInfo.Performance]==1);
            iscLog(obj.fiscTrials) = false;
            obj.obsvTrlSpikes = {ssnSpikes(:,:,iscLog)};
            obj.obsvIDvects = {ssnID(:,:,iscLog)};
        end
        %% Set Likelihoods as SubSampled ISC
        function SetLikes_SubSample(obj)
            if isempty(obj.ssType)
                error('Specify sub-sampling method for selecting observations');
            end
            % Set likelihoods & observations using sets of subsampled ISC trials
            [ssnSpikes, ssnID] = obj.PP_ConcatTrialData;
            % Determine number of trials used for subsampled likelihoods
            likeliSize = floor(min(min(obj.isTrialNums(:,:,1)))*obj.ssProportion);
            trlCounts = reshape(obj.isTrialNums(:,:,1)',[numel(obj.isTrialNums(:,:,1)),1]);
            odrIDs = reshape(obj.odrSeqs',[numel(obj.odrSeqs),1]);
            likeTrials = nan(length(trlCounts), likeliSize, obj.numPerms);
            for oip = 1:length(trlCounts)
                tempPool = randi(trlCounts(oip), [likeliSize,obj.numPerms]);
                for r = 1:obj.numPerms
                    while length(unique(tempPool(:,r)))~=size(tempPool,1)
                        tempPool(:,r) = [unique(tempPool(:,r)); randi(trlCounts(oip), [size(tempPool,1)-length(unique(tempPool(:,r))),1])];
                    end
                end
                while size(unique(tempPool', 'rows'),1)~=obj.numPerms
                    tempPool = [unique(tempPool', 'rows')', randi(trlCounts(oip), [likeliSize,obj.numPerms-size(unique(tempPool', 'rows'),1)])];
                    for r = 1:obj.numPerms
                        while length(unique(tempPool(:,r)))~=size(tempPool,1)
                            tempPool(:,r) = [unique(tempPool(:,r)); randi(trlCounts(oip), [size(tempPool,1)-length(unique(tempPool(:,r))),1])];
                        end
                    end
                end
                if size(unique(tempPool', 'rows'),1)~=obj.numPerms
                    error('Apparently I need to code for this eventuality too!?');
                end
                likeTrials(oip,:,:) = tempPool;
            end
            odrTrlSpikes = cell(size(odrIDs));
            odrTrlIDs = cell(size(odrIDs));
            for odr = 1:length(odrIDs)
                odrIDlog = [obj.trialInfo.Odor]==odrIDs(odr);
                odrPosLog = [obj.trialInfo.Position]==find(sum(obj.odrSeqs==odrIDs(odr),1));
                odrTrlSpikes{odr} = ssnSpikes(:,:,(odrIDlog & odrPosLog & [obj.trialInfo.Performance]==1));
                odrTrlIDs{odr} = ssnID(:,:,(odrIDlog & odrPosLog & [obj.trialInfo.Performance]==1));
            end
            obj.likeTrlSpikes = repmat({nan(size(ssnSpikes,1)*size(likeTrials,1), size(ssnSpikes,2), size(likeTrials,2))}, [1,obj.numPerms]);
            obj.obsvTrlSpikes = cell(1,obj.numPerms);
            obj.likeIDvects = repmat({nan(size(ssnID,1)*size(likeTrials,1), size(ssnID,2), size(likeTrials,2))}, [1,obj.numPerms]);
            obj.obsvIDvects = cell(1,obj.numPerms);
            for perm = 1:obj.numPerms
                tempOdrTrlSpikes = odrTrlSpikes;
                tempOdrTrlIDs = odrTrlIDs;
                for seq = 1:size(likeTrials,2)
                    tempLikeSpikes = cell(length(odrIDs),1);
                    tempLikeIDs = cell(length(odrIDs),1);
                    for odr = 1:size(likeTrials,1)
                        tempLikeSpikes{odr} = tempOdrTrlSpikes{odr}(:,:,likeTrials(odr,seq));
                        tempOdrTrlSpikes{odr}(:,:,likeTrials(odr,seq)) = nan;
                        tempLikeIDs{odr} = tempOdrTrlIDs{odr}(:,:,likeTrials(odr,seq));
                        tempOdrTrlIDs{odr}(:,:,likeTrials(odr,seq)) = nan;
                    end               
                    obj.likeTrlSpikes{perm}(:,:,seq) = cell2mat(tempLikeSpikes);
                    obj.likeIDvects{perm}(:,:,seq) = cell2mat(tempLikeIDs);
                end
                if obj.ssType == 1
                    for odr = 1:size(tempOdrTrlSpikes,1)
                        tempPool = randi(size(tempOdrTrlSpikes{odr},3), [1,likeliSize]);
                        tempPool = unique(tempPool);
                        tempPool(reshape(isnan(tempOdrTrlSpikes{odr}(1,1,tempPool)), [1,length(tempPool)])) = [];
                        while length(tempPool)~=likeliSize
                            tempPool = [tempPool, randi(size(tempOdrTrlSpikes{odr},3), [1, likeliSize-length(tempPool)])]; %#ok<AGROW>
                            tempPool = unique(tempPool);
                            tempPool(reshape(isnan(tempOdrTrlSpikes{odr}(1,1,tempPool)), [1,length(tempPool)])) = [];
                        end
                        tempOdrTrlSpikes{odr} = tempOdrTrlSpikes{odr}(:,:,tempPool);
                        tempOdrTrlIDs{odr} = tempOdrTrlIDs{odr}(:,:,tempPool);
                    end
                end                    
                tempObsvSpikes = cell2mat(reshape(tempOdrTrlSpikes, [1,1,numel(tempOdrTrlSpikes)]));
                tempObsvIDs = cell2mat(reshape(tempOdrTrlIDs, [1,1,numel(tempOdrTrlIDs)]));
                obsvLog = sum(sum(isnan(tempObsvSpikes)))==0;
                obj.obsvTrlSpikes{perm} = tempObsvSpikes(:,:,obsvLog);
                obj.obsvIDvects{perm} = tempObsvIDs(:,:,obsvLog);
            end
        end
    end
    methods % MLB Processing Methods
        %% Process via Leave-1-Out
        function Process_LikelyLOO(obj)
            obj.post = repmat({nan(size(obj.likeTrlSpikes{1},1), size(obj.likeTrlSpikes{1},1), size(obj.likeTrlSpikes{1},3))}, size(obj.likeTrlSpikes));
            obj.postIDvects = repmat({nan(size(obj.likeIDvects{1}))}, size(obj.likeIDvects));
            for perm = 1:length(obj.likeTrlSpikes)
                for seq = 1:size(obj.likeTrlSpikes{perm},3)
                    tempObsv = obj.likeTrlSpikes{perm}(:,:,seq);
                    tempLike = obj.likeTrlSpikes{perm};
                    tempLike(:,:,seq) = [];
                    if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                        obj.post{perm}(:,:,seq) = obj.CalcStaticBayesPost_Poisson(mean(tempLike,3), tempObsv);
                    elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                        obj.post{perm}(:,:,seq) = obj.CalcStaticBayesPost_Bernoulli(mean(tempLike,3), tempObsv);
                    elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                        obj.post{perm}(:,:,seq) = CalcStaticBayesPost_Gaussian(mean(tempLike,3), std(tempLike,0,3), tempObsv);
                    end
                    obj.postIDvects{perm}(:,:,seq) = obj.likeIDvects{perm}(:,:,seq);
                end
            end
        end
        %% Process all Observations
        function Process_Observes(obj)
            obj.post = cell(size(obj.obsvTrlSpikes));
            obj.postIDvects = cell(size(obj.obsvTrlSpikes));
            for perm = 1:length(obj.likeTrlSpikes)
                if obj.bayesType == 1 || strcmp(obj.bayesType, 'Poisson') || strcmp(obj.bayesType, 'poisson') || strcmp(obj.bayesType, 'P') || strcmp(obj.bayesType, 'p')
                    obj.post{perm} = obj.CalcStaticBayesPost_Poisson(mean(obj.likeTrlSpikes{perm},3), obj.obsvTrlSpikes{perm});
                elseif obj.bayesType == 2 || strcmp(obj.bayesType, 'Bernoulli') || strcmp(obj.bayesType, 'bernoulli') || strcmp(obj.bayesType, 'B') || strcmp(obj.bayesType, 'b')
                    obj.post{perm} = obj.CalcStaticBayesPost_Bernoulli(mean(obj.likeTrlSpikes{perm},3), obj.obsvTrlSpikes{perm});
                elseif obj.bayesType == 3 || strcmp(obj.bayesType, 'Gaussian') || strcmp(obj.bayesType, 'gaussian') || strcmp(obj.bayesType, 'G') || strcmp(obj.bayesType, 'g')
                    obj.post{perm} = obj.CalcStaticBayesPost_Gaussian(mean(obj.likeTrlSpikes{perm},3), std(obj.likeTrlSpikes{perm},0,3), obj.obsvTrlSpikes{perm});
                end
                obj.postIDvects{perm} = obj.obsvTrlSpikes{perm};
            end
        end
    end
    methods % MLB Algorithms
        %% Calculate Static MLB (Poisson... kept in here just so all previous analyses won't break... will remove eventually as the rest of the code base gets revised).
        function post = CalcStaticBayesPost(obj,likely, obsv)
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
                    tempPost = pp.*ee;
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end
            %             toc
        end
        %% Calculate Static MLB Poisson
        function post = CalcStaticBayesPost_Poisson(obj,likely, obsv)
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
                    tempPost = pp.*ee;
                    post(t,:,trl) = tempPost./sum(tempPost);
                end
            end
            %             toc
        end
        %% Calculate Static MLB Bernoulli
        function post = CalcStaticBayesPost_Bernoulli(~,likely, obsv)
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
                    indiv_neurons = nan(size(obsv,1), size(likely,1));
                    for u = 1:size(likely,2)
                        indiv_neurons(:,u) = (likely(:,u).^obsv(t,u,trl)).*((1-likely(:,u)).^(1-obsv(t,u,trl)));
                    end
                    post(t,:,trl) = prod(indiv_neurons,2);
                end
            end
        end
        %% Calculate Static MLB Gaussian
        function post = CalcStaticBayesPost_Gaussian(~,meanLikely, varLikely, obsv)
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
            post = nan(size(obsv,1), size(meanLikely,1), size(obsv,3));
            for trl = 1:size(obsv,3)
                for t = 1:size(obsv,1)
                    first_term = prod(1./(varLikely.*((2*pi)^0.5)),2);
                    second_term = nan(size(likely));
                    for u = 1:size(likely,2)
                        second_term(:,u) = -0.5.*(((obsv(t,u)-meanLikely(:,u))./varLikely(:,u)).^2);
                    end
                    second_term = exp(sum(second_term,2));
                    post(t,:,trl) = first_term.*second_term;
                end
            end                        
        end
        %% Calculate Iterative MLB Poisson
        % **** TO CREATE ****
        function post = CalcIterativeBayesPost_Poisson(obj, likely, obsv, depVar, grpVar)
            % Calculate posterior probability using a dependant variable and a grouping variable
            %   The dependant variable (mainly time) is the variable that is iterated across
            %   The grouping variable (mainly odor or position) is the variable used as the likelihoods at each level of the dependant variable
            % For example, typical use case here is to look for temporally invariant coding where
            %   depVar = time, i.e. the algorithm will step through different time points
            %   grpVar = odor/position, i.e. the algorithm will then choose likelihoods from different levels of odor or position
            %             tic;
            rateScalar = obj.binSize/obj.sampleRate;
            lvlsDepVar = unique(depVar);
            lvlsGrpVar = unique(grpVar);
            post = nan(size(obsv,1), size(likely,1)/length(lvlsGrpVar), length(lvlsGrpVar), size(obsv,3));
            for dv = 1:length(lvlsDepVar)
                tempLikely = nan(length(lvlsGrpVar), size(likely,2));
                for gv = 1:length(lvlsGrpVar)
                    tempLikely(gv,:) = likely(depVar==lvlsDepVar(dv) & grpVar==lvlsGrpVar(gv),:);
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
                        tempTempPost = pp.*ee;
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
        %... don't really remember what this's doing... or why I made it...
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