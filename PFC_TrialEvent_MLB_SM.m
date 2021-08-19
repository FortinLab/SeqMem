classdef PFC_TrialEvent_MLB_SM < MLB_SM
    %% Properties
    properties % Parameters
        beginTrialAlignment = 'PokeIn'
        endTrialAlignment = 'PokeOut'
        beginTrialWindow = [-500 500]
        endTrialWindow = [-500 500]
        thetaBand = [4 12]
        betaBand = [16 32]
    end
    properties % Organized Spiking Data
        beginTrialMtx
        endTrialMtx
        trialSpikeMtx
        fisSeqSpikeMtx
    end
    properties % Organized LFP Data
        beginTrialThetaPowerMtx
        beginTrialThetaPhaseMtx
        beginTrialBetaPowerMtx
        beginTrialBetaPhaseMtx
        
        endTrialThetaPowerMtx
        endTrialThetaPhaseMtx
        endTrialBetaPowerMtx
        endTrialBetaPhaseMtx
    end
    properties % Data Organization Vectors
        beginTrialTime
        endTrialTime
        fisSeqTrlNums
        osTrlNums
        skpTrlNums
        repTrlNums
        taoTrlNums
        taRepTrlNums
        taSkpTrlNums
    end
    properties % Posterior Decoding Logicals
        fisSeqSpikeTimeLog
        fisSeqSpikeOdorLog
        fisSeqSpikeTrialPeriodLog
    end
    properties % Trial Logicals
        trialPeriodTimeLog
        trialPeriodSize
        trialTimeLog
    end
    properties % Posterior Distributions
        fisSeqPosts
        osPosts
        skpPosts
        repPosts
        taoPosts
        taRepPosts
        taSkpPosts
    end
    properties % Decodings
        fisL1OdecodeOdrTime
        fisL1OdecodeTimeTime
        fisL1OdecodeOdr
        fisL1OdecodeOdr_TrlPrd
        fisL1OdecodeTime
        fisL1OdecodeTime_TrlPrd
        fisL1OdecodePosts
        fisL1OdecodePosts_TrlPrd
        
        osDecode
        osDecode_TrlPrd
        repDecode
        repDecode_TrlPrd
        skpDecode
        skpDecode_TrlPrd
        
        taoDecode
        taoDecode_TrlPrd
        taRepDecode
        taRepDecode_TrlPrd
        taSkpDecode
        taSkpDecode_TrlPrd
        taoDecodeTime
        taoDecodeTime_TrlPrd
    end
    
    %% Methods
    methods
        %% Object Creation
        %#ok<*EXIST>
        function obj = PFC_TrialEvent_MLB_SM(fileDir, binSize, dsRate, beginWindow, endWindow)
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
                end
                if exist('beginWindow')==1
                    obj.beginTrialWindow = beginWindow;
                end
                if exist('endWindow')==1
                    obj.endTrialWindow = endWindow;
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
            obj.CompileMLBmtx(obj.beginTrialAlignment);
            obj.CompileMLBmtx(obj.endTrialAlignment);
            obj.CompileTrialObsvs;
            obj.CompileFISlikes;
            obj.DecodeFIS_L1O;
            obj.DecodeOS;
            obj.DecodeTAO;
        end
        %% Decode Trial After OutSeq from FIS
        function DecodeTAO(obj)
            if isempty(obj.trialSpikeMtx)
                obj.CompileTrialObsvs;
            end
            obj.taoPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.taoTrlNums));
            obj.taRepPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.taRepTrlNums));
            obj.taSkpPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.taSkpTrlNums));
            obj.taoDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.taoTrlNums));
            obj.taoDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.taoTrlNums), 4);
            obj.taRepDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.taRepTrlNums));
            obj.taRepDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.taRepTrlNums), 4);
            obj.taSkpDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.taSkpTrlNums));
            obj.taSkpDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.taSkpTrlNums), 4);
            obj.taoDecodeTime = nan(size(obj.trialSpikeMtx{1},1), length(obj.taoTrlNums));
            obj.taoDecodeTime_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.taoTrlNums), 4);
            for taO = 1:length(obj.taoTrlNums)
                curObsv = obj.trialSpikeMtx{obj.taoTrlNums(taO)};
                obj.taoPosts(:,:,taO) = obj.CalcStaticBayesPost(mean(obj.fisSeqSpikeMtx,3), curObsv);
                [tempDecodeOdr, tempDecodeOdrPost] = obj.DecodeBayesPost(obj.taoPosts(:,:,taO), obj.fisSeqSpikeOdorLog);
                tempDecodeOdr(tempDecodeOdr==obj.trialInfo(obj.taoTrlNums(taO)).Position) = 0;
                tempDecodeOdr(tempDecodeOdr==obj.trialInfo(obj.taoTrlNums(taO)-1).Position) = -1;
                tempDecodeOdr(tempDecodeOdr==obj.trialInfo(obj.taoTrlNums(taO)-1).Odor) = -2;
                tempDecodeOdr(tempDecodeOdr>0) = 1;
                
                obj.taoDecode(:,taO) = tempDecodeOdr;
                
                obj.taoDecodeTime(:,taO) = obj.DecodeBayesPost(obj.taoPosts(:,:,taO), obj.fisSeqSpikeTimeLog) - obj.trialTimeLog;
                for prd = 1:4
                    obj.taoDecode_TrlPrd(1:obj.trialPeriodSize(prd),taO,prd) = tempDecodeOdr(obj.trialPeriodTimeLog==prd);
                    obj.taoDecodeTime_TrlPrd(1:obj.trialPeriodSize(prd),taO,prd) = obj.taoDecodeTime(obj.trialPeriodTimeLog==prd,taO);
                end
                
                if sum(obj.taoTrlNums(taO)==obj.taSkpTrlNums)==1
                    obj.taSkpPosts(:,:,obj.taoTrlNums(taO)==obj.taSkpTrlNums) = obj.taoPosts(:,:,taO);
                    obj.taSkpDecode(:,obj.taoTrlNums(taO)==obj.taSkpTrlNums) = obj.taoDecode(:,taO);
                    obj.taSkpDecode_TrlPrd(:,obj.taoTrlNums(taO)==obj.taSkpTrlNums,:) = obj.taoDecode_TrlPrd(:,taO,:);
                elseif sum(obj.taoTrlNums(taO)==obj.taRepTrlNums)==1
                    obj.taRepPosts(:,:,obj.taoTrlNums(taO)==obj.taRepTrlNums) = obj.taoPosts(:,:,taO);
                    obj.taRepDecode(:,obj.taoTrlNums(taO)==obj.taRepTrlNums) = obj.taoDecode(:,taO);
                    obj.taRepDecode_TrlPrd(:,obj.taoTrlNums(taO)==obj.taRepTrlNums,:) = obj.taoDecode_TrlPrd(:,taO,:);
                end
            end
        end
        %% Decode OutSeq from FIS
        function DecodeOS(obj)
            if isempty(obj.trialSpikeMtx)
                obj.CompileTrialObsvs;
            end
            obj.osPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.osTrlNums));
            obj.osDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.osTrlNums));
            obj.osDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.osTrlNums), 4);
            obj.skpPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.skpTrlNums));
            obj.skpDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.skpTrlNums));
            obj.skpDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.skpTrlNums), 4);
            obj.repPosts = nan(size(obj.trialSpikeMtx{1},1), size(obj.fisSeqSpikeMtx,1), length(obj.repTrlNums));
            obj.repDecode = nan(size(obj.trialSpikeMtx{1},1), length(obj.repTrlNums));
            obj.repDecode_TrlPrd = nan(max(obj.trialPeriodSize), length(obj.repTrlNums), 4);
            for osT = 1:length(obj.osTrlNums)
                curObsv = obj.trialSpikeMtx{obj.osTrlNums(osT)};
                obj.osPosts(:,:,osT) = obj.CalcStaticBayesPost(mean(obj.fisSeqSpikeMtx,3), curObsv);
                [tempDecode, tempDecodePost] = obj.DecodeBayesPost(obj.osPosts(:,:,osT), obj.fisSeqSpikeOdorLog);
                tempDecode(tempDecode==obj.trialInfo(obj.osTrlNums(osT)).Position) = 0;
                tempDecode(tempDecode==obj.trialInfo(obj.osTrlNums(osT)).Odor) = -1;
                tempDecode(tempDecode>0) = 1;
                
                obj.osDecode(:,osT) = tempDecode;
                for prd = 1:4
                    obj.osDecode_TrlPrd(1:obj.trialPeriodSize(prd),osT,prd) = tempDecode(obj.trialPeriodTimeLog==prd);
                end
                
                if sum(obj.osTrlNums(osT)==obj.skpTrlNums)==1
                    obj.skpPosts(:,:,obj.osTrlNums(osT)==obj.skpTrlNums) = obj.osPosts(:,:,osT);
                    obj.skpDecode(:,obj.osTrlNums(osT)==obj.skpTrlNums) = obj.osDecode(:,osT);
                    obj.skpDecode_TrlPrd(:,obj.osTrlNums(osT)==obj.skpTrlNums,:) = obj.osDecode_TrlPrd(:,osT,:);
                elseif sum(obj.osTrlNums(osT)==obj.repTrlNums)==1
                    obj.repPosts(:,:,obj.osTrlNums(osT)==obj.repTrlNums) = obj.osPosts(:,:,osT);
                    obj.repDecode(:,obj.osTrlNums(osT)==obj.repTrlNums) = obj.osDecode(:,osT);
                    obj.repDecode_TrlPrd(:,obj.osTrlNums(osT)==obj.repTrlNums,:) = obj.osDecode_TrlPrd(:,osT,:);
                end
            end
        end
        %% Decode FIS via Leave-1-Out
        function DecodeFIS_L1O(obj)
            if isempty(obj.fisSeqSpikeMtx)
                obj.CompileFISlikes;
            end
            obj.fisSeqPosts = nan(size(obj.fisSeqSpikeMtx,1), size(obj.fisSeqSpikeMtx,1), size(obj.fisSeqSpikeMtx,3));
            for s = 1:size(obj.fisSeqSpikeMtx,3)
                tempISS = obj.fisSeqSpikeMtx;
                tempISS(:,:,s) = [];
                obj.fisSeqPosts(:,:,s) = obj.CalcStaticBayesPost(mean(tempISS,3), obj.fisSeqSpikeMtx(:,:,s));
            end
            [tempDecodeOdrTm, tempDecodeFISpost] = obj.DecodeBayesPost(obj.fisSeqPosts,obj.fisSeqSpikeOdorLog);
            obj.fisL1OdecodeOdrTime = nan(size(tempDecodeOdrTm,1),4);
            for o = 1:4
                obj.fisL1OdecodeOdrTime(:,o) = sum(tempDecodeOdrTm==o,2)./sum(~isnan(tempDecodeOdrTm),2);
            end
            [tempDecodeTime, ~] = obj.DecodeBayesPost(obj.fisSeqPosts, obj.fisSeqSpikeTimeLog);
            for s = 1:size(tempDecodeTime,2)
                tempDecodeTime(:,s) = tempDecodeTime(:,s)-obj.fisSeqSpikeTimeLog;
            end
            obj.fisL1OdecodeTimeTime = tempDecodeTime;
            
            obj.fisL1OdecodeOdr = nan(4);
            obj.fisL1OdecodeTime = nan(4);
            obj.fisL1OdecodePosts = nan(4);
            obj.fisL1OdecodeOdr_TrlPrd = nan(4,4,4);
            obj.fisL1OdecodeTime_TrlPrd = nan(4,4,4);
            obj.fisL1OdecodePosts_TrlPrd = nan(4,4,4);
            for o1 = 1:4
                for o2 = 1:4
                    odrMtx = tempDecodeOdrTm(obj.fisSeqSpikeOdorLog==o1,:);
                    timeMtx = tempDecodeTime(obj.fisSeqSpikeOdorLog==o1,:);
                    postMtx = tempDecodeFISpost(obj.fisSeqSpikeOdorLog==o1,:);
                    obj.fisL1OdecodeOdr(o2,o1) = mean(mean(odrMtx==o2));
                    obj.fisL1OdecodeTime(o2,o1) = mean(timeMtx(odrMtx==o2));
                    obj.fisL1OdecodePosts(o2,o1) = mean(postMtx(odrMtx==o2));
                    for prd = 1:4
                        obj.fisL1OdecodeOdr_TrlPrd(o2,o1,prd) = mean(mean(odrMtx(obj.trialPeriodTimeLog==prd,:)==o2));
                        obj.fisL1OdecodeTime_TrlPrd(o2,o1,prd) = mean(timeMtx(odrMtx(obj.trialPeriodTimeLog==prd,:)==o2));
                        obj.fisL1OdecodePosts_TrlPrd(o2,o1,prd) = mean(postMtx(odrMtx(obj.trialPeriodTimeLog==prd,:)==o2));
                    end
                end
            end
        end
        %% Compile Trial Observations (Combine PokeIn and PokeOut across trials)
        function CompileTrialObsvs(obj)
            obj.trialPeriodTimeLog = [[obj.beginTrialTime<0; false(size(obj.endTrialTime))],...
                [obj.beginTrialTime>0; false(size(obj.endTrialTime))],...
                [false(size(obj.beginTrialTime)); obj.endTrialTime<0],...
                [false(size(obj.beginTrialTime)); obj.endTrialTime>0]];
            obj.trialPeriodTimeLog = obj.trialPeriodTimeLog*(1:4)';
            obj.trialPeriodTimeLog(obj.trialPeriodTimeLog==0) = nan;
            obj.trialPeriodSize = nan(1,4);
            for prd = 1:4
                obj.trialPeriodSize(prd) = sum(obj.trialPeriodTimeLog==prd);
            end
                
            obj.trialTimeLog = [obj.beginTrialTime; obj.endTrialTime+1+obj.dsRate/obj.sampleRate];
            obj.osTrlNums = [obj.trialInfo([obj.trialInfo.TranspositionDistance]~=0).TrialNum];
            obj.osTrlNums([obj.trialInfo(obj.osTrlNums).Performance]==0) = [];
            obj.skpTrlNums = obj.osTrlNums([obj.trialInfo(obj.osTrlNums).TranspositionDistance]<0);
            obj.repTrlNums = obj.osTrlNums([obj.trialInfo(obj.osTrlNums).TranspositionDistance]>0);
            obj.taoTrlNums = obj.osTrlNums(([obj.trialInfo(obj.osTrlNums).Position]+1 == [obj.trialInfo(obj.osTrlNums+1).Position])...
                & [obj.trialInfo(obj.osTrlNums+1).Performance]==1)+1;
            obj.taRepTrlNums = obj.osTrlNums(([obj.trialInfo(obj.osTrlNums).Position]+1 == [obj.trialInfo(obj.osTrlNums+1).Position])...
                & [obj.trialInfo(obj.osTrlNums+1).Performance]==1 ...
                & [obj.trialInfo(obj.osTrlNums).TranspositionDistance]>0)+1;
            obj.taSkpTrlNums = obj.osTrlNums(([obj.trialInfo(obj.osTrlNums).Position]+1 == [obj.trialInfo(obj.osTrlNums+1).Position])...
                & [obj.trialInfo(obj.osTrlNums+1).Performance]==1 ...
                & [obj.trialInfo(obj.osTrlNums).TranspositionDistance]<0)+1;
            obj.trialSpikeMtx = cell(size(obj.trialInfo));
            for t = 1:length(obj.trialInfo)
                obj.trialSpikeMtx{t} = [obj.beginTrialMtx(:,:,t); obj.endTrialMtx(:,:,t)];
            end
        end
        %% Compile Fully InSeq Likelihoods (trial PSTHs)
        function CompileFISlikes(obj)
            if isempty(obj.beginTrialMtx)
                error('Compile Poke In data');
            elseif isempty(obj.endTrialMtx)
                error('Complie Poke Out data');
            end
            if isempty(obj.trialPeriodTimeLog)
                obj.CompileTrialObsvs;
            end
            fisStarts = find(conv([obj.trialInfo.Odor], 1:4, 'valid')==20 &...
                conv([obj.trialInfo.Position], 1:4, 'valid')==20 &...
                conv([obj.trialInfo.Performance], ones(1,4), 'valid')==4);
            obj.fisSeqTrlNums = nan(4,length(fisStarts));
            for iS = 1:length(fisStarts)
                obj.fisSeqTrlNums(1,iS) = fisStarts(iS);
                obj.fisSeqTrlNums(2,iS) = fisStarts(iS) + 1;
                obj.fisSeqTrlNums(3,iS) = fisStarts(iS) + 2;
                obj.fisSeqTrlNums(4,iS) = fisStarts(iS) + 3;
            end
            obj.fisSeqSpikeMtx = nan((size(obj.beginTrialMtx,1) + size(obj.endTrialMtx,1))*4, size(obj.beginTrialMtx,2), length(fisStarts));
            for seq = 1:length(fisStarts)
                obj.fisSeqSpikeMtx(:,:,seq) = [obj.beginTrialMtx(:,:,obj.fisSeqTrlNums(1,seq)); obj.endTrialMtx(:,:,obj.fisSeqTrlNums(1,seq));...
                    obj.beginTrialMtx(:,:,obj.fisSeqTrlNums(2,seq)); obj.endTrialMtx(:,:,obj.fisSeqTrlNums(2,seq));...
                    obj.beginTrialMtx(:,:,obj.fisSeqTrlNums(3,seq)); obj.endTrialMtx(:,:,obj.fisSeqTrlNums(3,seq));...
                    obj.beginTrialMtx(:,:,obj.fisSeqTrlNums(4,seq)); obj.endTrialMtx(:,:,obj.fisSeqTrlNums(4,seq))];
            end
            obj.fisSeqSpikeOdorLog = [ones(1,size(obj.beginTrialMtx,1)+size(obj.endTrialMtx,1))*1,...
                ones(1,size(obj.beginTrialMtx,1)+size(obj.endTrialMtx,1))*2,...
                ones(1,size(obj.beginTrialMtx,1)+size(obj.endTrialMtx,1))*3,...
                ones(1,size(obj.beginTrialMtx,1)+size(obj.endTrialMtx,1))*4]';
            obj.fisSeqSpikeTimeLog = repmat([obj.beginTrialTime; obj.endTrialTime+1+obj.dsRate/obj.sampleRate], [4,1]);
            obj.fisSeqSpikeTrialPeriodLog = repmat(obj.trialPeriodTimeLog, [4,1]);
        end
        %% Extract Matrix
        function CompileMLBmtx(obj, alignment)
            switch alignment
                case obj.beginTrialAlignment
                    [obj.beginTrialMtx, obj.beginTrialTime] = obj.PP_TrialMatrix_Spiking(obj.beginTrialWindow, alignment);
                    [obj.beginTrialThetaPhaseMtx, obj.beginTrialThetaPowerMtx] = obj.PP_TrialMatrix_LFP(obj.thetaBand, obj.beginTrialWindow, alignment);
                    [obj.beginTrialBetaPhaseMtx, obj.beginTrialBetaPowerMtx] = obj.PP_TrialMatrix_LFP(obj.betaBand, obj.beginTrialWindow, alignment);
                case obj.endTrialAlignment
                    [obj.endTrialMtx, obj.endTrialTime] = obj.PP_TrialMatrix_Spiking(obj.endTrialWindow, alignment);
                    [obj.endTrialThetaPhaseMtx, obj.endTrialThetaPowerMtx] = obj.PP_TrialMatrix_LFP(obj.thetaBand, obj.endTrialWindow, alignment);
                    [obj.endTrialBetaPhaseMtx, obj.endTrialBetaPowerMtx] = obj.PP_TrialMatrix_LFP(obj.betaBand, obj.endTrialWindow, alignment);
            end
        end
    end
end