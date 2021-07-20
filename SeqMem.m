classdef SeqMem < handle
    properties (SetAccess = protected)
        pathDir
        behavMatFile
        nsmblMatFile
        smFiles
        seqLength
        sampleRate
        tsVect
        behavMatrix
        behavMatrixColIDs
        ensembleMatrix
        ensembleMatrixColIDs
        lfpMatrix
        lfpMatrixColIDs
        trialInfo
        unitInfo
        rejectedLFPmask
    end
    properties (Constant)
        PositionColors = [44/255, 168/255, 224/255;...
            154/255, 133/255, 122/255;...
            9/255, 161/255, 74/255;...
            128/255, 66/255, 151/255];
    end
    %% Object Creation Method
    methods
        function obj = SeqMem(path)
            fprintf('Compiling StatMatrix Data\n');
            if nargin == 0
                path = uigetdir;
            end
            files = dir(path);
            fileNames = {files.name};
            obj.pathDir = path;
            obj.behavMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
            behavMatrix = load([obj.pathDir '\' obj.behavMatFile], 'behavMatrix');
            obj.tsVect = behavMatrix.behavMatrix(:,1);
            obj.behavMatrix = behavMatrix.behavMatrix(:,2:end);
            behavMatrixColIDs = load([obj.pathDir '\' obj.behavMatFile], 'behavMatrixColIDs');
            obj.behavMatrixColIDs = behavMatrixColIDs.behavMatrixColIDs(2:end);
            obj.nsmblMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))};
            ensembleMatrix = load([obj.pathDir '\' obj.nsmblMatFile]);
            obj.ensembleMatrix = ensembleMatrix.ensembleMatrix(:,2:end);
            obj.ensembleMatrixColIDs = ensembleMatrix.ensembleMatrixColIDs(2:end);
            obj.unitInfo = ensembleMatrix.ensembleUnitSummaries;
            obj.smFiles = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';
            obj.OrganizeTrialInfo;
%             obj.CompileLFPmatrix;
        end
    end
    %% Data Organization & Extraction
    methods
        %% OrganizeTrialInfo
        function OrganizeTrialInfo(obj)
            obj.sampleRate = 1/mode(diff(obj.tsVect));
            
            %% Extract Trial Indexes & Poke Events
            % Separate Odor & Position columns
            odorTrlMtx = obj.behavMatrix(:,cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'Odor')));
            odorVals = cell2mat(cellfun(@(b)str2double(b(5:end)), obj.behavMatrixColIDs(cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'Odor'))), 'uniformoutput', 0));
            positionTrlMtx = obj.behavMatrix(:,cellfun(@(a)~isempty(a), regexp(obj.behavMatrixColIDs, 'Position[1-9]$')));
            % Sum them on 2d to extract trial indices
            trialVect = sum(odorTrlMtx,2);
            trialIndices = find(trialVect);
            numTrials = sum(trialVect);
            % Pull out Poke events and identify pokeIn/Out indices
            pokeVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'PokeEvents')));
            pokeInNdxs = find(pokeVect==1);
            pokeOutNdxs = find(pokeVect==-1);
            % Pull out Front Reward Indices
            frontRwrdVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'FrontReward')));
            frontRwrdNdxs = find(frontRwrdVect);
            % Pull out Rear Reward Indices
            rearRwrdVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'BackReward')));
            rearRwrdNdxs = find(rearRwrdVect);
            % Pull out Error Indices
            errorSigVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'ErrorSignal')));
            errorSigNdxs = find(errorSigVect);
            % Pull out Reward Signal Indices
            rwdSigVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'RewardSignal')));
            rewardSigNdxs = find(rwdSigVect);
            % Identify trial performance
            trialPerfVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'Performance')));
            % Identify InSeq logical
            inSeqLog = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'InSeqLog')));
            %% Create Trial Info Structure            
            seqNum = cell(1,numTrials);
            trialOdor = cell(1,numTrials);
            trialPosition = cell(1,numTrials);
            trialPerf = cell(1,numTrials);
            trialTransDist = cell(1,numTrials);
            trialItmItmDist = cell(1,numTrials);
            trialPokeInNdx = repmat({nan}, [1, numTrials]);
            trialOdorNdx = repmat({nan}, [1, numTrials]);
            trialPokeOutNdx = repmat({nan}, [1, numTrials]);
            trialRewardNdx = repmat({nan}, [1, numTrials]);
            trialErrorNdx = repmat({nan}, [1, numTrials]);
            trialRwdSigNdx = repmat({nan}, [1, numTrials]);
            trialRearRwdNdx = repmat({nan}, [1, numTrials]);
            trialNum = cell(1,numTrials);
            pokeDuration = cell(1,numTrials);
            withdrawLat = cell(1,numTrials);
            seq = 0;
            %% Go through each trial and pull out trial information and create a logical vector for that trial's time periods specified by the input trialLims
            for trl = 1:numTrials
                trialNum{trl} = trl;
                % Identify Trial/Position/Descriptors
                curTrlOdor = odorVals(odorTrlMtx(trialIndices(trl),:)==1);
                curTrlPos = find(positionTrlMtx(trialIndices(trl),:)==1);
                curTrlPerf = trialPerfVect(trialIndices(trl))==1;
                curTrlInSeqLog = inSeqLog(trialIndices(trl))==1;
                
                trialOdor{trl} = curTrlOdor;
                trialPosition{trl} = curTrlPos;
                trialPerf{trl} = curTrlPerf;
                % Increment the sequence counter as necessary
                if curTrlPos==1                                                         % Increment if trial is in the first position
                    seq = seq+1;
                elseif trl==1 && curTrlPos ~= 1                                         % Also increment if it's the first trial in the session but the position is not 1 (happens when first trial is curated out)
                    seq = seq+1;
                elseif curTrlPos <= trialPosition{trl-1}                                % Only gets here if the first trial in a sequence was curated out
                    seq = seq+1;
                end
                % Identify temporal context feature
                if curTrlInSeqLog
                    trialItmItmDist{trl} = 1;
                    trialTransDist{trl} = 0;
                else
                    if curTrlOdor < 10
                        trialTransDist{trl} = curTrlPos - curTrlOdor;
                        trialItmItmDist{trl} = curTrlOdor - curTrlPos + 1;
                    else
                        trialTransDist{trl} = curTrlPos+10 - curTrlOdor;
                        trialItmItmDist{trl} = curTrlOdor - (curTrlPos+10) + 1;
                    end
                end
                seqNum{trl} = seq;
                
                % Create trial logical vector
                curPokeIn = pokeInNdxs(find(pokeInNdxs<=trialIndices(trl)==1,1, 'last'));
                curPokeOut = pokeOutNdxs(find(pokeOutNdxs>trialIndices(trl)==1,1, 'first'));
                
                trialPokeInNdx{trl} = curPokeIn;
                trialOdorNdx{trl} = trialIndices(trl);
                trialPokeOutNdx{trl} = curPokeOut;
                curFrontRwrdNdx = frontRwrdNdxs(find(frontRwrdNdxs>trialIndices(trl)==1,1, 'first'));
                if  isempty(curFrontRwrdNdx) || trl==numTrials || curFrontRwrdNdx<trialIndices(trl+1)
                    trialRewardNdx{trl} = curFrontRwrdNdx;
                    if isempty(curFrontRwrdNdx)
                        trialRewardNdx{trl} = nan;
                    end
                else
                    trialRewardNdx{trl} = nan;
                end
                curErrSigNdx = errorSigNdxs(find(errorSigNdxs>trialIndices(trl)==1,1,'first'));
                if isempty(curErrSigNdx) || trl==numTrials || curErrSigNdx<trialIndices(trl+1)
                    trialErrorNdx{trl} = curErrSigNdx;
                    if isempty(curErrSigNdx)
                        trialErrorNdx{trl} = nan;
                    end
                else
                    trialErrorNdx{trl} = nan;
                end
                curRwdSigNdx = rewardSigNdxs(find(rewardSigNdxs>trialIndices(trl)==1,1,'first'));
                if isempty(curRwdSigNdx) || trl==numTrials || curRwdSigNdx<trialIndices(trl+1)
                    trialRwdSigNdx{trl} = curRwdSigNdx;
                    if isempty(curRwdSigNdx)
                        trialRwdSigNdx{trl} = nan;
                    end
                else
                    trialRwdSigNdx{trl} = nan;
                end
                curRearRwdNdx = rearRwrdNdxs(find(rearRwrdNdxs>trialIndices(trl)==1,1,'first'));
                if isempty(curRearRwdNdx) || trl==numTrials || curRearRwdNdx<trialIndices(trl+1)
                    trialRearRwdNdx{trl} = curRearRwdNdx;
                    if isempty(curRearRwdNdx)
                        trialRearRwdNdx{trl} = nan;
                    end
                else
                    trialRearRwdNdx{trl} = nan;
                end
                pokeDuration{trl} = (curPokeOut-curPokeIn)/obj.sampleRate;
                withdrawLat{trl} = (curPokeOut-trialRwdSigNdx{trl})/obj.sampleRate;
            end
            %% Create obj.behavMatrixTrialStruct
            obj.trialInfo = struct( 'TrialNum', trialNum, 'SequenceNum', seqNum,...
                'Odor', trialOdor, 'Position', trialPosition, 'Performance', trialPerf,...
                'PokeDuration', pokeDuration, 'WithdrawLatency', withdrawLat,...
                'PokeInIndex', trialPokeInNdx, 'OdorIndex', trialOdorNdx, 'PokeOutIndex', trialPokeOutNdx,...
                'RewardIndex', trialRewardNdx, 'RewardSignalIndex', trialRwdSigNdx,...
                'RearRewardIndex', trialRearRwdNdx, 'ErrorIndex', trialErrorNdx,...
                'TranspositionDistance', trialTransDist, 'ItemItemDistance', trialItmItmDist);
            obj.seqLength = size(positionTrlMtx,2);
        end
        %% Create LFP matrix
        function CompileLFPmatrix(obj)
            fprintf('Compiling LFP\n');
            obj.lfpMatrix = nan(size(obj.behavMatrix,1),length(obj.smFiles));
            obj.lfpMatrixColIDs = cell(1,length(obj.smFiles));
            for file = 1:length(obj.smFiles)
                sm = load([obj.pathDir '\' obj.smFiles{file}], 'statMatrix');
                colIDs = load([obj.pathDir '\' obj.smFiles{file}], 'statMatrixColIDs');
                idSplit = strsplit(colIDs.statMatrixColIDs{2}, '_');
                obj.lfpMatrixColIDs{file} = idSplit{1};
                obj.lfpMatrix(:,file) = sm.statMatrix(:,2);                
            end
        end
        %% Extract Trial Matrix
        function [trialMatrix] = ExtractTrialMatrix(obj, dataMatrix, window, alignment)
            if nargin == 1
                dataMatrix = obj.ensembleMatrix;
                window = [-500 500];
                alignment = 'PokeIn';
            end
            trialOrgData = cell(1,1,length(obj.trialInfo));
            for t = 1:length(obj.trialInfo)
                tempLogVect = false(size(obj.behavMatrix,1),1);
                switch alignment
                    case 'Odor'
                        curIndex = obj.trialInfo(t).OdorIndex;
                    case 'PokeIn'
                        curIndex = obj.trialInfo(t).PokeInIndex;
                    case 'PokeOut'
                        curIndex = obj.trialInfo(t).PokeOutIndex;
                    case 'FrontReward'
                        curIndex = obj.trialInfo(t).RewardIndex;
                    case 'RewardSignal'
                        curIndex = obj.trialInfo(t).RewardSignalIndex;
                    case 'RearReward'
                        curIndex = obj.trialInfo(t).RearRewardIndex;
                    case 'ErrorSignal'
                        curIndex = obj.trialInfo(t).ErrorIndex;
                    otherwise
                        error('invalid alignment used');
                end                
                if ~isnan(curIndex)
                    curWindow = curIndex + window;
                    tempLogVect(curWindow(1):curWindow(2)) = true;
                end
                trialOrgData{t} = dataMatrix(tempLogVect,:);                
            end
            trialMatrix = cell2mat(trialOrgData);
        end
    end
    %% Data Pre-Processing
    methods
        %% Simple Filtering & Instantaneous Phase Extraction
        function [filtSig, hilbPhase] = SimpleFilter(obj, signal, freqWin)
            Wn_FRange = [freqWin(1)/(obj.sampleRate/2) freqWin(2)/(obj.sampleRate/2)];
            [bFRange, aFRange] = butter(3, Wn_FRange);
            filtSig = filtfilt(bFRange, aFRange, signal);
            hilbPhase = atan2(imag(hilbert(filtSig)), filtSig);
        end
        %% LFP Artifact Rejection
    end
end