classdef SeqMem < handle
    properties % Files used to create the object
        pathDir
        behavMatFile
        nsmblMatFile
        smFiles
    end
    properties % Properties about the session files
        seqLength
        numSeqs
        odrSeqs
        odrVect
        lagVect
        sampleRate
    end
    properties % Raw data compiled from session files
        tsVect
        behavMatrix
        behavMatrixColIDs
        plxData
        ensembleMatrix
        ensembleMatrixColIDs
        lfpMatrix
        lfpMatrixColIDs
    end
    properties % Organized data from session files
        trialInfo
        unitInfo
        lfpRefTet
        numUniPerTet        
    end
    properties % Calculate behavioral variables
        isTrialNums
        % Hit/FP/FN/Miss tables
        responseMatrix
        responseMatrixSFP
        responseMatrixByPos
        responseMatrixByOdr
        % Raw Accuracy: Trial Type & By Lag
        transMatAcc
        transMatPerf
        lagAccVect
        lagAccVectSFP
        % Response Latency: Trial Type & By Lag
        transMatLatRaw
        transMatLatMean
        lagLatVectRaw
        lagLatVectMean
        lagLatVectSFPraw        
        lagLatVectSFPmean 
        % Trial After OutSeq: Accuracy & Latency
        taoAcc
        taoLatRaw
        taoLatMean
        % d' Accuracy
        dPrime
        dPrimeSFP
        dPrimeByPos
        dPrimeByOdr
        % smi Accuracy
        smi
        smiSFP
        smiByPos
        smiByOdr
        % ri Response Bias
        ri
        riSFP
        riByPos
        riByOdr
        % OrdX Trial Analysis
        ordXpokeLats
    end
    properties % Masks to select data
        rejectedLFPmask
        popVectIncludeLog
    end
    properties (Constant)
        PositionColors = [44/255, 168/255, 224/255;...
            154/255, 133/255, 122/255;...
            9/255, 161/255, 74/255;...
            128/255, 66/255, 151/255;...
            241/255, 103/255, 36/255];
        PPFcolors = [255/255, 230/255, 16/255;...
            255/255, 128/255, 3/255;...
            195/255, 0/255, 23/255];
        Rosetta = [{'A'},{'B'},{'C'},{'D'},{'E'},{'F'},{'G'},{'H'},{'I'},{'J'},{'K'},{'L'},{'M'},{'N'},{'O'},{'P'},{'Q'},{'R'},{'S'},{'T'},{'U'},{'V'},{'W'},{'X'},{'Y'},{'Z'}];
        RosettaLC = [{'a'},{'b'},{'c'},{'d'},{'e'},{'f'},{'g'},{'h'},{'i'},{'j'},{'k'},{'l'},{'m'},{'n'},{'o'},{'p'},{'q'},{'r'},{'s'},{'t'},{'u'},{'v'},{'w'},{'x'},{'y'},{'z'}];
    end
    %% Object Creation Method
    methods
        function obj = SeqMem(path)
            fprintf('Compiling StatMatrix Data....');
            if nargin == 0
                path = uigetdir;
            end
            files = dir(path);
            fileNames = {files.name};
            obj.pathDir = path;
            obj.behavMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
            behavMatrix = load([obj.pathDir '\' obj.behavMatFile], 'behavMatrix');
            plxData = load([obj.pathDir '\' obj.behavMatFile], 'plxData');
            obj.tsVect = behavMatrix.behavMatrix(:,1);
            if ~isempty(fieldnames(plxData))
                obj.plxData = plxData.plxData;
            end
            obj.behavMatrix = behavMatrix.behavMatrix(:,2:end);
            behavMatrixColIDs = load([obj.pathDir '\' obj.behavMatFile], 'behavMatrixColIDs');
            obj.behavMatrixColIDs = behavMatrixColIDs.behavMatrixColIDs(2:end);
            nsmblMatFileLog = cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'));
            if sum(nsmblMatFileLog)~=0
                obj.nsmblMatFile = fileNames{nsmblMatFileLog};
                ensembleMatrix = load([obj.pathDir '\' obj.nsmblMatFile]);
                if isempty(obj.popVectIncludeLog)
                    obj.ensembleMatrix = ensembleMatrix.ensembleMatrix(:,2:end);
                    obj.ensembleMatrixColIDs = ensembleMatrix.ensembleMatrixColIDs(2:end);
                    obj.unitInfo = ensembleMatrix.ensembleUnitSummaries;
                else
                    obj.ensembleMatrix = ensembleMatrix.ensembleMatrix(:,[false; obj.popVectIncludeLog(:)]);
                    obj.ensembleMatrixColIDs = ensembleMatrix.ensembleMatrixColIDs([false; obj.popVectIncludeLog(:)]);
                    obj.unitInfo = ensembleMatrix.ensembleUnitSummaries(1,obj.popVectIncludeLog);
                end
                sortedPVvect = sortrows([mean(obj.ensembleMatrix); 1:size(obj.ensembleMatrix,2)]', 'descend');
                obj.ensembleMatrix = obj.ensembleMatrix(:,sortedPVvect(:,2));
                obj.ensembleMatrixColIDs = obj.ensembleMatrixColIDs(1,sortedPVvect(:,2));
                obj.unitInfo = obj.unitInfo(1,sortedPVvect(:,2));
            end
            obj.smFiles = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';
            obj.OrganizeTrialInfo;
            obj.SummarizeSessionBehavior;
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
            % Pull out port light events
            portLightVect = obj.behavMatrix(:, cellfun(@(a)~isempty(a), strfind(obj.behavMatrixColIDs, 'PortLight')));
            portLightNdxs = find(portLightVect);
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
            trialPortLightNdx = repmat({nan}, [1, numTrials]);
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
%                 curTrlInSeqLog = inSeqLog(trialIndices(trl))==1;
                
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
                if curTrlPos == 1 || trl==1
                    trialItmItmDist{trl} = 1;
                else
                    trialItmItmDist{trl} = curTrlOdor-trialOdor{trl-1};
                end
                if curTrlOdor < 10
                    trialTransDist{trl} = curTrlPos - curTrlOdor;
                else
                    trialTransDist{trl} = curTrlPos+10 - curTrlOdor;
                end
                seqNum{trl} = seq;
                
                % Create trial logical vector
                curPokeIn = pokeInNdxs(find((pokeInNdxs<=trialIndices(trl))==1,1, 'last'));
                curPokeOut = pokeOutNdxs(find((pokeOutNdxs>trialIndices(trl))==1,1, 'first'));

                curPortLightNdx = portLightNdxs(find((portLightNdxs<trialIndices(trl))==1,1,'last'));
                if isempty(curPortLightNdx)
                    trialPortLightNdx{trl} = nan;
                else
                    trialPortLightNdx{trl} = curPortLightNdx;
                end
                
                trialPokeInNdx{trl} = curPokeIn;
                trialOdorNdx{trl} = trialIndices(trl);
                trialPokeOutNdx{trl} = curPokeOut;
                curFrontRwrdNdx = frontRwrdNdxs(find((frontRwrdNdxs>trialIndices(trl))==1,1, 'first'));
                if  isempty(curFrontRwrdNdx) || trl==numTrials || curFrontRwrdNdx<trialIndices(trl+1)
                    trialRewardNdx{trl} = curFrontRwrdNdx;
                    if isempty(curFrontRwrdNdx)
                        trialRewardNdx{trl} = nan;
                    end
                else
                    trialRewardNdx{trl} = nan;
                end
                curErrSigNdx = errorSigNdxs(find((errorSigNdxs>trialIndices(trl))==1,1,'first'));
                if isempty(curErrSigNdx) || trl==numTrials || curErrSigNdx<trialIndices(trl+1)
                    trialErrorNdx{trl} = curErrSigNdx;
                    if isempty(curErrSigNdx)
                        if ~curTrlPerf
                            trialErrorNdx{trl} = curPokeOut;
                        else
                            trialErrorNdx{trl} = nan;
                        end
                    end
                else
                    trialErrorNdx{trl} = nan;
                end
                curRwdSigNdx = rewardSigNdxs(find((rewardSigNdxs>trialIndices(trl))==1,1,'first'));
                if isempty(curRwdSigNdx) || trl==numTrials || curRwdSigNdx<trialIndices(trl+1)
                    trialRwdSigNdx{trl} = curRwdSigNdx;
                    if isempty(curRwdSigNdx)
                        trialRwdSigNdx{trl} = nan;
                    end
                else
                    trialRwdSigNdx{trl} = nan;
                end
                curRearRwdNdx = rearRwrdNdxs(find((rearRwrdNdxs>trialIndices(trl))==1,1,'first'));
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
            if ~isempty(obj.plxData)
                obj.trialInfo = struct( 'TrialNum', trialNum, 'SequenceNum', seqNum,...
                    'Odor', trialOdor, 'Position', trialPosition, 'Performance', trialPerf,...
                    'TargetDuration', {obj.plxData.Raw.TargetDuration}, 'PokeDuration', pokeDuration, 'WithdrawLatency', withdrawLat,...
                    'PortLightIndex', trialPortLightNdx, 'PokeInIndex', trialPokeInNdx, 'OdorIndex', trialOdorNdx, 'PokeOutIndex', trialPokeOutNdx,...
                    'RewardIndex', trialRewardNdx, 'RewardSignalIndex', trialRwdSigNdx,...
                    'RearRewardIndex', trialRearRwdNdx, 'ErrorIndex', trialErrorNdx,...
                    'TranspositionDistance', trialTransDist, 'ItemItemDistance', trialItmItmDist);
            else
                obj.trialInfo = struct( 'TrialNum', trialNum, 'SequenceNum', seqNum,...
                    'Odor', trialOdor, 'Position', trialPosition, 'Performance', trialPerf,...
                    'PokeDuration', pokeDuration, 'WithdrawLatency', withdrawLat,...
                    'PortLightIndex', trialPortLightNdx, 'PokeInIndex', trialPokeInNdx, 'OdorIndex', trialOdorNdx, 'PokeOutIndex', trialPokeOutNdx,...
                    'RewardIndex', trialRewardNdx, 'RewardSignalIndex', trialRwdSigNdx,...
                    'RearRewardIndex', trialRearRwdNdx, 'ErrorIndex', trialErrorNdx,...
                    'TranspositionDistance', trialTransDist, 'ItemItemDistance', trialItmItmDist);
            end
            obj.seqLength = size(positionTrlMtx,2);
            obj.lagVect = (1:obj.seqLength*2-1)-obj.seqLength;
            obj.numSeqs = length(unique(cell2mat(trialOdor)))/length(unique(cell2mat(trialPosition)));
            obj.odrSeqs = reshape(unique(cell2mat(trialOdor)), [obj.seqLength, obj.numSeqs])';
            obj.odrVect = unique(cell2mat(trialOdor))';
        end
        %% Create LFP matrix
        function CompileLFPmatrix(obj)
            fprintf('Compiling LFP\n');
            obj.lfpMatrix = nan(size(obj.behavMatrix,1),length(obj.smFiles));
            obj.lfpMatrixColIDs = cell(1,length(obj.smFiles));
            obj.numUniPerTet = nan(1,length(obj.smFiles));
            for file = 1:length(obj.smFiles)
                sm = load([obj.pathDir '\' obj.smFiles{file}], 'statMatrix');
                colIDs = load([obj.pathDir '\' obj.smFiles{file}], 'statMatrixColIDs');
                obj.numUniPerTet(file) = sum(cellfun(@(a)~isempty(a), regexp(colIDs.statMatrixColIDs, '-U[0-9]*')));
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
                    curWindow = curIndex + window + obj.alignmentOffset;
                    tempLogVect(curWindow(1):curWindow(2)) = true;
                    trialOrgData{t} = dataMatrix(tempLogVect,:);
                else
                    trialOrgData{t} = nan(length(window(1):window(2)), size(dataMatrix,2));
                end                
            end
            trialMatrix = cell2mat(trialOrgData);
        end
        %% Organize Matrix
        function [dataMtx] = OrganizeDataCellMatrix(~, data, repDim, rowIDs, colIDs)
            rowFacts = unique(rowIDs);
            colFacts = unique(colIDs);
            dataMtx = cell(length(rowFacts), length(colFacts));
            for r = 1:length(rowFacts)
                rowLog = rowIDs==rowFacts(r);
                for c = 1:length(colFacts)
                    colLog = colIDs==colFacts(c);
                    if repDim == 1
                        dataMtx{r,c} = data(rowLog & colLog,:,:,:);
                    elseif repDim == 2
                        dataMtx{r,c} = data(:,rowLog & colLog,:,:);
                    elseif repDim == 3
                        dataMtx{r,c} = data(:,:,rowLog & colLog,:);
                    elseif repDim == 4
                        dataMtx{r,c} = data(:,:,:,rowLog & colLog);
                    end
                end
            end
        end
        %% Build_MtxMaskArray
        function [mask] = BuildMtxMaskArray(~,rowAnchors,rowRanges,rowSz,colAnchors,colRanges,colSz)
            % Inputs:
            %   rowAnchors = vector of length N with anchors for the masks
            %   rowRanges = Ox2 matrix where column 2 = early row anchor offset and column 2 = late anchor offset
            %   rowSz = number of rows in the matrix
            %   colAnchors = vector of length N with anchors for the masks
            %   colRanges = Ox2 matrix where column 2 = early column anchor offset and column 2 = late anchor offset
            %   colSz = number of rows in the matrix
            % Outpus:
            %   mask = M (rowSz) x N (colSz) x O (# anchors) logical array of O MxN masks
            if length(rowAnchors)~=length(colAnchors)
                error('Row and Col anchors of unequal sizes');
            end
            mask = false(rowSz,colSz,length(rowAnchors));
            for ndx = 1:length(colAnchors)
                if rowAnchors(ndx)+rowRanges(ndx,1) < 0
                    rowLow = 0;
                else
                    rowLow = rowAnchors(ndx)+rowRanges(ndx,1);
                end
                if rowAnchors(ndx)+rowRanges(ndx,2) > rowSz
                    rowHigh = rowSz;
                else
                    rowHigh = rowAnchors(ndx)+rowRanges(ndx,2);
                end
                if colAnchors(ndx)+colRanges(ndx,1) < 0
                    colLow = 0;
                else
                    colLow = colAnchors(ndx)+colRanges(ndx,1);
                end
                if colAnchors(ndx)+colRanges(ndx,2) > colSz
                    colHigh = colSz;
                else
                    colHigh = colAnchors(ndx)+colRanges(ndx,2);
                end
                mask(rowLow:rowHigh,colLow:colHigh,ndx) = true;
            end
        end
        %% Build_VectMaskMtx
        function [mask] = BuildVectMaskMtx(~,anchors,ranges,maskSz)
            % Inputs:
            %   anchors = vector of length N of anchor indices for the masks
            %   ranges = Nx2 matrix where column 1 = early anchor offset and column 2 = late anchor offset
            %   maskSz = the size of the data the mask will be applied to
            % Outputs:
            %   mask = MxN logical with M (maskSz) rows and N (# of anchors) columns
            mask = false(maskSz,length(anchors));
            for m = 1:length(anchors)
                if anchors(m)+ranges(m,1) <= 0
                    low = 1;
                else
                    low = anchors(m)+ranges(m,1);
                end
                if anchors(m)+ranges(m,2) > maskSz
                    high = maskSz;
                else
                    high = anchors(m)+ranges(m,2);
                end
                mask(low:high,m) = true;
            end
        end
        %% ArrayMaskExtract
        function [maskedData] = MaskArrayExtract(obj,data2mask,rowMask,colMask,varargin)
            % Apply row and column masks across pages (dim 3) in an array.
            % rowMask = mask applied over the rows, i.e. dim 1
            % colMask = mask applied over the columns, i.e. dim 2
            % Check the masks
            mask4rows = obj.CheckMask(rowMask,size(data2mask,1),size(data2mask,3));
            mask4cols = obj.CheckMask(colMask,size(data2mask,2),size(data2mask,3));
            % Extract data according to the masks
            maskedData = nan(size(data2mask));
            for rep = 1:size(data2mask,3)
                maskedData(mask4rows(:,rep),mask4cols(:,rep),rep) = data2mask(mask4rows(:,rep), mask4cols(:,rep),rep);
            end
            if strcmp(varargin,'extract') % I know this is super messy, so sue me
                extractedData = cell(1,1,size(data2mask,3));
                for rep = 1:size(data2mask,3)
                    temp_maskedData = maskedData(:,:,rep);
                    extractedData{rep} = temp_maskedData(~isnan(temp_maskedData));
                end
                maskedData = cell2mat(extractedData);
            end
        end
        %% RowColMaskExtract
        function [maskedData] = MaskVectExtract(obj,data2mask,mask,dataDim,repDim)
            mask = obj.CheckMask(mask,size(data2mask,dataDim),size(data2mask,repDim));
            if dataDim == 1
                maskedData = nan(size(data2mask,2),size(data2mask,repDim));
            else
                maskedData = nan(size(data2mask,1),size(data2mask,repDim));
            end
            for rep = 1:size(data2mask,repDim)
                if dataDim == 1
                    maskedData(:,rep,:) = data2mask(mask(:,rep),:);
                elseif dataDim == 2
                    maskedData(:,rep,:) = data2mask(:,mask(:,rep));
                end
            end
        end
        %% MaskChecker
        function [newMask] = CheckMask(obj,inputMask,dataLength,numReps)
            if isempty(inputMask)
                newMask = true(dataLength,numReps);
            end
            if isa(inputMask, 'double') || islogical(inputMask)
                newMask = obj.CheckMask_Double(inputMask,dataLength,numReps);
            elseif iscell(inputMask)
                temp_newMask = obj.CheckMask_Double(inputMask{1},dataLength,numReps);
                newMask = false(size(temp_newMask));
                anchor = regexp(inputMask{2}, '[a-zA-Z]*', 'match');
                range = regexp(inputMask{2}, '[+|-]\d', 'match');
                if isempty(range)
                    offset = [0 0];
                elseif length(range)==1
                    offset = [0 str2double(range)];
                elseif length(range)==2
                    offset = [str2double(range{1}), str2double(range{2})];
                end
                for c = 1:size(temp_newMask,2)
                    ndx = find(temp_newMask(:,c),1,anchor{1});
                    newMask(ndx+offset(1):ndx+offset(2),c) = true;
                end
            end
        end
        %% CheckDoubleMask
        function [doubleMask] = CheckMask_Double(~,inputMask,dataLength,numReps)
            if size(inputMask,2)==1 && size(inputMask,2)~=numReps
                doubleMask = repmat(inputMask, [1,numReps]);
            elseif size(inputMask,1) ~= dataLength
                error('Mask Dimension Mismatch');
            else
                doubleMask = inputMask;
            end
        end
    end
    %% Data Pre-Processing
    methods
        %% Simple Filtering & Instantaneous Phase Extraction
        function [filtSig, hilbPhase, hilbPower] = SimpleFilter(obj, signal, freqWin)
            Wn_FRange = [freqWin(1)/(obj.sampleRate/2) freqWin(2)/(obj.sampleRate/2)];
            [bFRange, aFRange] = butter(3, Wn_FRange);
            filtSig = filtfilt(bFRange, aFRange, signal);
            hilbPhase = atan2(imag(hilbert(filtSig)), filtSig);
            hilbPower = abs(hilbert(filtSig));
%             hilbPower = zscore(abs(hilbert(filtSig)));
        end
        %% Diagonal Lag Mean Collapse 
        function [dlmVect, tempLagVect, lagCounts] = DiagonalLagMeanCollapse(~,mtx)
            dlmVect = nan(size(mtx,3),size(mtx,1)*2-1);
            lagCounts = dlmVect;
            tempLagVect = (1-size(mtx,1)):(size(mtx,1)-1);
            mtxLogMask = true(size(mtx,1),size(mtx,2));
            for pg = 1:size(mtx,3)
                tempPage = mtx(:,:,pg);
                for lag = 1:length(tempLagVect)
                    dlmVect(pg,lag) = mean(tempPage(triu(mtxLogMask,tempLagVect(lag)) & tril(mtxLogMask,tempLagVect(lag))), 'omitnan');
                    lagCounts(pg,lag) = sum(sum(triu(mtxLogMask,tempLagVect(lag)) & tril(mtxLogMask,tempLagVect(lag))));
                end
            end            
            lagCounts = lagCounts(1,:);
        end
    end
    %% Behavioral Analyses
    methods
        %% Summarize Session Behavior
        function SummarizeSessionBehavior(obj)
            % Extract Odor & Position Lists
            odrs = unique([obj.trialInfo.Odor]);
            poss = unique([obj.trialInfo.Position]); 
            % Fill in properties with blank data
            obj.isTrialNums = zeros(size(obj.odrSeqs,1), size(obj.odrSeqs,2),2);
            obj.responseMatrix = zeros(2,2,obj.numSeqs);
            obj.responseMatrixSFP = zeros(2,2,obj.numSeqs);
            obj.responseMatrixByPos = repmat({zeros(2,2)}, [obj.numSeqs, obj.seqLength]);
            obj.responseMatrixByOdr = repmat({zeros(2,2)}, [obj.numSeqs, obj.seqLength]);
            obj.transMatAcc = nan(length(odrs), length(poss));
            obj.transMatPerf = nan(length(odrs), length(poss),2);
            obj.lagAccVect = zeros(obj.numSeqs,length(obj.lagVect),2);
            obj.lagAccVectSFP = zeros(obj.numSeqs,length(obj.lagVect),2);
            obj.transMatLatRaw = cell(length(odrs), length(poss), 2);
            obj.lagLatVectRaw = cell(obj.numSeqs,length(obj.lagVect),2);
            obj.lagLatVectSFPraw = cell(obj.numSeqs,length(obj.lagVect),2);
            % Run through session file to calculate things
            for odr = 1:length(odrs)
                odrTrlLog = [obj.trialInfo.Odor]==odrs(odr);
                curSeq = find(sum(obj.odrSeqs==odrs(odr),2));                
                curOdrNdx = find(odrs(odr)==obj.odrSeqs(curSeq,:));
                for pos = 1:length(poss)
                    posTrlLog = [obj.trialInfo.Position]==pos;
                    if curOdrNdx == pos
                        obj.isTrialNums(curSeq,curOdrNdx,1) = sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.isTrialNums(curSeq,curOdrNdx,2) = sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        % Correct Trials
                        obj.responseMatrix(1,1,curSeq) = obj.responseMatrix(1,1,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.responseMatrixByPos{curSeq,pos}(1,1) = obj.responseMatrixByPos{curSeq,pos}(1,1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.responseMatrixByOdr{curSeq,curOdrNdx}(1,1) = obj.responseMatrixByOdr{curSeq,curOdrNdx}(1,1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        % Incorrect Trials
                        obj.responseMatrix(1,2,curSeq) = obj.responseMatrix(1,2,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        obj.responseMatrixByPos{curSeq,pos}(1,2) = obj.responseMatrixByPos{curSeq,pos}(1,2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        obj.responseMatrixByOdr{curSeq,curOdrNdx}(1,2) = obj.responseMatrixByOdr{curSeq,curOdrNdx}(1,2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                    elseif find(odrs(odr)==obj.odrSeqs(curSeq,:)) ~= pos
                        % Correct Trials
                        obj.responseMatrix(2,2,curSeq) = obj.responseMatrix(2,2,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.responseMatrixByPos{curSeq,pos}(2,2) = obj.responseMatrixByPos{curSeq,pos}(2,2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.responseMatrixByOdr{curSeq,curOdrNdx}(2,2) = obj.responseMatrixByOdr{curSeq,curOdrNdx}(2,2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        % Incorrect Trials
                        obj.responseMatrix(2,1,curSeq) = obj.responseMatrix(2,1,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        obj.responseMatrixByPos{curSeq,pos}(2,1) = obj.responseMatrixByPos{curSeq,pos}(2,1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        obj.responseMatrixByOdr{curSeq,curOdrNdx}(2,1) = obj.responseMatrixByOdr{curSeq,curOdrNdx}(2,1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                    end
                    if pos ~= 1
                        if curOdrNdx == pos
                            obj.responseMatrixSFP(1,1,curSeq) = obj.responseMatrixSFP(1,1,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                            obj.responseMatrixSFP(1,2,curSeq) = obj.responseMatrixSFP(1,2,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        elseif curOdrNdx ~= pos
                            obj.responseMatrixSFP(2,2,curSeq) = obj.responseMatrixSFP(2,2,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                            obj.responseMatrixSFP(2,1,curSeq) = obj.responseMatrixSFP(2,1,curSeq) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        end
                        obj.lagAccVectSFP(curSeq,obj.lagVect==(pos-curOdrNdx),1) = obj.lagAccVectSFP(curSeq,obj.lagVect==(pos-curOdrNdx),1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                        obj.lagAccVectSFP(curSeq,obj.lagVect==(pos-curOdrNdx),2) = obj.lagAccVectSFP(curSeq,obj.lagVect==(pos-curOdrNdx),2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                        
                        obj.lagLatVectSFPraw{curSeq,obj.lagVect==(pos-curOdrNdx),1} = [obj.lagLatVectSFPraw{curSeq,obj.lagVect==(pos-curOdrNdx),1}; [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==1).PokeDuration]'];
                        obj.lagLatVectSFPraw{curSeq,obj.lagVect==(pos-curOdrNdx),2} = [obj.lagLatVectSFPraw{curSeq,obj.lagVect==(pos-curOdrNdx),2}; [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==0).PokeDuration]'];
                    end
                    
                    obj.transMatAcc(odr,pos) = mean([obj.trialInfo(odrTrlLog & posTrlLog).Performance]);
                    obj.transMatPerf(odr,pos,1) = sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                    obj.transMatPerf(odr,pos,2) = sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);
                    obj.lagAccVect(curSeq,obj.lagVect==(pos-curOdrNdx),1) = obj.lagAccVect(curSeq,obj.lagVect==(pos-curOdrNdx),1) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==1);
                    obj.lagAccVect(curSeq,obj.lagVect==(pos-curOdrNdx),2) = obj.lagAccVect(curSeq,obj.lagVect==(pos-curOdrNdx),2) + sum([obj.trialInfo(odrTrlLog & posTrlLog).Performance]==0);                  
                    
                    obj.transMatLatRaw{odr,pos,1} = [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==1).PokeDuration]';
                    obj.transMatLatRaw{odr,pos,2} = [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==0).PokeDuration]';
                    obj.lagLatVectRaw{curSeq,obj.lagVect==(pos-curOdrNdx),1} = [obj.lagLatVectRaw{curSeq,obj.lagVect==(pos-curOdrNdx),1}; [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==1).PokeDuration]'];
                    obj.lagLatVectRaw{curSeq,obj.lagVect==(pos-curOdrNdx),2} = [obj.lagLatVectRaw{curSeq,obj.lagVect==(pos-curOdrNdx),2}; [obj.trialInfo(odrTrlLog & posTrlLog & [obj.trialInfo.Performance]==0).PokeDuration]'];
                end
            end
            obj.transMatLatMean = cellfun(@(a)mean(a, 'omitnan'), obj.transMatLatRaw, 'uniformoutput', 0);
            obj.lagLatVectMean = cellfun(@(a)mean(a, 'omitnan'), obj.lagLatVectRaw, 'uniformoutput', 0);
            obj.lagLatVectSFPmean = cellfun(@(a)mean(a, 'omitnan'), obj.lagLatVectSFPraw, 'uniformoutput', 0);
            % d' Accuracy
            [obj.dPrime,~,~] = obj.CalculateDprime(obj.responseMatrix);
            [obj.dPrimeSFP,~,~] = obj.CalculateDprime(obj.responseMatrixSFP);
            % smi Accuracy
            obj.smi = obj.CalculateSMI(obj.responseMatrix);
            obj.smiSFP = obj.CalculateSMI(obj.responseMatrixSFP);
            % ri Response Bias
            [obj.ri,~,~] = obj.CalculateRI(obj.responseMatrix);
            [obj.riSFP,~,~] = obj.CalculateRI(obj.responseMatrixSFP);
            % Calculate pos/odr values
            obj.dPrimeByPos = nan(size(obj.odrSeqs));
            obj.dPrimeByOdr = nan(size(obj.odrSeqs));
            obj.smiByPos = nan(size(obj.odrSeqs));
            obj.smiByOdr = nan(size(obj.odrSeqs));
            obj.riByPos = nan(size(obj.odrSeqs));
            obj.riByOdr = nan(size(obj.odrSeqs));
            for seq = 1:obj.numSeqs
                for op = 1:obj.seqLength
                    [obj.dPrimeByPos(seq,op),~,~] = obj.CalculateDprime(obj.responseMatrixByPos{seq,op});
                    [obj.dPrimeByOdr(seq,op),~,~] = obj.CalculateDprime(obj.responseMatrixByOdr{seq,op});
                    obj.smiByPos(seq,op) = obj.CalculateSMI(obj.responseMatrixByPos{seq,op});
                    obj.smiByOdr(seq,op) = obj.CalculateSMI(obj.responseMatrixByOdr{seq,op});
                    [obj.riByPos(seq,op),~,~] = obj.CalculateRI(obj.responseMatrixByPos{seq,op});
                    [obj.riByOdr(seq,op),~,~] = obj.CalculateRI(obj.responseMatrixByOdr{seq,op});
                end
            end
            % Evaluate Trial After OutSeq (TAO) Trials Directly
            taoTrls = obj.trialInfo(find([obj.trialInfo.TranspositionDistance]~=0 & [obj.trialInfo.Performance]==1 & [obj.trialInfo.Position]~=obj.seqLength & (1:length(obj.trialInfo))~=length(obj.trialInfo))+1);
            obj.taoAcc = mean([taoTrls.Performance]);
            obj.taoLatRaw = [{[taoTrls([taoTrls.Performance]==1).PokeDuration]'}, {[taoTrls([taoTrls.Performance]==0).PokeDuration]'}];
            obj.taoLatMean = [mean([taoTrls([taoTrls.Performance]==1).PokeDuration]), mean([taoTrls([taoTrls.Performance]==0).PokeDuration])];
        end
        %% Calculate dPrime
        function [dPrm, h, fa] = CalculateDprime(~, responseMatrix)
            %%
            dPrm = nan(1,size(responseMatrix,3));
            h = nan(1,size(responseMatrix,3));
            fa = nan(1,size(responseMatrix,3));
            for grp = 1:size(responseMatrix,3)
                if sum(responseMatrix(1,:,grp)) == 1 || sum(responseMatrix(2,:,grp)) == 1
                    dPrm(grp) = nan;
                    h(grp) = nan;
                    fa(grp) = nan;
                else
                    h(grp) = responseMatrix(1,1,grp)/sum(responseMatrix(1,:,grp));
                    if h(grp) == 1
                        h(grp) = (sum(responseMatrix(1,:,grp))-1)/sum(responseMatrix(1,:,grp));
                    elseif h(grp) == 0
                        h(grp) = 1/sum(responseMatrix(1,:,grp));
                    end
                    fa(grp) = responseMatrix(2,1,grp)/sum(responseMatrix(2,:,grp));
                    if fa(grp) == 1
                        fa(grp) = (sum(responseMatrix(2,:,grp))-1)/sum(responseMatrix(2,:,grp));
                    elseif fa(grp) == 0
                        fa(grp) = 1/sum(responseMatrix(2,:,grp));
                    end
                    
                    dPrm(grp) = norminv(h(grp))-norminv(fa(grp));
                end
            end
        end
        %% Calculate SMI
        function [smi] = CalculateSMI(~, responseMatrix)
            smi = nan(1,size(responseMatrix,3));
            for grp = 1:size(responseMatrix,3)
                numerator = (responseMatrix(1,1,grp)*responseMatrix(2,2,grp))-(responseMatrix(2,1,grp)*responseMatrix(1,2,grp));
                denominator = sqrt((responseMatrix(1,1,grp)+responseMatrix(1,2,grp))*(responseMatrix(1,1,grp)+responseMatrix(2,1,grp))*(responseMatrix(2,1,grp)+responseMatrix(2,2,grp))*(responseMatrix(1,2,grp)+responseMatrix(2,2,grp)));
                smi(grp) = numerator/denominator;
            end
        end
        %% Calculate RI
        function [ri, h, fa] = CalculateRI(~, responseMatrix)
            % RI calculation taken from:
            % Talwar & Gerstein (2001) J Neurophysiol. 1555-1572
            ri = nan(1,size(responseMatrix,3));
            h = nan(1,size(responseMatrix,3));
            fa = nan(1,size(responseMatrix,3));
            for grp = 1:size(responseMatrix,3)
                if sum(responseMatrix(1,:,grp)) == 1 || sum(responseMatrix(2,:,grp)) == 1
                    ri(grp) = nan;
                    h(grp) = nan;
                    fa(grp) = nan;
                else
                    h(grp) = responseMatrix(1,1,grp)/sum(responseMatrix(1,:,grp));
                    if h(grp) == 1
                        h(grp) = (sum(responseMatrix(1,:,grp))-1)/sum(responseMatrix(1,:,grp));
                    elseif h(grp) == 0
                        h(grp) = 1/sum(responseMatrix(1,:,grp));
                    end
                    fa(grp) = responseMatrix(2,1,grp)/sum(responseMatrix(2,:,grp));
                    if fa(grp) == 1
                        fa(grp) = (sum(responseMatrix(2,:,grp))-1)/sum(responseMatrix(2,:,grp));
                    elseif fa(grp) == 0
                        fa(grp) = 1/sum(responseMatrix(2,:,grp));
                    end
                    
                    ri(grp) = (h(grp) + fa(grp) - 1) / (1 - (h(grp) - fa(grp))^2);
                end
            end
        end
    end
    %% Standard Plots
    methods
        %% Line Plot w/Mean, SEM & CI
        function [plt, mptValLog] = PlotMeanVarLine(obj,timeVect,data,repDim,pCrit,color)
            dtaMean = mean(data, repDim, 'omitnan');
            if size(data,repDim)>=2
                %             dtaMean = median(data, repDim, 'omitnan');
                dtaMean = dtaMean(:);
                dtaSEM = obj.SEMcalc(data,0,repDim);
                dtaSEM = dtaSEM(:);
                mptValLog = isnan(dtaMean) | isnan(dtaSEM);
                dtaMean(mptValLog) = [];
                dtaSEM(mptValLog) = [];
                if ~exist('pCrit','var')
                    dtaCI = tinv(0.975, size(data,repDim)-1).*dtaSEM;
                else
                    dtaCI = tinv(1-(pCrit/2), size(data,repDim)-1).*dtaSEM;
                end
                if ~exist('color','var')
                    color = 'k';
                end
                timeVect(mptValLog) = [];
                patch('XData', [timeVect(:); flipud(timeVect(:))],...
                    'YData', [(dtaMean+dtaSEM)', flipud(dtaMean-dtaSEM)'],...
                    'linestyle', 'none', 'facecolor', color, 'facealpha', 0.25);
                hold on;
            end
            plt = plot(timeVect, dtaMean, 'color', color, 'linewidth', 1);
            if pCrit ~= 0 
                patch('XData', [timeVect(:); flipud(timeVect(:))],...
                    'YData', [(dtaMean+dtaCI)', flipud(dtaMean-dtaCI)'],...
                    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', color, 'facealpha', 0);
            end
        end
        %% Bar & Swarm Plot w/Mean, SEM (and CI if updated)
        function [plt] = PlotMeanVarSwarmBar(obj,xVal,data,repDim,pCrit,color,varargin)
            dtaMean = mean(data, repDim, 'omitnan');
            dtaSEM = obj.SEMcalc(data,0,repDim);
            if ~exist('pCrit','var')
                dtaCI = tinv(0.975, size(data,repDim)-1).*dtaSEM;
            else
                dtaCI = tinv(1-(pCrit/2), size(data,repDim)-1).*dtaSEM;
            end
            if sum(strcmp(varargin, 'filled'))>=1
                plt = bar(xVal, dtaMean, 'facecolor', color, 'facealpha', 0.4);
            else
                plt = bar(xVal, dtaMean, 'facecolor', 'none');
            end
            hold on;
            swarmchart(zeros(numel(data),1)+xVal,data, 20, 'markerfacecolor', color,...
                'markeredgecolor', 'none', 'markerfacealpha', 0.1);
            if sum(strcmp(varargin, 'error'))>=1
                if strcmp(varargin{find(strcmp(varargin, 'error'))+1}, 'CI')           
                    errorbar(xVal, dtaMean, dtaCI, dtaCI, 'color', 'k', 'capsize', 0);
                elseif strcmp(varargin{find(strcmp(varargin, 'error'))+1}, 'SEM')           
                    errorbar(xVal, dtaMean, dtaSEM, dtaSEM, 'color', 'k', 'capsize', 0);
                end
            else
                errorbar(xVal, dtaMean, dtaSEM, dtaSEM, 'color', 'k', 'capsize', 0);
            end                            
        end
        %% Bar & Distribution plots
        function [plt] = PlotMeanVarViolin(obj,xVal,data,repDim,pCrit,color,varargin)
            if sum(strcmp(varargin, 'binNum'))==1
                numBins = varargin{find(strcmp(varargin, 'binNum'))+1};
            else
                numBins = 100;
            end
            if sum(strcmp(varargin, 'median'))>=1
                dtaMean = median(data, repDim, 'omitnan');
            else
                dtaMean = mean(data, repDim, 'omitnan');
            end
            dtaSEM = obj.SEMcalc(data,0,repDim);
             if ~exist('pCrit','var')
                dtaCI = tinv(0.975, size(data,repDim)-1).*dtaSEM;
            else
                dtaCI = tinv(1-(pCrit/2), size(data,repDim)-1).*dtaSEM;
            end
            if sum(strcmp(varargin, 'filled'))>=1
                plt = bar(xVal, dtaMean, 'edgecolor', color, 'facecolor', color, 'facealpha', 0.4);
            else
                plt = bar(xVal, dtaMean, 'edgecolor', color,  'facecolor', 'none');
            end
            hold on;
            if sum(isnan(data))~=numel(data)
                [counts,edges] = histcounts(data, linspace(min(data), max(data),numBins));
                hist = smooth(smooth(counts,ceil(numBins*.2)),ceil(numBins*.2));
                hist = hist./max(hist);
                binCenters = (edges(2:end)-mode(diff(edges))/2)';
                patch('XData', (([hist; flipud(hist.*-1)]).*0.45)+xVal,...
                    'YData', [binCenters;flipud(binCenters)],...
                    'linestyle','none', 'facecolor', color, 'facealpha', 0.25);
            end
            if sum(strcmp(varargin, 'error'))>=1
                if strcmp(varargin{find(strcmp(varargin, 'error'))+1}, 'CI')           
                    errorbar(xVal, dtaMean, dtaCI, dtaCI, 'color', color, 'capsize', 0);
                elseif strcmp(varargin{find(strcmp(varargin, 'error'))+1}, 'SEM')           
                    errorbar(xVal, dtaMean, dtaSEM, dtaSEM, 'color', color, 'capsize', 0);
                end
            else
                errorbar(xVal, dtaMean, dtaSEM, dtaSEM, 'color', 'k', 'capsize', 0);
            end
        end
        %% Correlation Scatter Plot
        function [rVal,pVal,scatPlot] = CorrScatPlot(~, xDta, yDta, xLbl, yLbl, varargin)
            % [rVal, pVal] = corrScatPlot(xDta, yDta, xLbl, yLbl, varargin);
            % xDta = X-data
            % yDta = Y-data
            % xLbl = X-axis Label
            % yLbl = Y-axis Label
            %%
            if nargin == 2
                xLbl = [];
                yLbl = [];
            end
            if strcmp(varargin, 'zroNaN')
                xNanLog = isnan(xDta);
                yNanLog = isnan(yDta);
                xDta(xNanLog | yNanLog) = 0;
                yDta(xNanLog | yNanLog) = 0;
            elseif strcmp(varargin, 'rmvNaN')
                xNanLog = isnan(xDta);
                yNanLog = isnan(yDta);
                xDta(xNanLog | yNanLog) = [];
                yDta(xNanLog | yNanLog) = [];
            elseif strcmp(varargin, 'rmvZro')
                zroLog = xDta==0 | yDta==0;
                xDta(zroLog) = [];
                yDta(zroLog) = [];
            end
            [r, p] = corr(xDta, yDta, 'rows', 'pairwise');
            scatPlot = scatter(xDta, yDta);
            if sum(strcmp(varargin, 'markerStyle'))>=1
                mrkrSpt = find(strcmp(varargin, 'markerStyle'));
                set(scatPlot, 'Marker', varargin{mrkrSpt+1});
            end
            if sum(strcmp(varargin, 'markerColor'))>=1
                clrSpt = find(strcmp(varargin, 'markerColor'));
                set(scatPlot, 'MarkerEdgeColor', varargin{clrSpt+1});
            end
            hold on
            data = [xDta,yDta];
            data(sum(isnan(data),2)>=1,:) = [];
            eq1 = polyfit(data(:,1),data(:,2),1);
            unqX = unique(data(:,1));
            vCalc1 = polyval(eq1, unqX);
            plot(unqX, vCalc1, 'red', 'linewidth', 2);
            if p<=0.05
                annotation(gcf, 'textbox', get(gca, 'position'), 'string',...
                    sprintf('r=%.02f; p = %.02i', r, p), 'linestyle', 'none', 'color', 'r');
            else
                annotation(gcf, 'textbox', get(gca, 'position'), 'string',...
                    sprintf('r=%.02f; p = %.03f', r, p), 'linestyle', 'none');
            end
            xlabel(xLbl);
            ylabel(yLbl);

            rVal = r;
            pVal = p;
            if sum(strcmp(varargin, 'coordOn'))==1
                hold on
                line(get(gca, 'xlim'), [0 0], 'linewidth', 0.5, 'linestyle', ':', 'color', 'black')
                line([0 0], get(gca, 'ylim'), 'linewidth', 0.5, 'linestyle', ':', 'color', 'black')
            end
            if sum(strcmp(varargin, 'gridOn'))==1
                grid on;
            end
        end
    end
    %% Misc
    methods
        %% SEM calculation
        function [semVal] = SEMcalc(~,data, nVal, dim)
            %% SEMcalc Calculates SEM
            % data = data set where SEM is being calculated
            % nVal = determination if denominator is n or n-1
            
            %%
            if isempty(data)
                semVal = nan;
                return
            end
            if nargin==1
                nVal = 0;
                dim = 1;
            elseif nargin==2
                dim = 1;
            end
            
            denom = sum(~isnan(data),dim);
            denom(denom==0) = 1;
            semVal = std(data,nVal,dim, 'omitnan')./sqrt(denom-1);
        end
    end
end