classdef PFC_UniSum_MLB_SM < MLB_SM
    %% Properties
    properties % Analysis parameters
        preTrialWindow = [-500 0]
        erlyTrialWindow = [0 500]
        lateTrialWindow = [-500 0]
        pstTrialWindow = [0 500]
        preRwdWindow = [-500 0]
        pstRwdWindow = [0 500]
        preErrWindow = [-500 0]
        pstErrWindow = [0 500]
        trialWindow = [-500 1500]
    end
    properties % Categorical Descriptors
        maxResponsePeriod
    end
    properties % Data Matrices
        wholeTrialMtx
        wholeTrialTimeVect
        
        preTrialMtx
        preTrialTime
        
        erlyTrialMtx
        erlyTrialTime
        
        lateTrialMtx
        lateTrialTime
        
        pstTrialMtx
        pstTrialTime   
        
        preRwdMtx
        preRwdTime 
        
        pstRwdMtx
        pstRwdTime
        
        preErrMtx
        preErrTime
        
        pstErrMtx
        pstErrTime
    end
    properties % Mean Firing Rate Values
        preTrialFR
        erlyTrialFR
        lateTrialFR
        pstTrialFR
        preRwdFR
        pstRwdFR
        preErrFR
        pstErrFR
        allTrialPeriodsFR
        ensembleTrialPeriodMeanFR
    end
    properties % Summary Stats
        epochCorrRhoMatrixC
        epochCorrRhoMatrixIC
        epochCorrSigMatrixC
        epochCorrSigMatrixIC
        
        trialEpochFc = cell(4,7)
        trialEpochFic = cell(4,7)
        
        trialEpochPerfF = cell(6,7)
        
        trialEpochOdrPosPerfF = cell(18,7)
        
        preTrialPosFc = cell(4,7)
        preTrialPosFic = cell(4,7)
        preTrialOdrFc = cell(4,7)
        preTrialOdrFic = cell(4,7)
        
        erlyTrialPosFc = cell(4,7)
        erlyTrialPosFic = cell(4,7)
        erlyTrialOdrFc = cell(4,7)
        erlyTrialOdrFic = cell(4,7)
        
        lateTrialPosFc = cell(4,7)
        lateTrialPosFic = cell(4,7)
        lateTrialOdrFc = cell(4,7)
        lateTrialOdrFic = cell(4,7)
        
        pstTrialPosFc = cell(4,7)
        pstTrialPosFic = cell(4,7)
        pstTrialOdrFc = cell(4,7)
        pstTrialOdrFic = cell(4,7)
        
        preRwdPosFc = cell(4,7)
        preRwdOdrFc = cell(4,7)
        
        pstRwdPosFc = cell(4,7)
        pstRwdOdrFc = cell(4,7)
        
        preErrPosFic = cell(4,7)
        preErrOdrFic = cell(4,7)
        
        pstErrPosFc = cell(4,7)
        pstErrPosFic = cell(4,7)
        pstErrOdrFc = cell(4,7)
        pstErrOdrFic = cell(4,7)
    end
        
    %% Methods
    methods
        %% Object Creation
        %#ok<*EXIST>
        function obj = PFC_UniSum_MLB_SM(fileDir, binSize, dsRate)
            if nargin == 0 || isempty(fileDir)
                fileDir = uigetdir;
            end
            obj@MLB_SM(fileDir);
            if nargin == 3
                obj.binSize = binSize;
                obj.dsRate = dsRate;
                obj.ExtractTrialMatrices;
                obj.ExtractTrialPeriods;
                obj.CalculateEpochCorrMtx;
                obj.TrialPeriodModulation;
            end
        end
    end
    %% Pre-Processing Methods
    methods
        %% Extract FIS 
        function ExtractTrialMatrices(obj)
            [obj.wholeTrialMtx, obj.wholeTrialTimeVect] = obj.PP_TrialMatrix(obj.trialWindow, 'PokeIn');            
        end
        %% Extract Poke In
        function ExtractTrialPeriods(obj)
            [obj.preTrialMtx, obj.preTrialTime] = obj.PP_TrialMatrix(obj.preTrialWindow, 'PokeIn');
            [obj.erlyTrialMtx, obj.erlyTrialTime] = obj.PP_TrialMatrix(obj.erlyTrialWindow, 'PokeIn');
            [obj.lateTrialMtx, obj.lateTrialTime] = obj.PP_TrialMatrix(obj.lateTrialWindow, 'PokeOut');
            [obj.pstTrialMtx, obj.pstTrialTime] = obj.PP_TrialMatrix(obj.pstTrialWindow, 'PokeOut');
            [obj.preRwdMtx, obj.preRwdTime] = obj.PP_TrialMatrix(obj.preRwdWindow, 'FrontReward');
            [obj.pstRwdMtx, obj.pstRwdTime] = obj.PP_TrialMatrix(obj.pstRwdWindow, 'FrontReward');
            [obj.preErrMtx, obj.preErrTime] = obj.PP_TrialMatrix(obj.preErrWindow, 'ErrorSignal');
            [obj.pstErrMtx, obj.pstErrTime] = obj.PP_TrialMatrix(obj.pstErrWindow, 'ErrorSignal');
            
            obj.preTrialFR = reshape(mean(obj.preTrialMtx,1), [size(obj.preTrialMtx,2),size(obj.preTrialMtx,3)])';
            obj.erlyTrialFR = reshape(mean(obj.erlyTrialMtx,1), [size(obj.erlyTrialMtx,2), size(obj.erlyTrialMtx,3)])';
            obj.lateTrialFR = reshape(mean(obj.lateTrialMtx,1), [size(obj.lateTrialMtx,2), size(obj.lateTrialMtx,3)])';
            obj.pstTrialFR = reshape(mean(obj.pstTrialMtx,1), [size(obj.pstTrialMtx,2), size(obj.pstTrialMtx,3)])';
            obj.preRwdFR = reshape(mean(obj.preRwdMtx,1), [size(obj.preRwdMtx,2), size(obj.preRwdMtx,3)])';
            obj.pstRwdFR = reshape(mean(obj.pstRwdMtx,1), [size(obj.pstRwdMtx,2), size(obj.pstRwdMtx,3)])';
            obj.preErrFR = reshape(mean(obj.preErrMtx,1), [size(obj.preErrMtx,2), size(obj.preErrMtx,3)])';
            obj.pstErrFR = reshape(mean(obj.pstErrMtx,1), [size(obj.pstErrMtx,2), size(obj.pstErrMtx,3)])';
            
            obj.allTrialPeriodsFR(:,1,:) = obj.preTrialFR;
            obj.allTrialPeriodsFR(:,2,:) = obj.erlyTrialFR;
            obj.allTrialPeriodsFR(:,3,:) = obj.lateTrialFR;
            obj.allTrialPeriodsFR(:,4,:) = obj.pstTrialFR;      
            obj.allTrialPeriodsFR(:,5,:) = obj.preRwdFR;
            obj.allTrialPeriodsFR(:,6,:) = obj.pstRwdFR;
            obj.allTrialPeriodsFR(:,7,:) = obj.preErrFR;
            obj.allTrialPeriodsFR(:,8,:) = obj.pstErrFR;
            
            obj.ensembleTrialPeriodMeanFR = reshape(nanmean(obj.allTrialPeriodsFR), [size(obj.allTrialPeriodsFR,2), size(obj.allTrialPeriodsFR,3)])';
            obj.maxResponsePeriod = nan(1,length(obj.ensembleMatrixColIDs));
            for tet = 1:length(obj.ensembleMatrixColIDs)
                obj.maxResponsePeriod(tet) = find(max(obj.ensembleTrialPeriodMeanFR(tet,1:4))==obj.ensembleTrialPeriodMeanFR(tet,1:4),1,'first');
            end
        end
    end
    %% Analysis Methods
    methods 
        %% Analyze Ensemble
        % Step through and run each analysis on individual cells
        function CalculateEpochCorrMtx(obj)
            obj.epochCorrRhoMatrixC = repmat(eye(4), [1,1,length(obj.ensembleMatrixColIDs)]);
            obj.epochCorrRhoMatrixIC = repmat(eye(4), [1,1,length(obj.ensembleMatrixColIDs)]);
            obj.epochCorrSigMatrixC = nan(4,4,length(obj.ensembleMatrixColIDs));
            obj.epochCorrSigMatrixIC = nan(4,4,length(obj.ensembleMatrixColIDs));
            perfLog = [obj.trialInfo.Performance]==1;
            for tet = 1:length(obj.ensembleMatrixColIDs)
                for prd1 = 1:size(obj.allTrialPeriodsFR,2)
                    for prd2 = 1:size(obj.allTrialPeriodsFR,2)
                        [obj.epochCorrRhoMatrixC(prd2,prd1,tet), obj.epochCorrSigMatrixC(prd2,prd1,tet)] = corr(obj.allTrialPeriodsFR(perfLog,prd1,tet),...
                            obj.allTrialPeriodsFR(perfLog,prd2,tet), 'rows', 'complete');
                        [obj.epochCorrRhoMatrixIC(prd2,prd1,tet), obj.epochCorrSigMatrixIC(prd2,prd1,tet)] = corr(obj.allTrialPeriodsFR(~perfLog,prd1,tet),...
                            obj.allTrialPeriodsFR(~perfLog,prd2,tet), 'rows', 'complete');
                    end
                end                    
%                 figure; subplot(1,2,1); imagesc(obj.epochCorrRhoMatrixC(:,:,tet), [0 1]); title(obj.ensembleMatrixColIDs{tet}); subplot(1,2,2); imagesc(obj.epochCorrRhoMatrixIC(:,:,tet));
            end            
        end
        %% Trial Period Modulation Analysis
        function TrialPeriodModulation(obj)
            odorLog = [obj.trialInfo.Odor];
            posLog = [obj.trialInfo.Position];
            perfLog = [obj.trialInfo.Performance];
            for tet = 1:length(obj.ensembleMatrixColIDs)
                % Do cells modulate their firing rate based on trial period?
                tempEpochFRc = obj.allTrialPeriodsFR(perfLog,1:4,tet);
                tempEpochLogC = ones(sum(perfLog),4).*repmat(1:4, [sum(perfLog),1]);
                [~, obj.trialEpochFc(:,:,tet)] = anovan(tempEpochFRc(:), tempEpochLogC(:), 'display', 'off');
                tempEpochFRic = obj.allTrialPeriodsFR(~perfLog,1:4,tet);
                tempEpochLogIC = ones(sum(~perfLog),4).*repmat(1:4, [sum(~perfLog),1]);
                [~, obj.trialEpochFic(:,:,tet)] = anovan(tempEpochFRic(:), tempEpochLogIC(:), 'display', 'off');
                
                % Do cells modulate their firing rate based on trial period and outcome?
                tempEpochAll = obj.allTrialPeriodsFR(:,1:4,tet);
                tempEpochAllLog = ones(length(perfLog),4).*repmat(1:4, [length(perfLog),1]);
                [~, obj.trialEpochPerfF(:,:,tet)] = anovan(tempEpochAll(:), [tempEpochAllLog(:), repmat(perfLog', 4,1)],...
                    'model', 'interaction', 'varnames', {'Epoch', 'Performance'}, 'display', 'off');
                
                % Do cells modulate their firing rate during trial periods based on trial position or odor?
                [~,obj.preTrialPosFc(:,:,tet)] = anovan(obj.preTrialFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.preTrialPosFic(:,:,tet)] = anovan(obj.preTrialFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.preTrialOdrFc(:,:,tet)] = anovan(obj.preTrialFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                [~,obj.preTrialOdrFic(:,:,tet)] = anovan(obj.preTrialFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                [~,obj.erlyTrialPosFc(:,:,tet)] = anovan(obj.erlyTrialFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.erlyTrialPosFic(:,:,tet)] = anovan(obj.erlyTrialFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.erlyTrialOdrFc(:,:,tet)] = anovan(obj.erlyTrialFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                [~,obj.erlyTrialOdrFic(:,:,tet)] = anovan(obj.erlyTrialFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                [~,obj.lateTrialPosFc(:,:,tet)] = anovan(obj.lateTrialFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.lateTrialPosFic(:,:,tet)] = anovan(obj.lateTrialFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.lateTrialOdrFc(:,:,tet)] = anovan(obj.lateTrialFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                [~,obj.lateTrialOdrFic(:,:,tet)] = anovan(obj.lateTrialFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                [~,obj.pstTrialPosFc(:,:,tet)] = anovan(obj.pstTrialFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.pstTrialPosFic(:,:,tet)] = anovan(obj.pstTrialFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.pstTrialOdrFc(:,:,tet)] = anovan(obj.pstTrialFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                [~,obj.pstTrialOdrFic(:,:,tet)] = anovan(obj.pstTrialFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                [~,obj.preRwdPosFc(:,:,tet)] = anovan(obj.preRwdFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.preRwdOdrFc(:,:,tet)] = anovan(obj.preRwdFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                
                [~,obj.pstRwdPosFc(:,:,tet)] = anovan(obj.pstRwdFR(perfLog,tet), posLog(perfLog)', 'display', 'off');
                [~,obj.pstRwdOdrFc(:,:,tet)] = anovan(obj.pstRwdFR(perfLog,tet), odorLog(perfLog)', 'display', 'off');
                
                [~,obj.preErrPosFic(:,:,tet)] = anovan(obj.preErrFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.preErrOdrFic(:,:,tet)] = anovan(obj.preErrFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                [~,obj.pstErrPosFic(:,:,tet)] = anovan(obj.pstErrFR(~perfLog,tet), posLog(~perfLog)', 'display', 'off');
                [~,obj.pstErrOdrFic(:,:,tet)] = anovan(obj.pstErrFR(~perfLog,tet), odorLog(~perfLog)', 'display', 'off');
                
                % Do cells modulate their firing rate based on trial position and outcome?
                [~, obj.trialEpochOdrPosPerfF(:,:,tet)] = anovan(tempEpochAll(:), [tempEpochAllLog(:), repmat(posLog', 4,1), repmat(perfLog', 4,1), repmat(odorLog', 4,1)],...
                    'model', 'full', 'sstype', 2, 'varnames', {'Epoch', 'Position', 'Performance', 'Odor'}, 'display', 'off');
            end
        % Do cells modulate their firing rate during trial periods based on
        %   odor-in-position conjunction?
        end
        %% Instantaneous trial modulation
        %% Convolve Half Gaussian
        function convdData = ConvolveHalfGauss(obj, data)
            convdData = nan(size(data));
            gauss = gausswin(obj.binSize);
            blankingVect = [zeros(floor(obj.binSize/2),1);repmat(1/max(gauss), [obj.binSize-floor(obj.binSize/2),1])];
            gauss = gauss.*blankingVect;
            % Assumes data is in Time x Unit x Trial
            for u = 1:size(data,2)
                for t = 1:size(data,3)
                    convdData(:,u,t) = conv(data(:,u,t), gauss, 'same');
                end
            end
        end
        %% Convolve Full Gaussian
        function convdData = ConvolveFullGauss(obj, data)
            convdData = nan(size(data));
            gauss = gausswin(obj.binSize);
            % Assumes data is in Time x Unit x Trial
            for u = 1:size(data,2)
                for t = 1:size(data,3)
                    convdData(:,u,t) = conv(data(:,u,t), gauss, 'same');
                end
            end
        end
    end
    %% Plotting Methods
    methods
        %% Plot UniSum Summary
        % Create the axes for each subplot to be associated in the final
        % figure
        %% Plot Unit Waveform & Features
        %% Plot Trial Epoch Spiking
        %% Plot Epoch Correlations
        %% Plot Trial Period Modulation Rates
        %% Plot Trial Rasters
        %% Plot Trial Firing Rates
        %% Plot Positional Modulation 
    end
end