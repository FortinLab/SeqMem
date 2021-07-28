classdef PFC_UniSum_MLB_SM < MLB_SM
    %% Properties
    properties % Analysis parameters
        preTrialWindow = [-500 0]
        erlyTrialWindow = [0 500]
        lateTrialWindow = [-500 0]
        pstTrialWindow = [0 500]
        preRwdWindow = [-500 0]
        pstRwdWindow = [0 500]
        errWindow = [250 250]
        trialWindow = [-500 1500]
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
        
        errPrdMtx
        errPrdTime
    end
    properties % Mean Firing Rate Values
        preTrialFR
        erlyTrialFR
        lateTrialFR
        pstTrialFR
        preRwdFR
        pstRwdFR
        errPrdFR
    end
    properties % Summary Stats
        epochCorrMatrix
        preTrialPosF
        erlyTrialPosF
        lateTrialPosF
        pstTrialPosF
        preRwdPosF
        pstRwdPosF
        errPrdPosF
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
                obj.ExtractTrialPeriods
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
            obj.preTrialFR = reshape(mean(obj.preTrialMtx,1), [size(obj.preTrialMtx,2),size(obj.preTrialMtx,3)])';
            [obj.erlyTrialMtx, obj.erlyTrialTime] = obj.PP_TrialMatrix(obj.erlyTrialWindow, 'PokeIn');
            obj.erlyTrialFR = reshape(mean(obj.erlyTrialMtx,1), [size(obj.erlyTrialMtx,2), size(obj.erlyTrialMtx,3)])';
            [obj.lateTrialMtx, obj.lateTrialTime] = obj.PP_TrialMatrix(obj.lateTrialWindow, 'PokeOut');
            obj.lateTrialFR = reshape(mean(obj.lateTrialMtx,1), [size(obj.lateTrialMtx,2), size(obj.lateTrialMtx,3)])';
            [obj.pstTrialMtx, obj.pstTrialTime] = obj.PP_TrialMatrix(obj.pstTrialWindow, 'PokeOut');
            obj.pstTrialFR = reshape(mean(obj.pstTrialMtx,1), [size(obj.pstTrialMtx,2), size(obj.pstTrialMtx,3)])';
            [obj.preRwdMtx, obj.preRwdTime] = obj.PP_TrialMatrix(obj.preRwdWindow, 'FrontReward');
            obj.preRwdFR = reshape(mean(obj.preRwdMtx,1), [size(obj.preRwdMtx,2), size(obj.preRwdMtx,3)])';
            [obj.pstRwdMtx, obj.pstRwdTime] = obj.PP_TrialMatrix(obj.pstRwdWindow, 'FrontReward');
            obj.pstRwdFR = reshape(mean(obj.pstRwdMtx,1), [size(obj.pstRwdMtx,2), size(obj.pstRwdMtx,3)])';
            [obj.errPrdMtx, obj.errPrdTime] = obj.PP_TrialMatrix(obj.pstTrialWindow, 'ErrorSignal');
            obj.errPrdFR = reshape(mean(obj.errPrdMtx,1), [size(obj.errPrdMtx,2), size(obj.errPrdMtx,3)])';
        end
    end
    %% Analysis Methods
    methods 
        %% Analyze Ensemble
        % Step through and run each analysis on individual cells
        %% Trial Period Modulation Analysis
        % Do cells modulate their firing rate based on trial period?
        % Do cells modulate their firing rate during trial periods based on
        %   trial position?
        % Do cells modulate their firing rate during trial periods based on
        %   the presented odor?
        % Do cells modulate their firing rate during trial periods based on
        %   odor-in-position conjunction?
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