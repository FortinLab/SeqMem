classdef PFC_UniSum_MLB_SM < MLB_SM
    %% Properties
    properties % Analysis parameters
        preErlyWindow = [-500 500]
        ltPstWindow = [-500 500]
        rwdWindow = [-500 500]
        errWindow = [-500 500]
        trialWindow = [-500 1500]
    end
    properties % Data Matrices
        wholeTrialMtx
        wholeTrialTimeVect
        
        preErlyMtx
        preErlyTime
        
        ltPstMtx
        ltPstTime
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
                obj.ExtractTrialMatrix;
                obj.ExtractTrialPeriods
            end
        end
    end
    %% Pre-Processing Methods
    methods
        %% Extract FIS 
        function ExtractTrialMatrix(obj)
            [obj.wholeTrialMtx, obj.wholeTrialTimeVect] = obj.PP_TrialMatrix(obj.trialWindow, 'PokeIn');            
        end
        %% Extract Poke In
        function ExtractTrialPeriods(obj)
            [obj.preErlyMtx, obj.preErlyTime] = obj.PP_TrialMatrix(obj.preErlyWindow, 'PokeIn');
            [obj.ltPstMtx, obj.ltPstTime] = obj.PP_TrialMatrix(obj.preErlyWindow, 'PokeOut');
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