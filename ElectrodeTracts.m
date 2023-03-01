classdef ElectrodeTracts < handle
    %% Properties
    properties (SetObservable, AbortSet)
        Subject
        SlideHeight = 25000     % Height of the slide the brain slice was mounted on in microns
        SlidePixelHeight        % Height of the slide the brain slice was mounted on in pixels
        SlidePPM                % Height of a pixel
        SliceWidth  = 40        % Width of the brain slice in microns
        DriveTravelDist = 250   % Distance the tetrode should have traveled in microns with each full rotation of the drive screw
        TractsByRater
        Tracts
        TractLMs
        TractsReconstructed
        DriveTable
        DistanceTable
        APTable

        raterColors = [0.9, 0.65, 0.1;...
            0.1, 0.1, 0.9;...
            0.9, 0.1, 0.1;...
            0.9, 0.1, 0.9];
    end
    %% Creation Method
    methods
        function obj = ElectrodeTracts(varargin)
            addlistener(obj,'SlidePixelHeight','PostSet',@obj.ListenerCalcPPM)
            if nargin == 1
                obj.Subject = varargin{1};
            elseif nargin == 2
                obj.Subject = varargin{1};
                obj.SlidePixelHeight = varargin{2};
            elseif nargin == 3
                obj.Subject = varargin{1};
                obj.SlidePixelHeight = varargin{2};
                obj.SliceWidth = varargin{3};
            end
        end
    end
    %% Import Methods
    methods
        %% Read Tract Annotations doc
        function ReadTractAnnotations(obj,fn,rater)
            % 'readcell' the excel file
            tractR = readcell(fn);
            tractTetNdx = find(cellfun(@(a)sum(ismissing(a))==0,tractR(1,:)));
            tractTetIDs = tractR(1,tractTetNdx);
            for t = 1:length(tractTetNdx)
                tempTractData = tractR(3:end,tractTetNdx(t):tractTetNdx(t)+3);
                misLog = cellfun(@(a)sum(ismissing(a))==1, tempTractData);
                tempTractData(misLog) = {nan};
                notDoubleLog = cellfun(@(a)~isa(a,'double'), tempTractData);
                notDoubleLog(:,end) = false;
                tempTractData(notDoubleLog) = {nan};
                
                nanLog = sum(cellfun(@(a)sum(isnan(a))==1, tempTractData),2)>=1;
                tempTractData(nanLog,:) = [];
                if ~isempty(tempTractData)
                    obj.TractsByRater.(rater).(tractTetIDs{t}) = tempTractData;
                end
            end
        end
        %% Read Driving Records
        function ReadDriveRecord(obj,fn)
            driveR = readcell(fn);
            driveTetIDs = arrayfun(@(a){sprintf('TET%i',a)},[driveR{4:end,1}]);
            driveColLogs = find(cellfun(@(a)sum(ismissing(a))==0,driveR(1,:)));
            driveRecs = [driveR(4:end,:)];

            sessions = driveR(1,driveColLogs)';
            dates = driveR(2,driveColLogs)';
            driveTable = table(sessions,dates);
            for t = 1:length(driveTetIDs)
                tempTetDriveRec = [driveRecs{t,driveColLogs}]';
                tempTetDriveRec(isnan(tempTetDriveRec)) = 0;
                driveTable.(driveTetIDs{t}) = tempTetDriveRec;
            end
            obj.DriveTable = driveTable;
        end
    end
    %% Listening Methods
    methods (Static)
        %% Calculate Pixels Per Micron
        function ListenerCalcPPM(~,evt)
            evt.AffectedObject.SlidePPM = evt.AffectedObject.SlidePixelHeight/evt.AffectedObject.SlideHeight;
        end
    end
    %% Processing Methods
    methods
        %% Compile Tracts Across Raters
        function CompileTracts(obj)
            raters = fieldnames(obj.TractsByRater);
            tets = struct;
            for r = 1:length(raters)
                measuredTets = fieldnames(obj.TractsByRater.(raters{r}));
                for t = 1:length(measuredTets)
                    if isfield(tets, measuredTets{t})
                        tets.(measuredTets{t}) = sortrows([tets.(measuredTets{t}); obj.TractsByRater.(raters{r}).(measuredTets{t})]);
                    else
                        tets.(measuredTets{t}) = obj.TractsByRater.(raters{r}).(measuredTets{t});
                    end
                end
            end
            tetNames = fieldnames(tets);
            for t = 1:length(tetNames)
                curTet = tets.(tetNames{t});
                if size(curTet,1)==length(unique(cell2mat(curTet(:,1))))
                    tempMeasure = curTet;
                else
                    unqTracts = cell(size(curTet));
                    for s = 1:size(curTet,1)
                        sliceLog = curTet{s,1}==[curTet{:,1}];
                        if find(sliceLog,1,'first')==s
                            unqTracts(s,1:end-1) = num2cell(mean(cell2mat(curTet(sliceLog,1:3)),1));
                            if length(unique(curTet(sliceLog,end))) ~= 1
                                error('Mismatched regional attribution');
                            else
                                unqTracts(s,end) = unique(curTet(sliceLog,end));
                            end
                        end
                    end
                    nonMptLog = sum(cellfun(@(a)isempty(a),unqTracts),2)==0;
                    tempMeasure = unqTracts(nonMptLog,:);
                end
                tets.(tetNames{t}) = [num2cell(cell2mat(tempMeasure(:,1))*obj.SliceWidth), num2cell(cell2mat(tempMeasure(:,2:end-1))/obj.SlidePPM), tempMeasure(:,end)];
            end
            obj.Tracts = tets;
        end
        %% Fit tracts to linear models
        function FitTractLM(obj)
            if isempty(obj.Tracts)
                obj.CompileTracts;
            end
            tetIDs = fieldnames(obj.Tracts);
            for t = 1:length(tetIDs)
                tempTetMeasure = cell2mat(obj.Tracts.(tetIDs{t})(:,1:end-1));
                if size(tempTetMeasure,1)>=2
                    % AP Regression
                    [obj.TractLMs.(tetIDs{t}).B(:,1), ~, apR] = regress(tempTetMeasure(:,1), tempTetMeasure(:,2:3));
                    obj.TractLMs.(tetIDs{t}).SS(1) = sum((apR*1000).^2);
                    obj.TractLMs.(tetIDs{t}).ML(:,1) = linspace(max(tempTetMeasure(:,2)), min(tempTetMeasure(:,2)),1000);
                    obj.TractLMs.(tetIDs{t}).DV(:,1) = linspace(min(tempTetMeasure(:,3)), max(tempTetMeasure(:,3)),1000);
                    obj.TractLMs.(tetIDs{t}).AP(:,1) = obj.TractLMs.(tetIDs{t}).ML(:,1)*obj.TractLMs.(tetIDs{t}).B(1,1)...
                        + obj.TractLMs.(tetIDs{t}).DV(:,1)*obj.TractLMs.(tetIDs{t}).B(2,1);
                    % DV Regression
                    [obj.TractLMs.(tetIDs{t}).B(:,2), ~, dvR] = regress(tempTetMeasure(:,3), tempTetMeasure(:,1:2));
                    obj.TractLMs.(tetIDs{t}).SS(2) = sum((dvR*1000).^2);
                    obj.TractLMs.(tetIDs{t}).AP(:,2) = linspace(min(tempTetMeasure(:,1)), max(tempTetMeasure(:,1)),1000); 
                    obj.TractLMs.(tetIDs{t}).ML(:,2) = linspace(max(tempTetMeasure(:,2)), min(tempTetMeasure(:,2)),1000);
                    obj.TractLMs.(tetIDs{t}).DV(:,2) = obj.TractLMs.(tetIDs{t}).AP(:,2)*obj.TractLMs.(tetIDs{t}).B(1,2)...
                        + obj.TractLMs.(tetIDs{t}).ML(:,2)*obj.TractLMs.(tetIDs{t}).B(2,2);
                    % ML Regression
                    [obj.TractLMs.(tetIDs{t}).B(:,3), ~, mlR] = regress(tempTetMeasure(:,2), tempTetMeasure(:,[1,3]));
                    obj.TractLMs.(tetIDs{t}).SS(3) = sum((mlR*1000).^2);
                    obj.TractLMs.(tetIDs{t}).AP(:,3) = linspace(min(tempTetMeasure(:,1)), max(tempTetMeasure(:,1)),1000); 
                    obj.TractLMs.(tetIDs{t}).DV(:,3) = linspace(min(tempTetMeasure(:,3)), max(tempTetMeasure(:,3)),1000);
                    obj.TractLMs.(tetIDs{t}).ML(:,3) = obj.TractLMs.(tetIDs{t}).AP(:,3)*obj.TractLMs.(tetIDs{t}).B(1,3)...
                        + obj.TractLMs.(tetIDs{t}).DV(:,3)*obj.TractLMs.(tetIDs{t}).B(2,3);
                end
            end
        end
        %% Integrate Driving Record
        function CompileDriving(obj)
            if isempty(obj.TractLMs) || isempty(obj.Tracts)
                obj.FitTractLM;
            end
            if isempty(obj.DriveTable)
                error('No drive history read in');
            end
            tetIDs = fieldnames(obj.TractLMs);
            distTable = obj.DriveTable;
            for t = 1:length(tetIDs)
                distTable.(tetIDs{t}) = cumsum(obj.DriveTable.(tetIDs{t})*obj.DriveTravelDist);
            end
            obj.DistanceTable = distTable;
        end
        %% Reconstruct Tracts
        function ReconstructTracts(obj)
            if isempty(obj.TractLMs) || isempty(obj.Tracts)
                obj.FitTractLM;
            end
            if isempty(obj.DriveTable)
                error('No drive history read in');
            end
            if isempty(obj.DistanceTable)
                obj.CompileDriving;
            end
            tetIDs = fieldnames(obj.TractLMs);
            reconTracts = struct;
            apTable = obj.DistanceTable;
            for t = 1:length(tetIDs)
                % Put in drive travel reconstruction here
                tempDriveDist = obj.DistanceTable.(tetIDs{t});
                lmLog = obj.TractLMs.(tetIDs{t}).SS==min(obj.TractLMs.(tetIDs{t}).SS);
                dvDiff = mean(diff(flipud(obj.TractLMs.(tetIDs{t}).DV(:,lmLog))));
                apDiff = mean(diff(flipud(obj.TractLMs.(tetIDs{t}).AP(:,lmLog))));
                mlDiff = mean(diff(flipud(obj.TractLMs.(tetIDs{t}).ML(:,lmLog))));
                diffPerNdx = sqrt(dvDiff^2+sqrt(apDiff^2+mlDiff^2));
                tractIndexLength = ceil(tempDriveDist(end)/diffPerNdx);
                tractIndexVect = (tractIndexLength:-1:1)';
                dvTract = obj.TractLMs.(tetIDs{t}).DV(end,lmLog)+(dvDiff*(tractIndexVect-1));
                mlTract = obj.TractLMs.(tetIDs{t}).ML(end,lmLog)+(mlDiff*(tractIndexVect-1));
                apTract = obj.TractLMs.(tetIDs{t}).AP(end,lmLog)+(apDiff*(tractIndexVect-1));
                reconTracts.(tetIDs{t}) = [apTract, mlTract, dvTract];
                
                tractLength = flipud(tractIndexVect)*diffPerNdx;
                tempTractAP = nan(size(tempDriveDist));
                for d = 1:length(tempTractAP)
                    lengthIndex = find(tempDriveDist(d)<tractLength,1,'first');
                    tempTractAP(d) = ceil(apTract(lengthIndex)/obj.SliceWidth);
                end
                apTable.(tetIDs{t}) = tempTractAP;
            end
            obj.APTable = apTable;
            obj.TractsReconstructed = reconTracts;
        end
    end
    %% Visualization Methods
    methods
        %% Plot Tracts By Rater
        function PlotTractsByRater(obj)
            figure;
            raters = fieldnames(obj.TractsByRater);
            rateLegs = nan(1,length(raters));
            for r = 1:length(raters)
                measuredTets = fieldnames(obj.TractsByRater.(raters{r}));
                for t = 1:length(measuredTets)
                    tempTet = cell2mat(obj.TractsByRater.(raters{r}).(measuredTets{t})(:,1:end-1));
                    rateLegs(r) = plot3(tempTet(:,2),...
                        tempTet(:,1),...
                        tempTet(:,3),...
                        'linewidth', 2, 'color', obj.raterColors(r,:), 'marker', '.', 'markersize', 10);
                    hold on;
                    text(tempTet(1,2),...
                        tempTet(1,1),...
                        tempTet(1,3),...
                        measuredTets{t});
                end
            end
            grid on;
            set(gca, 'ZDir', 'reverse');
            legend(rateLegs, raters)
            xlabel(gca,'ML (pixels)');
            ylabel(gca,'AP (pixels)');
            zlabel(gca,'DV (pixels)');
        end
        %% Plot All Tracts
        function PlotTracts(obj)
            if isempty(obj.TractLMs) || isempty(obj.Tracts)
                obj.FitTractLM;
            end

            tetIDs = fieldnames(obj.Tracts);
            figure;
            for t = 1:length(tetIDs)
                plot3([obj.Tracts.(tetIDs{t}){:,2}],...
                    [obj.Tracts.(tetIDs{t}){:,1}],...
                    [obj.Tracts.(tetIDs{t}){:,3}],...
                    'linewidth', 2, 'color', 'k', 'marker', '.', 'markersize', 10);
                hold on;
                text(obj.Tracts.(tetIDs{t}){1,2},...
                    obj.Tracts.(tetIDs{t}){1,1},...
                    obj.Tracts.(tetIDs{t}){1,3},...
                    tetIDs{t});
                if isfield(obj.TractLMs, tetIDs{t})
                    minSSlog = obj.TractLMs.(tetIDs{t}).SS==min(obj.TractLMs.(tetIDs{t}).SS);
                    tempReg = plot3(obj.TractLMs.(tetIDs{t}).ML(:,minSSlog),...
                        obj.TractLMs.(tetIDs{t}).AP(:,minSSlog),...
                        obj.TractLMs.(tetIDs{t}).DV(:,minSSlog));
                    tempReg.Color = obj.raterColors(find(minSSlog),:);
                end
            end
            grid on;
            set(gca, 'XLim', [0 3200], 'YLim', [0 3200], 'ZLim', [0 3200], 'ZDir', 'reverse');
            xlabel(gca,'ML (microns)');
            ylabel(gca,'AP (microns)');
            zlabel(gca,'DV (microns)');
        end
        %% Plot Reconstructed Tracts
        function PlotReconstruction(obj)
            if isempty(obj.TractsReconstructed)
                obj.ReconstructTracts;
            end
            tetIDs = fieldnames(obj.TractsReconstructed);
            figure;
            m2Tracts = [];
            accTracts = [];
            prlTracts = [];
            for t = 1:length(tetIDs)
                plot3([obj.Tracts.(tetIDs{t}){:,2}],...
                    [obj.Tracts.(tetIDs{t}){:,1}],...
                    [obj.Tracts.(tetIDs{t}){:,3}],...
                    'linewidth', 2, 'color', 'k');
                hold on;
                text(obj.Tracts.(tetIDs{t}){1,2},...
                    obj.Tracts.(tetIDs{t}){1,1},...
                    obj.Tracts.(tetIDs{t}){1,3},...
                    tetIDs{t});
                %  TractsReconstructed = [apTract, mlTract, dvTract]
                plot3(obj.TractsReconstructed.(tetIDs{t})(:,2),...
                    obj.TractsReconstructed.(tetIDs{t})(:,1),...
                    obj.TractsReconstructed.(tetIDs{t})(:,3), 'r');
                
                % Identify potential regional boundaries
                m2Tracts = [m2Tracts; cell2mat(obj.Tracts.(tetIDs{t})(contains(obj.Tracts.(tetIDs{t})(:,end), 'M2'),1:end-1))]; %#ok<AGROW> 
                accTracts = [accTracts; cell2mat(obj.Tracts.(tetIDs{t})(contains(obj.Tracts.(tetIDs{t})(:,end), 'ACC'),1:end-1))]; %#ok<AGROW> 
                prlTracts = [prlTracts; cell2mat(obj.Tracts.(tetIDs{t})(contains(obj.Tracts.(tetIDs{t})(:,end), 'PrL'),1:end-1))]; %#ok<AGROW> 
            end
            m2 = scatter3(m2Tracts(:,2), m2Tracts(:,1), m2Tracts(:,3), 100, [0.5,0.5,0.1], 'filled', 'o');
            acc = scatter3(accTracts(:,2), accTracts(:,1), accTracts(:,3), 100, [0.5,0.1,0.1], 'filled', 'o');
            prl = scatter3(prlTracts(:,2), prlTracts(:,1), prlTracts(:,3), 100, [0.1,0.5,0.5], 'filled', 'o');
            legend([m2, acc, prl], [{'M2'}, {'ACC'}, {'PrL'}]);
            grid on;
            set(gca, 'XLim', [-1000 3200], 'ZDir', 'reverse');
            xlabel(gca,'ML (microns)');
            ylabel(gca,'AP (microns)');
            zlabel(gca,'DV (microns)');


        end
    end
end