classdef ElectrodeTracts < handle
    %% Properties
    properties
        Subject
        SlideHeight % in microns
        SliceWidth  % in microns
        TetrodeIDs
        TractsByRater
        Tracts
        TractLMs
    end
    %% Creation Method
    methods
        function obj = ElectrodeTracts(varargin)
            if nargin == 1
                obj.Subject = varargin{1};
            elseif nargin == 3
                obj.Subject = varargin{1};
                obj.SlideHeight = varargin{2};
                obj.SliceWidth = varargin{3};
            end
        end
    end
    %% Import Methods
    methods
        %% Read Tract Annotations doc
        function ReadTractAnnotations(obj,fn,rater)
            % 'readcell' the excel file
            % Pull out the tetrode IDs from the file
            % Identify the measurement columns
            % Compile everything into structure
            % Store structure as 'tractsByRater'
        end
        %% Read Driving Records
        function ReadDriveRecord(obj,fn)

        end
    end
    %% Processing Methods
    methods
        %% Compile Tracts Across Raters
        function CompileTracts(obj)
        end
        %% Fit tracts to linear models
        function FitTractLM(obj)
            % use 'regress' to create the tractLMs 
        end
        %% 
    end
    %% Visualization Methods
    methods
        function PlotTracts(obj)
        end
    end