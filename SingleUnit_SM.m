classdef SingleUnit_SM < SeqMem
    % Single unit analyses from the Sequence Memory project
    properties % Analysis Variables
        binSize
        dsRate
        binType
        window
        alignment
    end
    properties % General Variables
        unitIDs
        tetIDs
        sessionSpikes
    end
    %% Creation Method
    methods
        %% Main Creation Method
        function obj = SingleUnit_SM(fileDir)
            if nargin == 0
                fileDir = uigetdir;
            end
            obj@SeqMem(fileDir);        
        end
    end
end