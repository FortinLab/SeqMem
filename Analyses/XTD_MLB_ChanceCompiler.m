function XTD_MLB_ChanceCompiler(fileDir,analysisWindow,numPerms)
    %% XTD_MLB_ChanceCompiler
    % m-file to compile the mlb cross-temporal data as well as chance data
    if nargin == 1 || isempty(analysisWindow)
        analysisWindow = [-1400 1400];
    end
    %% Create MLB object
    mlb = MLB_SM(fileDir);
    mlb.binSize = 200;
    mlb.dsRate = 50;
    mlb.windows = {[-3000 3000]};
    mlb.bayesType = 1;
    mlb.alignments = {'PokeIn'};
    mlb.SetLikes_ISC;
    fprintf('Processing Iterative Likelihoods via Leave-1-Out...')
    mlb.Process_IterativeLikelyL1O;
    mlb.Create_TrialInfoMasks;
    fprintf('complete\n');
    % Extract out HR, Decodability and TrialInfo structures in the trial history transMat organization
    [xtdHR_Real, xtdD_Real, trialInfo] = mlb.OrganizeDecodabilityTrialHistoryTransMat;
    % Extract out the cross-interval decodings
    [xInt_prePreInt_Real,...
        xInt_pstPstInt_Real,...
        xInt_prePstInt_Real,...
        xInt_pstPreInt_Real] = mlb.CalcXintDiff(xtdHR_Real,trialInfo);
    % Calculate the symmetry values using the decodability values
    [symTR_Real, symDEC_Real, symMN_Real] = mlb.CalcSymmetry_D(xtdD_Real);
    % Calculate persistence model fits
    [persTrialFit_Real, persIntervalFit_Real] = mlb.CalcModelPersistenceFit_XTD(xtdHR_Real,trialInfo,'latency');
    % Calculate Decoding Peaks w/in analysis window
    [decPeaksNdx_Real,...
        decPeaksVal_Real,...
        decPeaksWids_Real] = mlb.CalcDecodablityPeaks_XTD_Windowed(xtdHR_Real,analysisWindow);
    % Extract Exemplar XTD Model Decodings
    [modInt_PosDec_Real,...
        modPreDec_PosDec_Real,...
        modPstDec_PosDec_Real, offsets] = mlb.ExtractExemplarXTDmodels(xtdHR_Real,[],trialInfo,'windowed',analysisWindow);
    % Extract out Exemplar XTD Model Decoding Peaks
    modInt_PeakNdx_Real = cell(size(decPeaksNdx_Real));
    modInt_PeakVals_Real = cell(size(decPeaksNdx_Real));
    modInt_PeakWids_Real = cell(size(decPeaksNdx_Real));
    modPreDec_PeakNdx_Real = cell(size(decPeaksNdx_Real));
    modPreDec_PeakVals_Real = cell(size(decPeaksNdx_Real));
    modPreDec_PeakWids_Real = cell(size(decPeaksNdx_Real));
    modPstDec_PeakNdx_Real = cell(size(decPeaksNdx_Real));
    modPstDec_PeakVals_Real = cell(size(decPeaksNdx_Real));
    modPstDec_PeakWids_Real = cell(size(decPeaksNdx_Real));
    for p = 1:numel(decPeaksNdx_Real)
        if ~isempty(decPeaksNdx_Real{p})
            temp_trialIDs = [trialInfo{p}.TrialNum];
            modInt_PeakNdx_Real{p} = mlb.MaskVectExtract(decPeaksNdx_Real{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
            modInt_PeakVals_Real{p} = mlb.MaskVectExtract(decPeaksVal_Real{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
            modInt_PeakWids_Real{p} = mlb.MaskVectExtract(decPeaksWids_Real{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
            modPreDec_PeakNdx_Real{p} = mlb.MaskVectExtract(decPeaksNdx_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
            modPreDec_PeakVals_Real{p} = mlb.MaskVectExtract(decPeaksVal_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
            modPreDec_PeakWids_Real{p} = mlb.MaskVectExtract(decPeaksWids_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
            modPstDec_PeakNdx_Real{p} = mlb.MaskVectExtract(decPeaksNdx_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
            modPstDec_PeakVals_Real{p} = mlb.MaskVectExtract(decPeaksVal_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
            modPstDec_PeakWids_Real{p} = mlb.MaskVectExtract(decPeaksWids_Real{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
        end
    end
    
    %% Run Chance Permutations
    for perm = 1:numPerms
        fprintf('Starting perm %i @%s..... ',perm, datetime('now', 'Format', 'HH:mm:ss'));
        %% Run FULL (trial & time ID) permutation
        mlb.Process_IterativeLikelyL1O('Full');
        [cr_xtdHR, cr_xtdD, ~] = mlb.OrganizeDecodabilityTrialHistoryTransMat;
        if perm == 1
            xtdHR_ChanceFull = cr_xtdHR;
            xtdHR_ChanceFullCount = cellfun(@(a){~isnan(a)},cr_xtdHR);
            xtdD_ChanceFull = cr_xtdD;
            xtdD_ChanceFullCount = cellfun(@(a){~isnan(a)},cr_xtdD);
        else
            [xtdHR_ChanceFull, xtdHR_ChanceFullCount] = UpdateChanceVar(xtdHR_ChanceFull, xtdHR_ChanceFullCount, cr_xtdHR);
            [xtdD_ChanceFull, xtdD_ChanceFullCount] = UpdateChanceVar(xtdD_ChanceFull, xtdD_ChanceFullCount, cr_xtdD);
            if perm == numPerms
                [xtdHR_ChanceFull] = AverageChanceVar(xtdHR_ChanceFull, xtdHR_ChanceFullCount);
                [xtdD_ChanceFull] = AverageChanceVar(xtdD_ChanceFull, xtdD_ChanceFullCount);
                clear xtdHR_ChanceFullCount xtdD_ChanceFullCount
            end
        end
        % Extract out the cross-interval decodings
        [cr_xInt_prePreInt,...
            cr_xInt_pstPstInt,...
            cr_xInt_prePstInt,...
            cr_xInt_pstPreInt] = mlb.CalcXintDiff(cr_xtdHR,trialInfo);
        if perm == 1
            xInt_prePreInt_ChanceFull = cr_xInt_prePreInt;
            xInt_prePreInt_ChanceFullCount = cellfun(@(a){~isnan(a)},xInt_prePreInt_ChanceFull);
            xInt_pstPstInt_ChanceFull = cr_xInt_pstPstInt;
            xInt_pstPstInt_ChanceFullCount = cellfun(@(a){~isnan(a)},xInt_pstPstInt_ChanceFull);
            xInt_prePstInt_ChanceFull = cr_xInt_prePstInt;
            xInt_prePstInt_ChanceFullCount = cellfun(@(a){~isnan(a)},xInt_prePstInt_ChanceFull);
            xInt_pstPreInt_ChanceFull = cr_xInt_pstPreInt;
            xInt_pstPreInt_ChanceFullCount = cellfun(@(a){~isnan(a)},xInt_pstPreInt_ChanceFull);
        else
            [xInt_prePreInt_ChanceFull, xInt_prePreInt_ChanceFullCount] = UpdateChanceVar(xInt_prePreInt_ChanceFull, xInt_prePreInt_ChanceFullCount, cr_xInt_prePreInt);
            [xInt_pstPstInt_ChanceFull, xInt_pstPstInt_ChanceFullCount] = UpdateChanceVar(xInt_pstPstInt_ChanceFull, xInt_pstPstInt_ChanceFullCount, cr_xInt_pstPstInt);
            [xInt_prePstInt_ChanceFull, xInt_prePstInt_ChanceFullCount] = UpdateChanceVar(xInt_prePstInt_ChanceFull, xInt_prePstInt_ChanceFullCount, cr_xInt_prePstInt);
            [xInt_pstPreInt_ChanceFull, xInt_pstPreInt_ChanceFullCount] = UpdateChanceVar(xInt_pstPreInt_ChanceFull, xInt_pstPreInt_ChanceFullCount, cr_xInt_pstPreInt);
            if perm == numPerms
                [xInt_prePreInt_ChanceFull] = AverageChanceVar(xInt_prePreInt_ChanceFull, xInt_prePreInt_ChanceFullCount);
                [xInt_pstPstInt_ChanceFull] = AverageChanceVar(xInt_pstPstInt_ChanceFull, xInt_pstPstInt_ChanceFullCount);
                [xInt_prePstInt_ChanceFull] = AverageChanceVar(xInt_prePstInt_ChanceFull, xInt_prePstInt_ChanceFullCount);
                [xInt_pstPreInt_ChanceFull] = AverageChanceVar(xInt_pstPreInt_ChanceFull, xInt_pstPreInt_ChanceFullCount);
                clear xInt_prePreInt_ChanceFullCount xInt_pstPstInt_ChanceFullCount xInt_prePstInt_ChanceFullCount xInt_pstPreInt_ChanceFullCount
            end
        end
    
        % Calculate the symmetry values using the decodability values
        [cr_symTR, cr_symDEC, cr_symMN] = mlb.CalcSymmetry_D(cr_xtdD);
        if perm == 1
            symTR_ChanceFull = cr_symTR;
            symTR_ChanceFullCount = cellfun(@(a){~isnan(a)},symTR_ChanceFull);
            symDEC_ChanceFull = cr_symDEC;
            symDEC_ChanceFullCount = cellfun(@(a){~isnan(a)},symDEC_ChanceFull);
            symMN_ChanceFull = cr_symMN;
            symMN_ChanceFullCount = cellfun(@(a){~isnan(a)},symMN_ChanceFull);
        else
            [symTR_ChanceFull, symTR_ChanceFullCount] = UpdateChanceVar(symTR_ChanceFull, symTR_ChanceFullCount, cr_symTR);
            [symDEC_ChanceFull, symDEC_ChanceFullCount] = UpdateChanceVar(symDEC_ChanceFull, symDEC_ChanceFullCount, cr_symDEC);
            [symMN_ChanceFull, symMN_ChanceFullCount] = UpdateChanceVar(symMN_ChanceFull, symMN_ChanceFullCount, cr_symMN);
            if perm == numPerms
                [symTR_ChanceFull] = AverageChanceVar(symTR_ChanceFull, symTR_ChanceFullCount);
                [symDEC_ChanceFull] = AverageChanceVar(symDEC_ChanceFull, symDEC_ChanceFullCount);
                [symMN_ChanceFull] = AverageChanceVar(symMN_ChanceFull, symMN_ChanceFullCount);
                clear symTR_ChanceFullCount symDEC_ChanceFullCount symMN_ChanceFullCount
            end
        end

        % Calculate persistence model fits
        [cr_persTrialFit, cr_persIntervalFit] = mlb.CalcModelPersistenceFit_XTD(cr_xtdHR,trialInfo,'latency');
        if perm == 1
            persTrialFit_ChanceFull = cr_persTrialFit;
            persTrialFit_ChanceFullCount = cellfun(@(a){~isnan(a)}, persTrialFit_ChanceFull);
            persIntervalFit_ChanceFull = cr_persIntervalFit;
            persIntervalFit_ChanceFullCount = cellfun(@(a){~isnan(a)}, persIntervalFit_ChanceFull);
        else
            [persTrialFit_ChanceFull, persTrialFit_ChanceFullCount] = UpdateChanceVar(persTrialFit_ChanceFull, persTrialFit_ChanceFullCount,cr_persTrialFit);
            [persIntervalFit_ChanceFull, persIntervalFit_ChanceFullCount] = UpdateChanceVar(persIntervalFit_ChanceFull, persIntervalFit_ChanceFullCount,cr_persIntervalFit);
            if perm == numPerms
                [persTrialFit_ChanceFull] = AverageChanceVar(persTrialFit_ChanceFull, persTrialFit_ChanceFullCount);
                [persIntervalFit_ChanceFull] = AverageChanceVar(persIntervalFit_ChanceFull, persIntervalFit_ChanceFullCount);
                clear persTrialFit_ChanceFullCount persIntervalFit_ChanceFullCount
            end
        end

        % Calculate Decoding Peaks w/in analysis window
        [cr_decPeaksNdx,...
            cr_decPeaksVal,...
            cr_decPeaksWid] = mlb.CalcDecodablityPeaks_XTD_Windowed(cr_xtdHR,analysisWindow);
        if perm == 1
            decPeaksNdx_ChanceFull = cr_decPeaksNdx;
            decPeaksNdx_ChanceFullCount = cellfun(@(a){~isnan(a)}, decPeaksNdx_ChanceFull);
            decPeaksVal_ChanceFull = cr_decPeaksVal;
            decPeaksVal_ChanceFullCount = cellfun(@(a){~isnan(a)}, decPeaksVal_ChanceFull);
            decPeaksWid_ChanceFull = cr_decPeaksWid;
            decPeaksWid_ChanceFullCount = cellfun(@(a){~isnan(a)}, decPeaksWid_ChanceFull);
        else
            [decPeaksNdx_ChanceFull, decPeaksNdx_ChanceFullCount] = UpdateChanceVar(decPeaksNdx_ChanceFull, decPeaksNdx_ChanceFullCount, cr_decPeaksNdx);
            [decPeaksVal_ChanceFull, decPeaksVal_ChanceFullCount] = UpdateChanceVar(decPeaksVal_ChanceFull, decPeaksVal_ChanceFullCount, cr_decPeaksVal);
            [decPeaksWid_ChanceFull, decPeaksWid_ChanceFullCount] = UpdateChanceVar(decPeaksWid_ChanceFull, decPeaksWid_ChanceFullCount, cr_decPeaksWid);
            if perm == numPerms
                [decPeaksNdx_ChanceFull] = AverageChanceVar(decPeaksNdx_ChanceFull, decPeaksNdx_ChanceFullCount);
                [decPeaksVal_ChanceFull] = AverageChanceVar(decPeaksVal_ChanceFull, decPeaksVal_ChanceFullCount);
                [decPeaksWid_ChanceFull] = AverageChanceVar(decPeaksWid_ChanceFull, decPeaksWid_ChanceFullCount);
                clear decPeaksNdx_ChanceFullCount decPeaksWid_ChanceFullCount
            end
        end
        % Extract Exemplar XTD Model Decodings
        [cr_modInt_PosDec,...
            cr_modPreDec_PosDec,...
            cr_modPstDec_PosDec] = mlb.ExtractExemplarXTDmodels(cr_xtdHR,[],trialInfo,'windowed',analysisWindow);
        if perm == 1
            modInt_PosDec_ChanceFull = cr_modInt_PosDec;
            modInt_PosDec_ChanceFullCount = cellfun(@(a){~isnan(a)},modInt_PosDec_ChanceFull);
            modPreDec_PosDec_ChanceFull = cr_modPreDec_PosDec;
            modPreDec_PosDec_ChanceFullCount = cellfun(@(a){~isnan(a)},modPreDec_PosDec_ChanceFull);
            modPstDec_PosDec_ChanceFull = cr_modPstDec_PosDec;
            modPstDec_PosDec_ChanceFullCount = cellfun(@(a){~isnan(a)},modPstDec_PosDec_ChanceFull);
        else
            [modInt_PosDec_ChanceFull, modInt_PosDec_ChanceFullCount] = UpdateChanceVar(modInt_PosDec_ChanceFull, modInt_PosDec_ChanceFullCount, cr_modInt_PosDec);
            [modPreDec_PosDec_ChanceFull, modPreDec_PosDec_ChanceFullCount] = UpdateChanceVar(modPreDec_PosDec_ChanceFull, modPreDec_PosDec_ChanceFullCount, cr_modPreDec_PosDec);
            [modPstDec_PosDec_ChanceFull, modPstDec_PosDec_ChanceFullCount] = UpdateChanceVar(modPstDec_PosDec_ChanceFull, modPstDec_PosDec_ChanceFullCount, cr_modPstDec_PosDec);
            if perm == numPerms
                [modInt_PosDec_ChanceFull] = AverageChanceVar(modInt_PosDec_ChanceFull, modInt_PosDec_ChanceFullCount);
                [modPreDec_PosDec_ChanceFull] = AverageChanceVar(modPreDec_PosDec_ChanceFull, modPreDec_PosDec_ChanceFullCount);
                [modPstDec_PosDec_ChanceFull] = AverageChanceVar(modPstDec_PosDec_ChanceFull, modPstDec_PosDec_ChanceFullCount);
                clear modInt_PosDec_ChanceFullCount modPreDec_PosDec_ChanceFullCount modPstDec_PosDec_ChanceFullCount
            end
        end
        % Extract out Exemplar XTD Model Decoding Peaks
        cr_modInt_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modInt_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modInt_PeakWids = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakWids = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakWids = cell(size(cr_decPeaksNdx));
        for p = 1:numel(cr_decPeaksNdx)
            if ~isempty(cr_decPeaksNdx{p})
                temp_trialIDs = [trialInfo{p}.TrialNum];
                cr_modInt_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modInt_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modInt_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modPreDec_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPreDec_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPreDec_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPstDec_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
                cr_modPstDec_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
                cr_modPstDec_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
            end
        end
        if perm == 1
            modInt_PeakNdx_ChanceFull = cr_modInt_PeakNdx;
            modInt_PeakNdx_ChanceFullCount = cellfun(@(a){~isnan(a)},modInt_PeakNdx_ChanceFull);
            modInt_PeakVals_ChanceFull = cr_modInt_PeakVals;
            modInt_PeakVals_ChanceFullCount = cellfun(@(a){~isnan(a)},modInt_PeakVals_ChanceFull);
            modInt_PeakWids_ChanceFull = cr_modInt_PeakWids;
            modInt_PeakWids_ChanceFullCount = cellfun(@(a){~isnan(a)},modInt_PeakWids_ChanceFull);
            modPreDec_PeakNdx_ChanceFull = cr_modPreDec_PeakNdx;
            modPreDec_PeakNdx_ChanceFullCount = cellfun(@(a){~isnan(a)},modPreDec_PeakNdx_ChanceFull);
            modPreDec_PeakVals_ChanceFull = cr_modPreDec_PeakVals;
            modPreDec_PeakVals_ChanceFullCount = cellfun(@(a){~isnan(a)},modPreDec_PeakVals_ChanceFull);
            modPreDec_PeakWids_ChanceFull = cr_modPreDec_PeakWids;
            modPreDec_PeakWids_ChanceFullCount = cellfun(@(a){~isnan(a)},modPreDec_PeakWids_ChanceFull);
            modPstDec_PeakNdx_ChanceFull = cr_modPstDec_PeakNdx;
            modPstDec_PeakNdx_ChanceFullCount = cellfun(@(a){~isnan(a)},modPstDec_PeakNdx_ChanceFull);
            modPstDec_PeakVals_ChanceFull = cr_modPstDec_PeakVals;
            modPstDec_PeakVals_ChanceFullCount = cellfun(@(a){~isnan(a)},modPstDec_PeakVals_ChanceFull);
            modPstDec_PeakWids_ChanceFull = cr_modPstDec_PeakWids;
            modPstDec_PeakWids_ChanceFullCount = cellfun(@(a){~isnan(a)},modPstDec_PeakWids_ChanceFull);
        else
            [modInt_PeakNdx_ChanceFull, modInt_PeakNdx_ChanceFullCount] = UpdateChanceVar(modInt_PeakNdx_ChanceFull, modInt_PeakNdx_ChanceFullCount, cr_modInt_PeakNdx);
            [modInt_PeakVals_ChanceFull, modInt_PeakVals_ChanceFullCount] = UpdateChanceVar(modInt_PeakVals_ChanceFull, modInt_PeakVals_ChanceFullCount, cr_modInt_PeakVals);
            [modInt_PeakWids_ChanceFull, modInt_PeakWids_ChanceFullCount] = UpdateChanceVar(modInt_PeakWids_ChanceFull, modInt_PeakWids_ChanceFullCount, cr_modInt_PeakWids);
            [modPreDec_PeakNdx_ChanceFull, modPreDec_PeakNdx_ChanceFullCount] = UpdateChanceVar(modPreDec_PeakNdx_ChanceFull, modPreDec_PeakNdx_ChanceFullCount, cr_modPreDec_PeakNdx);
            [modPreDec_PeakVals_ChanceFull, modPreDec_PeakVals_ChanceFullCount] = UpdateChanceVar(modPreDec_PeakVals_ChanceFull, modPreDec_PeakVals_ChanceFullCount, cr_modPreDec_PeakVals);
            [modPreDec_PeakWids_ChanceFull, modPreDec_PeakWids_ChanceFullCount] = UpdateChanceVar(modPreDec_PeakWids_ChanceFull, modPreDec_PeakWids_ChanceFullCount, cr_modPreDec_PeakWids);
            [modPstDec_PeakNdx_ChanceFull, modPstDec_PeakNdx_ChanceFullCount] = UpdateChanceVar(modPstDec_PeakNdx_ChanceFull, modPstDec_PeakNdx_ChanceFullCount, cr_modPstDec_PeakNdx);
            [modPstDec_PeakVals_ChanceFull, modPstDec_PeakVals_ChanceFullCount] = UpdateChanceVar(modPstDec_PeakVals_ChanceFull, modPstDec_PeakVals_ChanceFullCount, cr_modPstDec_PeakVals);
            [modPstDec_PeakWids_ChanceFull, modPstDec_PeakWids_ChanceFullCount] = UpdateChanceVar(modPstDec_PeakWids_ChanceFull, modPstDec_PeakWids_ChanceFullCount, cr_modPstDec_PeakWids);
            if perm == numPerms
                [modInt_PeakNdx_ChanceFull] = AverageChanceVar(modInt_PeakNdx_ChanceFull, modInt_PeakNdx_ChanceFullCount);
                [modInt_PeakVals_ChanceFull] = AverageChanceVar(modInt_PeakVals_ChanceFull, modInt_PeakVals_ChanceFullCount);
                [modInt_PeakWids_ChanceFull] = AverageChanceVar(modInt_PeakWids_ChanceFull, modInt_PeakWids_ChanceFullCount);
                [modPreDec_PeakNdx_ChanceFull] = AverageChanceVar(modPreDec_PeakNdx_ChanceFull, modPreDec_PeakNdx_ChanceFullCount);
                [modPreDec_PeakVals_ChanceFull] = AverageChanceVar(modPreDec_PeakVals_ChanceFull, modPreDec_PeakVals_ChanceFullCount);
                [modPreDec_PeakWids_ChanceFull] = AverageChanceVar(modPreDec_PeakWids_ChanceFull, modPreDec_PeakWids_ChanceFullCount);
                [modPstDec_PeakNdx_ChanceFull] = AverageChanceVar(modPstDec_PeakNdx_ChanceFull, modPstDec_PeakNdx_ChanceFullCount);
                [modPstDec_PeakVals_ChanceFull] = AverageChanceVar(modPstDec_PeakVals_ChanceFull, modPstDec_PeakVals_ChanceFullCount);
                [modPstDec_PeakWids_ChanceFull] = AverageChanceVar(modPstDec_PeakWids_ChanceFull, modPstDec_PeakWids_ChanceFullCount);
                clear modInt_PeakNdx_ChanceFullCount modInt_PeakVals_ChanceFullCount modInt_PeakWids_ChanceFullCount...
                    modPreDec_PeakNdx_ChanceFullCount modPreDec_PeakVals_ChanceFullCount modPreDec_PeakWids_ChanceFullCount...
                    modPstDec_PeakNdx_ChanceFullCount modPstDec_PeakVals_ChanceFullCount modPstDec_PeakWids_ChanceFullCount
            end
        end

        %% Run TRIAL ONLY permutation
        mlb.Process_IterativeLikelyL1O('Trial');
        [cr_xtdHR, cr_xtdD, ~] = mlb.OrganizeDecodabilityTrialHistoryTransMat;
        if perm == 1
            xtdHR_ChancePos = cr_xtdHR;
            xtdHR_ChancePosCount = cellfun(@(a){~isnan(a)},cr_xtdHR);
            xtdD_ChancePos = cr_xtdD;
            xtdD_ChancePosCount = cellfun(@(a){~isnan(a)},cr_xtdD);
        else
            [xtdHR_ChancePos, xtdHR_ChancePosCount] = UpdateChanceVar(xtdHR_ChancePos, xtdHR_ChancePosCount, cr_xtdHR);
            [xtdD_ChancePos, xtdD_ChancePosCount] = UpdateChanceVar(xtdD_ChancePos, xtdD_ChancePosCount, cr_xtdD);
            if perm == numPerms
                [xtdHR_ChancePos] = AverageChanceVar(xtdHR_ChancePos, xtdHR_ChancePosCount);
                [xtdD_ChancePos] = AverageChanceVar(xtdD_ChancePos, xtdD_ChancePosCount);
                clear xtdHR_ChancePosCount xtdD_ChancePosCount
            end
        end
        % Extract out the cross-interval decodings
        [cr_xInt_prePreInt,...
            cr_xInt_pstPstInt,...
            cr_xInt_prePstInt,...
            cr_xInt_pstPreInt] = mlb.CalcXintDiff(cr_xtdHR,trialInfo);
        if perm == 1
            xInt_prePreInt_ChancePos = cr_xInt_prePreInt;
            xInt_prePreInt_ChancePosCount = cellfun(@(a){~isnan(a)},xInt_prePreInt_ChancePos);
            xInt_pstPstInt_ChancePos = cr_xInt_pstPstInt;
            xInt_pstPstInt_ChancePosCount = cellfun(@(a){~isnan(a)},xInt_pstPstInt_ChancePos);
            xInt_prePstInt_ChancePos = cr_xInt_prePstInt;
            xInt_prePstInt_ChancePosCount = cellfun(@(a){~isnan(a)},xInt_prePstInt_ChancePos);
            xInt_pstPreInt_ChancePos = cr_xInt_pstPreInt;
            xInt_pstPreInt_ChancePosCount = cellfun(@(a){~isnan(a)},xInt_pstPreInt_ChancePos);
        else
            [xInt_prePreInt_ChancePos, xInt_prePreInt_ChancePosCount] = UpdateChanceVar(xInt_prePreInt_ChancePos, xInt_prePreInt_ChancePosCount, cr_xInt_prePreInt);
            [xInt_pstPstInt_ChancePos, xInt_pstPstInt_ChancePosCount] = UpdateChanceVar(xInt_pstPstInt_ChancePos, xInt_pstPstInt_ChancePosCount, cr_xInt_pstPstInt);
            [xInt_prePstInt_ChancePos, xInt_prePstInt_ChancePosCount] = UpdateChanceVar(xInt_prePstInt_ChancePos, xInt_prePstInt_ChancePosCount, cr_xInt_prePstInt);
            [xInt_pstPreInt_ChancePos, xInt_pstPreInt_ChancePosCount] = UpdateChanceVar(xInt_pstPreInt_ChancePos, xInt_pstPreInt_ChancePosCount, cr_xInt_pstPreInt);
            if perm == numPerms
                [xInt_prePreInt_ChancePos] = AverageChanceVar(xInt_prePreInt_ChancePos, xInt_prePreInt_ChancePosCount);
                [xInt_pstPstInt_ChancePos] = AverageChanceVar(xInt_pstPstInt_ChancePos, xInt_pstPstInt_ChancePosCount);
                [xInt_prePstInt_ChancePos] = AverageChanceVar(xInt_prePstInt_ChancePos, xInt_prePstInt_ChancePosCount);
                [xInt_pstPreInt_ChancePos] = AverageChanceVar(xInt_pstPreInt_ChancePos, xInt_pstPreInt_ChancePosCount);
                clear xInt_prePreInt_ChancePosCount xInt_pstPstInt_ChancePosCount xInt_prePstInt_ChancePosCount xInt_pstPreInt_ChancePosCount
            end
        end
    
        % Calculate the symmetry values using the decodability values
        [cr_symTR, cr_symDEC, cr_symMN] = mlb.CalcSymmetry_D(cr_xtdD);
        if perm == 1
            symTR_ChancePos = cr_symTR;
            symTR_ChancePosCount = cellfun(@(a){~isnan(a)},symTR_ChancePos);
            symDEC_ChancePos = cr_symDEC;
            symDEC_ChancePosCount = cellfun(@(a){~isnan(a)},symDEC_ChancePos);
            symMN_ChancePos = cr_symMN;
            symMN_ChancePosCount = cellfun(@(a){~isnan(a)},symMN_ChancePos);
        else
            [symTR_ChancePos, symTR_ChancePosCount] = UpdateChanceVar(symTR_ChancePos, symTR_ChancePosCount, cr_symTR);
            [symDEC_ChancePos, symDEC_ChancePosCount] = UpdateChanceVar(symDEC_ChancePos, symDEC_ChancePosCount, cr_symDEC);
            [symMN_ChancePos, symMN_ChancePosCount] = UpdateChanceVar(symMN_ChancePos, symMN_ChancePosCount, cr_symMN);
            if perm == numPerms
                [symTR_ChancePos] = AverageChanceVar(symTR_ChancePos, symTR_ChancePosCount);
                [symDEC_ChancePos] = AverageChanceVar(symDEC_ChancePos, symDEC_ChancePosCount);
                [symMN_ChancePos] = AverageChanceVar(symMN_ChancePos, symMN_ChancePosCount);
                clear symTR_ChancePosCount symDEC_ChancePosCount symMN_ChancePosCount
            end
        end

        % Calculate persistence model fits
        [cr_persTrialFit, cr_persIntervalFit] = mlb.CalcModelPersistenceFit_XTD(cr_xtdHR,trialInfo,'latency');
        if perm == 1
            persTrialFit_ChancePos = cr_persTrialFit;
            persTrialFit_ChancePosCount = cellfun(@(a){~isnan(a)}, persTrialFit_ChancePos);
            persIntervalFit_ChancePos = cr_persIntervalFit;
            persIntervalFit_ChancePosCount = cellfun(@(a){~isnan(a)}, persIntervalFit_ChancePos);
        else
            [persTrialFit_ChancePos, persTrialFit_ChancePosCount] = UpdateChanceVar(persTrialFit_ChancePos, persTrialFit_ChancePosCount,cr_persTrialFit);
            [persIntervalFit_ChancePos, persIntervalFit_ChancePosCount] = UpdateChanceVar(persIntervalFit_ChancePos, persIntervalFit_ChancePosCount,cr_persIntervalFit);
            if perm == numPerms
                [persTrialFit_ChancePos] = AverageChanceVar(persTrialFit_ChancePos, persTrialFit_ChancePosCount);
                [persIntervalFit_ChancePos] = AverageChanceVar(persIntervalFit_ChancePos, persIntervalFit_ChancePosCount);
                clear persTrialFit_ChancePosCount persIntervalFit_ChancePosCount
            end
        end

        % Calculate Decoding Peaks w/in analysis window
        [cr_decPeaksNdx,...
            cr_decPeaksVal,...
            cr_decPeaksWid] = mlb.CalcDecodablityPeaks_XTD_Windowed(cr_xtdHR,analysisWindow);
        if perm == 1
            decPeaksNdx_ChancePos = cr_decPeaksNdx;
            decPeaksNdx_ChancePosCount = cellfun(@(a){~isnan(a)}, decPeaksNdx_ChancePos);
            decPeaksVal_ChancePos = cr_decPeaksVal;
            decPeaksVal_ChancePosCount = cellfun(@(a){~isnan(a)}, decPeaksVal_ChancePos);
            decPeaksWid_ChancePos = cr_decPeaksWid;
            decPeaksWid_ChancePosCount = cellfun(@(a){~isnan(a)}, decPeaksWid_ChancePos);
        else
            [decPeaksNdx_ChancePos, decPeaksNdx_ChancePosCount] = UpdateChanceVar(decPeaksNdx_ChancePos, decPeaksNdx_ChancePosCount, cr_decPeaksNdx);
            [decPeaksVal_ChancePos, decPeaksVal_ChancePosCount] = UpdateChanceVar(decPeaksVal_ChancePos, decPeaksVal_ChancePosCount, cr_decPeaksVal);
            [decPeaksWid_ChancePos, decPeaksWid_ChancePosCount] = UpdateChanceVar(decPeaksWid_ChancePos, decPeaksWid_ChancePosCount, cr_decPeaksWid);
            if perm == numPerms
                [decPeaksNdx_ChancePos] = AverageChanceVar(decPeaksNdx_ChancePos, decPeaksNdx_ChancePosCount);
                [decPeaksVal_ChancePos] = AverageChanceVar(decPeaksVal_ChancePos, decPeaksVal_ChancePosCount);
                [decPeaksWid_ChancePos] = AverageChanceVar(decPeaksWid_ChancePos, decPeaksWid_ChancePosCount);
                clear decPeaksNdx_ChancePosCount decPeaksWid_ChancePosCount
            end
        end
        % Extract Exemplar XTD Model Decodings
        [cr_modInt_PosDec,...
            cr_modPreDec_PosDec,...
            cr_modPstDec_PosDec] = mlb.ExtractExemplarXTDmodels(cr_xtdHR,[],trialInfo,'windowed',analysisWindow);
        if perm == 1
            modInt_PosDec_ChancePos = cr_modInt_PosDec;
            modInt_PosDec_ChancePosCount = cellfun(@(a){~isnan(a)},modInt_PosDec_ChancePos);
            modPreDec_PosDec_ChancePos = cr_modPreDec_PosDec;
            modPreDec_PosDec_ChancePosCount = cellfun(@(a){~isnan(a)},modPreDec_PosDec_ChancePos);
            modPstDec_PosDec_ChancePos = cr_modPstDec_PosDec;
            modPstDec_PosDec_ChancePosCount = cellfun(@(a){~isnan(a)},modPstDec_PosDec_ChancePos);
        else
            [modInt_PosDec_ChancePos, modInt_PosDec_ChancePosCount] = UpdateChanceVar(modInt_PosDec_ChancePos, modInt_PosDec_ChancePosCount, cr_modInt_PosDec);
            [modPreDec_PosDec_ChancePos, modPreDec_PosDec_ChancePosCount] = UpdateChanceVar(modPreDec_PosDec_ChancePos, modPreDec_PosDec_ChancePosCount, cr_modPreDec_PosDec);
            [modPstDec_PosDec_ChancePos, modPstDec_PosDec_ChancePosCount] = UpdateChanceVar(modPstDec_PosDec_ChancePos, modPstDec_PosDec_ChancePosCount, cr_modPstDec_PosDec);
            if perm == numPerms
                [modInt_PosDec_ChancePos] = AverageChanceVar(modInt_PosDec_ChancePos, modInt_PosDec_ChancePosCount);
                [modPreDec_PosDec_ChancePos] = AverageChanceVar(modPreDec_PosDec_ChancePos, modPreDec_PosDec_ChancePosCount);
                [modPstDec_PosDec_ChancePos] = AverageChanceVar(modPstDec_PosDec_ChancePos, modPstDec_PosDec_ChancePosCount);
                clear modInt_PosDec_ChancePosCount modPreDec_PosDec_ChancePosCount modPstDec_PosDec_ChancePosCount
            end
        end
        % Extract out Exemplar XTD Model Decoding Peaks
        cr_modInt_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modInt_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modInt_PeakWids = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modPreDec_PeakWids = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakNdx = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakVals = cell(size(cr_decPeaksNdx));
        cr_modPstDec_PeakWids = cell(size(cr_decPeaksNdx));
        for p = 1:numel(cr_decPeaksNdx)
            if ~isempty(cr_decPeaksNdx{p})
                temp_trialIDs = [trialInfo{p}.TrialNum];
                cr_modInt_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modInt_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modInt_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},mlb.mask_IntMid(:,temp_trialIDs),2,3);
                cr_modPreDec_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPreDec_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPreDec_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(1),offsets(1))}],2,3);
                cr_modPstDec_PeakNdx{p} = mlb.MaskVectExtract(cr_decPeaksNdx{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
                cr_modPstDec_PeakVals{p} = mlb.MaskVectExtract(cr_decPeaksVal{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
                cr_modPstDec_PeakWids{p} = mlb.MaskVectExtract(cr_modPstDec_PosDec{p},[{mlb.mask_Trial(:,temp_trialIDs)}, {sprintf('first+%i+%i',offsets(2),offsets(2))}],2,3);
            end
        end
        if perm == 1
            modInt_PeakNdx_ChancePos = cr_modInt_PeakNdx;
            modInt_PeakNdx_ChancePosCount = cellfun(@(a){~isnan(a)},modInt_PeakNdx_ChancePos);
            modInt_PeakVals_ChancePos = cr_modInt_PeakVals;
            modInt_PeakVals_ChancePosCount = cellfun(@(a){~isnan(a)},modInt_PeakVals_ChancePos);
            modInt_PeakWids_ChancePos = cr_modInt_PeakWids;
            modInt_PeakWids_ChancePosCount = cellfun(@(a){~isnan(a)},modInt_PeakWids_ChancePos);
            modPreDec_PeakNdx_ChancePos = cr_modPreDec_PeakNdx;
            modPreDec_PeakNdx_ChancePosCount = cellfun(@(a){~isnan(a)},modPreDec_PeakNdx_ChancePos);
            modPreDec_PeakVals_ChancePos = cr_modPreDec_PeakVals;
            modPreDec_PeakVals_ChancePosCount = cellfun(@(a){~isnan(a)},modPreDec_PeakVals_ChancePos);
            modPreDec_PeakWids_ChancePos = cr_modPreDec_PeakWids;
            modPreDec_PeakWids_ChancePosCount = cellfun(@(a){~isnan(a)},modPreDec_PeakWids_ChancePos);
            modPstDec_PeakNdx_ChancePos = cr_modPstDec_PeakNdx;
            modPstDec_PeakNdx_ChancePosCount = cellfun(@(a){~isnan(a)},modPstDec_PeakNdx_ChancePos);
            modPstDec_PeakVals_ChancePos = cr_modPstDec_PeakVals;
            modPstDec_PeakVals_ChancePosCount = cellfun(@(a){~isnan(a)},modPstDec_PeakVals_ChancePos);
            modPstDec_PeakWids_ChancePos = cr_modPstDec_PeakWids;
            modPstDec_PeakWids_ChancePosCount = cellfun(@(a){~isnan(a)},modPstDec_PeakWids_ChancePos);
        else
            [modInt_PeakNdx_ChancePos, modInt_PeakNdx_ChancePosCount] = UpdateChanceVar(modInt_PeakNdx_ChancePos, modInt_PeakNdx_ChancePosCount, cr_modInt_PeakNdx);
            [modInt_PeakVals_ChancePos, modInt_PeakVals_ChancePosCount] = UpdateChanceVar(modInt_PeakVals_ChancePos, modInt_PeakVals_ChancePosCount, cr_modInt_PeakVals);
            [modInt_PeakWids_ChancePos, modInt_PeakWids_ChancePosCount] = UpdateChanceVar(modInt_PeakWids_ChancePos, modInt_PeakWids_ChancePosCount, cr_modInt_PeakWids);
            [modPreDec_PeakNdx_ChancePos, modPreDec_PeakNdx_ChancePosCount] = UpdateChanceVar(modPreDec_PeakNdx_ChancePos, modPreDec_PeakNdx_ChancePosCount, cr_modPreDec_PeakNdx);
            [modPreDec_PeakVals_ChancePos, modPreDec_PeakVals_ChancePosCount] = UpdateChanceVar(modPreDec_PeakVals_ChancePos, modPreDec_PeakVals_ChancePosCount, cr_modPreDec_PeakVals);
            [modPreDec_PeakWids_ChancePos, modPreDec_PeakWids_ChancePosCount] = UpdateChanceVar(modPreDec_PeakWids_ChancePos, modPreDec_PeakWids_ChancePosCount, cr_modPreDec_PeakWids);
            [modPstDec_PeakNdx_ChancePos, modPstDec_PeakNdx_ChancePosCount] = UpdateChanceVar(modPstDec_PeakNdx_ChancePos, modPstDec_PeakNdx_ChancePosCount, cr_modPstDec_PeakNdx);
            [modPstDec_PeakVals_ChancePos, modPstDec_PeakVals_ChancePosCount] = UpdateChanceVar(modPstDec_PeakVals_ChancePos, modPstDec_PeakVals_ChancePosCount, cr_modPstDec_PeakVals);
            [modPstDec_PeakWids_ChancePos, modPstDec_PeakWids_ChancePosCount] = UpdateChanceVar(modPstDec_PeakWids_ChancePos, modPstDec_PeakWids_ChancePosCount, cr_modPstDec_PeakWids);
            if perm == numPerms
                [modInt_PeakNdx_ChancePos] = AverageChanceVar(modInt_PeakNdx_ChancePos, modInt_PeakNdx_ChancePosCount);
                [modInt_PeakVals_ChancePos] = AverageChanceVar(modInt_PeakVals_ChancePos, modInt_PeakVals_ChancePosCount);
                [modInt_PeakWids_ChancePos] = AverageChanceVar(modInt_PeakWids_ChancePos, modInt_PeakWids_ChancePosCount);
                [modPreDec_PeakNdx_ChancePos] = AverageChanceVar(modPreDec_PeakNdx_ChancePos, modPreDec_PeakNdx_ChancePosCount);
                [modPreDec_PeakVals_ChancePos] = AverageChanceVar(modPreDec_PeakVals_ChancePos, modPreDec_PeakVals_ChancePosCount);
                [modPreDec_PeakWids_ChancePos] = AverageChanceVar(modPreDec_PeakWids_ChancePos, modPreDec_PeakWids_ChancePosCount);
                [modPstDec_PeakNdx_ChancePos] = AverageChanceVar(modPstDec_PeakNdx_ChancePos, modPstDec_PeakNdx_ChancePosCount);
                [modPstDec_PeakVals_ChancePos] = AverageChanceVar(modPstDec_PeakVals_ChancePos, modPstDec_PeakVals_ChancePosCount);
                [modPstDec_PeakWids_ChancePos] = AverageChanceVar(modPstDec_PeakWids_ChancePos, modPstDec_PeakWids_ChancePosCount);
                clear modInt_PeakNdx_ChancePosCount modInt_PeakVals_ChancePosCount modInt_PeakWids_ChancePosCount...
                    modPreDec_PeakNdx_ChancePosCount modPreDec_PeakVals_ChancePosCount modPreDec_PeakWids_ChancePosCount...
                    modPstDec_PeakNdx_ChancePosCount modPstDec_PeakVals_ChancePosCount modPstDec_PeakWids_ChancePosCount
            end
        end
    
    end
    [hrSymTR_Full, hrSymDEC_Full, hrSymMN_Full] = mlb.CalcSymmetry_HR(xtdD_Real,xtdHR_ChanceFull);
    [hrSymTR_Pos, hrSymDEC_Pos, hrSymMN_Pos] = mlb.CalcSymmetry_HR(xtdD_Real,xtdHR_ChancePos);
    
    save(sprintf('%s\\xtd_mlb_wChance.mat', fileDir), '-v7.3');
    fprintf('File Saved @%s\n', datetime('now', 'Format', 'HH:mm:ss'));
end


%% Helper Functions
% Update Chance
function [chanceData, chanceCounter] = UpdateChanceVar(chanceData, chanceCounter, permData)
    for d = 1:numel(chanceData)
        chanceData{d} = chanceData{d} + permData{d};
        chanceCounter{d} = chanceCounter{d} + cellfun(@(a){~isnan(a)},permData);
    end
end
% Average Chance
function [chanceData] = AverageChanceVar(chanceData, chanceCounter)
    for d = 1:numel(chanceData)
        chanceData{d} = chanceData{d}./chanceCounter{d};
    end
end