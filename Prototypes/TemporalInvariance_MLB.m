function [odrPost, odrDecode, trialInfo, fisc, meanPostOdr, meanDecodeOdr, trlLFPphase, trlLFPpower, trlTimeVect] = TemporalInvariance_MLB(smDir, binSize, dsRate, trlWindow)
%     clc;
    alignment = 'PokeIn';
    
    if nargin==0
        smDir = uigetdir;
    end
    mlb = MLB_SM(smDir);
    if nargin == 0
        mlb.binSize = 100;
        mlb.dsRate = 50;    
        trlWindow = [-800 1200];
    else
        mlb.binSize = binSize;
        mlb.dsRate = dsRate;
    end
    
    rng(sum(clock));

    %%
    [trlSpikes, trlTimeVect] = mlb.PP_TrialMatrix_Spiking(trlWindow, alignment);
    [trlLFPphase, trlLFPpower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow, alignment);
    mlb.PP_IdentifyFISCseqs;
    trialInfo = mlb.trialInfo;
    fisc = mlb.fiscTrials;

    odrPost = cell(size(trlSpikes,3),1);
    odrDecode = nan(length(trlTimeVect),length(trlTimeVect),length(mlb.trialInfo));
    % For each trial
    for trl = 1:size(trlSpikes,3)
        %% Create the likelihoods
        % Pull out the FISC trials...
        fiscTrls = mlb.fiscTrials;
        % ... and remove the whole Sequence, if the current trial is a part of it
        if sum(sum(trl==mlb.fiscTrials))>=1
            seqLog = sum(trl==mlb.fiscTrials)==1;
            fiscTrls(:,seqLog) = [];
        else
        end
        %... then create templates for each odor based on the FISC trials (i.e. calculate their average PSTHs)
        fiscLikes = cell(1,4);
        for odr = 1:size(fiscTrls,1)
            tempLike = nan(size(trlSpikes,1), size(trlSpikes,2), size(fiscTrls,2));
            for seq = 1:size(fiscTrls,2)
                tempLike(:,:,seq) = trlSpikes(:,:,fiscTrls(odr,seq));
            end
            fiscLikes{odr} = mean(tempLike,3,'omitnan');
        end
        %% Use the likelihoods
        % Step through each trial
        curObsv = trlSpikes(:,:,trl);
        tempPost = nan(size(trlSpikes,1),size(trlSpikes,1),mlb.seqLength);
        for likeTime = 1:size(trlTimeVect,1)
            tempLike = cell2mat(cellfun(@(a){a(likeTime,:)}, fiscLikes)');
            curPosts = mlb.CalcStaticBayesPost(tempLike, curObsv);
            tempPost(:,likeTime,:) = curPosts;
            for t = 1:size(curPosts,1)
                decode = find(curPosts(t,:)==max(curPosts(t,:)));
                select = rand(1,length(decode));
                odrDecode(t,likeTime,trl) = decode(select==max(select));
            end
        end
        odrPost{trl} = tempPost;        
        fprintf('%i\n', trl);
    end
    %% Now Plot stuff I guess?
    odorVect = [mlb.trialInfo.Odor];
    posVect = [mlb.trialInfo.Position];
    isLog = odorVect==posVect;
    perfLog = logical([mlb.trialInfo.Performance]);

    meanPostOdr = cell(mlb.seqLength);
    meanDecodeOdr = cell(mlb.seqLength);
    for trlOdr = 1:mlb.seqLength
        tempFig1 = figure;
        tempFig2 = figure;
        trlVect = mlb.fiscTrials(trlOdr,:);
%         trlVect = odorVect==trlOdr & isLog & perfLog;
        for decodeOdr = 1:mlb.seqLength
            figure(tempFig1);
            subplot(3,3,decodeOdr);
            meanPost = mean(cell2mat(reshape(cellfun(@(a){a(:,:,decodeOdr)}, odrPost(trlVect)), [1,1,sum(trlVect~=0)])),3,'omitnan');
            meanPostOdr{trlOdr, decodeOdr} = meanPost;
            imagesc(trlTimeVect,trlTimeVect,meanPost, [0 0.75]);
            title(sprintf('Posteriors: Decode Pos%i During Pos%i', decodeOdr, trlOdr));
            xlabel('Template Time');
            ylabel('Decode Time');
            set(gca, 'ydir', 'normal');
            figure(tempFig2)
            subplot(3,3,decodeOdr);
            meanDecode = mean(odrDecode(:,:,trlVect)==decodeOdr,3,'omitnan');
            meanDecodeOdr{trlOdr,decodeOdr} = meanDecode;
            imagesc(trlTimeVect, trlTimeVect, meanDecode, [0 0.75]);
            title(sprintf('Decoding: Decode Pos%i During Pos%i', decodeOdr, trlOdr));
            xlabel('Template Time');
            ylabel('Decode Time');
            set(gca, 'ydir', 'normal');
        end
        drawnow
    end
end