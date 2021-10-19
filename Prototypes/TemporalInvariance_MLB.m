function TemporalInvariance_MLB(smDir)
trlWindow = [-800 1200];

if nargin==0
    smDir = uigetdir;
end
mlb = MLB_SM(smDir);
mlb.binSize = 100;
mlb.dsRate = 5;
rng(sum(clock));

%%
[trlSpikes, trlTimeVect] = mlb.PP_TrialMatrix_Spiking(trlWindow, 'PokeIn');
mlb.PP_IdentifyFISCseqs;

odrPost = cell(size(trlSpikes,3),1);
% For each trial
for trl = 1:size(trlSpikes,3)
    tempPost = nan(size(trlSpikes,1),size(trlSpikes,1),4);
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
    for likeTime = 1:size(trlTimeVect,1)
        tempLike = cell2mat(cellfun(@(a){a(likeTime,:)}, fiscLikes)');
        tempPost(:,likeTime,:) = mlb.CalcStaticBayesPost(tempLike, curObsv);
    end
    odrPost{trl} = tempPost;
%     fprintf('%i\n', trl);
    %% Now Plot stuff I guess?
    odorVect = [mlb.trialInfo.Odor];
    posVect = [mlb.trialInfo.Position];
    isLog = odorVect==posVect;
    perfLog = logical([mlb.trialInfo.Performance]);
    
    for trlOdr = 1:4
        figure;
%         tempTrls = odrPost(odorVect==trlOdr & isLog & perfLog);
        tempTrls = odrPost(mlb.fiscTrials(trlOdr,:));
        for decodeOdr = 1:4
            subplot(2,2,decodeOdr);
            imagesc(trlTimeVect,trlTimeVect,mean(cell2mat(reshape(cellfun(@(a){a(:,:,decodeOdr)}, tempTrls), [1,1,length(tempTrls)])),3,'omitnan'), [0 0.75]);
            title(sprintf('Decode Pos%i During Pos%i', decodeOdr, trlOdr));
            xlabel('Train Time');
            ylabel('Decode Time');
            set(gca, 'ydir', 'normal');
        end
    end
end


    

end