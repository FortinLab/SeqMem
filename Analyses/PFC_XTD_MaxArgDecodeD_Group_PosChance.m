% PFC_XTD_TrialWindows_Group_PosChance

%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\GE24_Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'}
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
%     {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];
% % CA1 Data
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% tets = [1,22,17,18,17]; % Lateral/Distal
% % tets = [7,3,1,5,5]; % Medial/Proximal

fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'}];
binSize = 200;
dsRate = 50;
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts
% 
% alignments = [{'PokeIn'}, {'PokeOut'}];
% trlWindows = [{[-1200 2000]}, {[-2000 1200]}];

alignments = {'PokeIn'};
% trlWindows = {[-1200 2500]}; %% Overlapping Pre/Post ITD
trlWindows = {[-700 2000]}; %% Non-Overlapping Pre/Post ITD

% alignments = {'PokeOut'};
% trlWindows = {[-1500 1500]};
    
lfpWindow = [16 32];
numChancePerms = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Chance Perm Nums

%% Create the mlb objects
realTic = tic;
mlb = cell(size(fileDirs));
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb{ani} = MLB_SM(fileDirs{ani});
end

%% Run Cross-Temporal Decoding using MaxArg decoding
xtdTrlPosts = cell(length(alignments),length(fileDirs));
xtdTrlDecodes = cell(length(alignments),1,length(fileDirs));
xtdTrlIDs = cell(1,length(fileDirs));
xtdTSvects = cell(1,length(alignments));
for ani = 1:length(fileDirs)  
    mlb{ani}.binSize = binSize;
    mlb{ani}.dsRate = dsRate;
    mlb{ani}.bayesType = bayesType;
    for al = 1:length(alignments)
        mlb{ani}.alignments = alignments(al);
        mlb{ani}.windows = trlWindows(al);
        mlb{ani}.SetLikes_ISC;
        mlb{ani}.Process_IterativeLikelyL1O;
        xtdTrlPosts{al,ani} = permute(mlb{ani}.post, [1,3,2]);
        tempDecodes = cellfun(@(a){mlb{ani}.DecodeBayesPostNew(a,3)}, mlb{ani}.post);
        xtdTrlDecodes{al,1,ani} = cell(size(mlb{ani}.post,1),1);
        for pos = 1:size(mlb{ani}.post,1)
            xtdTrlDecodes{al,1,ani}{pos} = cell2mat(tempDecodes(pos,:,~isnan(mlb{ani}.postTrlIDs(pos,:))));
        end
        xtdTSvects{al} = mlb{ani}.obsvTimeVect;
    end
    xtdTrlIDs{ani} = mlb{ani}.postTrlIDs;   
    fprintf('Ani#%i done\n',ani);  
end

%% Collapse Trial Across Animals
trialXTD = cell(1,mlb{ani}.seqLength);
trialHFR = cell(mlb{ani}.seqLength);
for p = 1:mlb{ani}.seqLength
    trialXTD{p} = cell2mat(cellfun(@(a)a(p),xtdTrlDecodes));
    for pp = 1:mlb{ani}.seqLength
        trialHFR{p,pp} = trialXTD{p}==pp;
    end        
end
figure;
trialD = cell(mlb{ani}.seqLength);
for p = 1:mlb{ani}.seqLength
    posDecode = trialHFR(:,p);
    for pp = 1:mlb{ani}.seqLength
        tempHR = mean(posDecode{pp},3);
        tempFAR = mean(cell2mat(permute(posDecode(1:mlb{ani}.seqLength~=pp),[2,3,1])),3);
        probs = sort([tempHR(:);tempFAR(:)]);
        minProb = probs(find(probs>0,1,'first'));
        maxProb = probs(find(probs<1,1,'last'));
        tempHR(tempHR==0) = minProb;
        tempHR(tempHR==1) = maxProb;
        tempFAR(tempFAR==0) = minProb;
        tempFAR(tempFAR==1) = maxProb;
        minFAR = min(sort(tempFAR(:)));
        trialD{p,pp} = arrayfun(@(a,b)norminv(a)-norminv(b),tempHR,tempFAR);
        subplot(mlb{ani}.seqLength, mlb{ani}.seqLength, sub2ind([mlb{ani}.seqLength, mlb{ani}.seqLength],pp,p));
%         imagesc(zscore(trialD{p,pp}',0, 'all'),[-1 1]);
        imagesc(trialD{p,pp}',[-1 1]);
        ylabel('Training Time');
        xlabel('Testing Time');
        title(sprintf('Pos=%i, D=%i',p,pp));
        set(gca,'ydir', 'normal');
    end
end

        
%% Calculate 
xtdTrlDpos = cell(1,length(alignments));
xtdTrlDodr = cell(1,length(alignments));
xtdAcc = cell(size(mlb{ani}.post,1), size(mlb{ani}.post,1),2);
for al = 1:length(alignments)
    tempDecodes = [xtdTrlDecodes{al,:,:}];
    posLogs = cell(size(tempDecodes));
    odrLog = cell(size(tempDecodes));
    for p = 1:numel(posLogs)
        [r,~] = ind2sub(size(tempDecodes),p);
        if r > 4
            posLogs{p} = ones(size(tempDecodes{p},3),1)*(r-4);
        else
            posLogs{p} = ones(size(tempDecodes{p},3),1)*r;
        end
        odrLog{p} = ones(size(tempDecodes{p},3),1)*r;
    end
    xtdTrlDpos{al}= mlb{ani}.CalcDprmMtxFromDecode(cell2mat(permute(tempDecodes(:), [2,3,1])), cell2mat(posLogs(:)));
    xtdTrlDodr{al}= mlb{ani}.CalcDprmMtxFromDecode(cell2mat(permute(tempDecodes(:), [2,3,1])), cell2mat(odrLog(:)));
%     for posT = 1:mlb{ani}.seqLength
%         curPosTrl = cell2mat(permute(tempDecodes(posT,:), [1,3,2]));
%         for posD = 1:size(mlb{ani}.post,1)
%             xtdAcc{posT,posD,al} = mean(curPosTrl==posD,3,'omitnan');
%         end
%     end    
%     figure;
%     tempACC = xtdAcc(:,:,al);
%     for p = 1:numel(tempACC)
%         subplot(5,5,p); 
%         imagesc(xtdTSvects{al},xtdTSvects{al},tempACC {p}');
%         set(gca, 'ydir', 'normal'); 
%         xlabel('Obsv Time'); 
%         ylabel('Train Time'); 
%     end
end

for al = 1:length(alignments)
    figure;
    for p1 = 1:size(xtdTrlDpos{al},1)
        for p2 = 1:size(xtdTrlDpos{al},2)
            subplot(size(xtdTrlDpos{al},1),size(xtdTrlDpos{al},2),sub2ind([size(xtdTrlDpos{al},1),size(xtdTrlDpos{al},2)],p2,p1));
            imagesc(xtdTSvects{al},xtdTSvects{al},xtdTrlDpos{al}{p1,p2}', [-1 1]);
            set(gca, 'ydir', 'normal');
            xlabel('Obsv Time');
            ylabel('Train Time');
        end
    end
    figure;
    for p1 = 1:size(xtdTrlDodr{al},1)
        for p2 = 1:size(xtdTrlDodr{al},2)
            subplot(size(xtdTrlDodr{al},1),size(xtdTrlDodr{al},2),sub2ind([size(xtdTrlDodr{al},1),size(xtdTrlDodr{al},2)],p2,p1));
            imagesc(xtdTSvects{al},xtdTSvects{al},xtdTrlDodr{al}{p1,p2}', [-1 1]);
            set(gca, 'ydir', 'normal');
            xlabel('Obsv Time');
            ylabel('Train Time');
        end
    end
end

diff = cell(1,8);
for pos = 1:8
    if pos>4
        diff{pos} = xtdTrlDodr{al}{pos,pos}-xtdTrlDodr{al}{pos,pos-4};
    else
        diff{pos} = xtdTrlDodr{al}{pos,pos}-xtdTrlDodr{al}{pos,pos+4};
    end
    subplot(1,8,pos)
    imagesc(xtdTSvects{al},xtdTSvects{al},diff{pos},[-0.1 0.1]);
    set(gca,'ydir','normal');
end

