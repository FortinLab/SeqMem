% PFC_SingleUnit_Modulation
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];

% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
%
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

%%
su = cell(1,length(fileDirs));
trlTAB = cell(1,length(fileDirs));
errTAB = cell(1,length(fileDirs));
rwdTAB = cell(1,length(fileDirs));
trlSTATS = cell(1,length(fileDirs));
errSTATS = cell(1,length(fileDirs));
rwdSTATS = cell(1,length(fileDirs));
trlDATA = cell(1,length(fileDirs));
errDATA = cell(1,length(fileDirs));
rwdDATA = cell(1,length(fileDirs));
trlIDS = cell(1,length(fileDirs));
errIDS = cell(1,length(fileDirs));
rwdIDS = cell(1,length(fileDirs));
errT = cell(1,length(fileDirs));
errTABseq = cell(1,length(fileDirs));
errSTATSseq = cell(1,length(fileDirs));
posInfo = cell(1,length(fileDirs));
for ani = 1:length(fileDirs)
    su{ani} = SingleUnit_SM(fileDirs{ani});
    [trlTAB{ani}, trlSTATS{ani}, trlDATA{ani}, trlIDS{ani}] = su{ani}.TrialPeriodSpiking;
    [errTAB{ani}, errSTATS{ani}, errDATA{ani}, errIDS{ani}] = su{ani}.ErrorSpiking;
    [rwdTAB{ani}, rwdSTATS{ani}, rwdDATA{ani}, rwdIDS{ani}] = su{ani}.RewardSpiking;
    [errT{ani}, errTABseq{ani}, errSTATSseq{ani}] = su{ani}.ErrorSpikingSequential;
    su{ani}.PlotUnitSummary;
    su{ani}.dsRate = 10;
    [posInfo{ani}, tsVect] = su{ani}.QuantPosInfo([-800 2000], 'PokeIn', 'pos');
end

%% Colormap setup
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);


%%
ensmblPosInfo = cell2mat(posInfo);
normPI = nan(size(ensmblPosInfo));
for u = 1:size(ensmblPosInfo,2)
    normPI(:,u) = ensmblPosInfo(:,u)./max(ensmblPosInfo(:,u));
end
normPI = normPI';
loc = nan(size(normPI,1),1);
for u = 1:size(normPI,1)
    loc(u) = find(normPI(u,:)==1,1,'first');
end
sortNormPI = sortrows([loc, normPI]);
sortNormPI = sortNormPI(:,2:end);
figure; 
imagesc(tsVect, 1:size(sortNormPI,1), sortNormPI);
colormap(cMap);

%%
% trlPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},trlTAB));
% errPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},errTAB));
% rwdPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2:4,end,:)))},rwdTAB));
% errTPvals = cell2mat(cellfun(@(a){[a.p]}, errT));
% errSeqPvals = cell2mat(cellfun(@(a){squeeze(cell2mat(a(2,end,:)))},errTABseq)')';