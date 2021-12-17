clear all; %#ok<CLALL>
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Dual_List\GE11_Session146'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE13_Session103'},...
    {'D:\WorkBigDataFiles\PFC\Dual_List\GE17_Session110'}];

% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

binSize = 200;
dsRate = 50;
trlWindow = {[-800 500]; [-500 500]};
alignment = [{'PokeIn'}, {'PokeOut'}];
lfpWindow = [16 32];
numPerms = 10;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];