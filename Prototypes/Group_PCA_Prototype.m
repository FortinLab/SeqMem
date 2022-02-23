
fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
binSize = 10;
dsRate = 10;
trlWindow = {[-500 500]};
alignment = {'PokeOut'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
ssProportion = 0.5;
ssType = 0; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

windowSteps = 100;

for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.ssProportion = ssProportion;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    mlb.SetLikes_FISC;
    mlb.SetObserves_Session;
    mlb.PFC_PCA;
end