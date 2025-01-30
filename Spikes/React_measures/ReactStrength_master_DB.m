%%%% ReactStrength_master_DB
%
% This script calculates reactivation strength of neural ensembles
% derived from principal components calculated from correlation matrices
%
% References:
% - Peyrache et al., 2009, Nat Neuro
% - Peyrache et al., 2010, J Comput Neurosci
% - Tingley & Peyrache, 2020, Phil.Trans.B
%
% Please go inside the script and check the parameters

%% Parameters

nmouse = 994;
% Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = PathForExperimentsERC_DimaMAC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% Parameters of cross-correlograms
binsize = 0.1*1e4; % (measured in tsd units!!!);

% What to do?
Temp_Name = {'UMazeEpoch'};
Match_Name = {'PreSleep', 'PostSleep'};

% Do you want to save the figures?
sav = false; %%% Does not work now

%% Preallocate data
% Data
N = cell(length(Dir.path),1); % neuronal data
B = cell(length(Dir.path),1); % behavioral data
R = cell(length(Dir.path),1); % ripple data
S = cell(length(Dir.path),1); % sleep data
% Epochs
UMazeEpoch = cell(length(Dir.path),1);
CondEpoch = cell(length(Dir.path),1);
TaskEpoch = cell(length(Dir.path),1);
AfterConditioningEpoch = cell(length(Dir.path),1);
PreSleep = cell(length(Dir.path),1);
PostSleep = cell(length(Dir.path),1);
% Spike histograms
Q = cell(length(Dir.path),1);
DatTemplate = cell(length(Dir.path),length(Temp_Name));
DatMatch = cell(length(Dir.path),length(Match_Name));
TimeTemplate = cell(length(Dir.path),length(Temp_Name));
TimeMatch = cell(length(Dir.path),length(Match_Name));
% Results
R = cell(length(Dir.path), length(Temp_Name), length(Match_Name));
phi = cell(length(Dir.path), length(Temp_Name), length(Match_Name));


%% Get Data
for imouse=1:length(Dir.path)
    N{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'],'S','PlaceCells');
    % If there are less than 2 PCs - don't do
    if isfield(N{imouse}.PlaceCells,'idx')
        if length(N{imouse}.PlaceCells.idx)>2
            B{imouse} = load([Dir.path{imouse}{1} 'behavResources.mat'],'SessionEpoch', 'CleanVtsd', 'FreezeAccEpoch');
            R{imouse} = load([Dir.path{imouse}{1} 'Ripples.mat'],'ripples');
            if strcmp(Dir.name{imouse}, 'Mouse906') || strcmp(Dir.name{imouse}, 'Mouse977') % Mice with bad OB-based sleep scoring
                S{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % Sleep is not used
            else
                S{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');  % Sleep is not used
            end
            
        end
    end
end

%% Create auxiliary epochs (ERC-MOBS specific)
for imouse = 1:length(Dir.path)
    
    [~, UMazeEpoch{imouse}, CondEpoch{imouse}, TaskEpoch{imouse}, AfterConditioningEpoch{imouse}] = ReturnMnemozyneEpochs(B{imouse}.SessionEpoch,...
        'SpeedData', B{imouse}.CleanVtsd, 'SpeedThresh', 5);
    PreSleep{imouse} = and(B{imouse}.SessionEpoch.PreSleep, S{imouse}.SWSEpoch);
    PostSleep{imouse} = and(B{imouse}.SessionEpoch.PostSleep, S{imouse}.SWSEpoch);
end

%% Create template and match epochs

for imouse = 1:length(Dir.path)
    for itemp = 1:length(Temp_Name)
        eval(['template{imouse, itemp} = ' Temp_Name{itemp} '{imouse};']);
    end
    for imatch = 1:length(Match_Name)
        eval(['match{imouse, imatch} = ' Match_Name{imatch} '{imouse};']);
    end
end

%% Make binned spike histograms

for imouse = 1:length(Dir.path)
    Q{imouse} = MakeQfromS(N{imouse}.S,binsize);
    
    % Template epochs
    for itemp = 1:length(Temp_Name)
        QTemplate = Restrict(Q{imouse},template{imouse, itemp});
        DatTemplate{imouse, itemp} = zscore(full(Data(QTemplate)));
        TimeTemplate{imouse, itemp} = Range(QTemplate);
    end
    
    % Matching epochs
    for imatch = 1:length(Match_Name)
        QMatch = Restrict(Q{imouse}, match{imouse, imatch});
        DatMatch{imouse, imatch} = zscore(full(Data(QMatch)));
        TimeMatch{imouse, imatch} = Range(QMatch);
    end
    
end

%% Calculate reactivation strength
for imouse = 1:length(Dir.path)
    for itemp = 1:length(Temp_Name)
        for imatch = 1:length(Match_Name)
            [R{imouse, itemp, imatch},phi{imouse, itemp, imatch}] = ReactStrength(DatTemplate{imouse, itemp},DatMatch{imouse, imatch},'Method','ICA');
        end
    end
end

%% Visally check for patterns (first two are not reactivated maybe)

%% Mean reactivation strength Pre/Post
MeanPre = zeros(size(R{1,1,1},2),1);
MeanPost = zeros(size(R{1,1,2},2),1);
SkewPre = zeros(size(R{1,1,1},2),1);
SkewPost = zeros(size(R{1,1,2},2),1);
for i = 1:length(MeanPre)
    MeanPre(i) = mean(R{1,1,1}(:,i));
    SkewPre(i) = skewness(R{1,1,1}(:,i));
end
for i = 1:length(MeanPost)
    MeanPost(i) = mean(R{1,1,2}(:,i));
    SkewPost(i) = skewness(R{1,1,2}(:,i));
end

subplot(121)
scatter(MeanPre(2:end), MeanPost(2:end))
hold on
line([0 0.5], [0 0.5])
hold off
title('Mean PreVsPost')
subplot(122)
scatter(SkewPre, SkewPost)
% hold on
% line([0 0.35], [0 0.5])
title('Skewness PreVsPost')
xlim([0 20])
hold on
line([0 20], [0 20])
hold off

%% Over time put shaded SWSEpochs

PreSleepInt = [Start(PreSleep{imouse}) End(PreSleep{imouse})];
% PreSleepInt = PreSleepInt - PreSleepInt(1,1);
PostSleepInt = [Start(PostSleep{imouse}) End(PostSleep{imouse})];
% PostSleepInt = PostSleepInt - PostSleepInt(1,1);

figure
subplot(121)
plot(TimeMatch{1,1}, R{1,1,1}(:,3));
yl = [0 100];
ylim(yl);
hold on
for iep = 1:length(PreSleepInt)
    patch([PreSleepInt(iep,1) PreSleepInt(iep,2)], [yl(2)*0.1 yl(2)*0.9], 'g', 'FaceAlpha', 0.3)
end
hold off
xlim([1.8e7 6e7])
subplot(122);
plot(TimeMatch{1,2}, R{1,1,2}(:,3));
ylim([0 100])
hold on
for iep = 1:length(PreSleepInt)
    patch([PostSleepInt(iep,1) PostSleepInt(iep,2)], [yl(2)*0.1 yl(2)*0.9], 'g', 'FaceAlpha', 0.3)
end
hold off


%% Distribution - check skewness (plot skewness pre/post)
figure
subplot(121)
histogram(R{1,1,1}(:,4), 'Normalization', 'pdf');
set(gca, 'YScale', 'log')
xlim([-50 250])
subplot(122);
histogram(R{1,1,2}(:,4), 'Normalization', 'pdf');
xlim([-50 250])
set(gca, 'YScale', 'log')