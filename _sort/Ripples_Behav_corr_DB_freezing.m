function EV_4corr = Ripples_Behav_corr_DB_freezing(mice, expe, varargin)
%
% The function calculates Spearman rank correlation between explained
% variance and two behavioral parameters:
%   - Shock zone occupancy
%   - Latency to enter the shock zone
%
% INPUT
%
%   mice                array with mice numbers that will harvested from ERC
%                       PathForExperiments. Each mouse should contain
%                       PlaceCells structure in its SpikeData.mat
%   expe                type of experiment to analyse ('PAG' or 'MFB')
%   IsSave              if true, saves figures in dropbox (default = false)(optional)
%
% OUTPUT
%
%   EV_4corr            array of EV values used for analysis
%   Diff_occup          array of SZ occupancy values used for analysis
%   Diff_entrytime      array of SZ entry latency values used for analysis
%
% LIST OF EPOCHS PERMITTED
%
%   PAG expe           'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'RipplesEpoch', 'PostTests'
%   MFB expe           'Explo', 'CondMov', 'FullTask', 'RipplesEpoch', 'PostTests'
%   Novel expe         'Explo', 'CondMov', 'FullTask', 'RipplesEpoch'
%
% EXAMPLE
%
%   [EV_4corr, Diff_occup, Diff_entrytime] = EV_Behav_corr_DB_new(mice, expe, EV_epoch);
%   [EV_4corr, Diff_occup, Diff_entrytime] = EV_Behav_corr_DB_new(mice, expe, EV_epoch, 'SleepType', 'REM');
%
% By Dima Bryzgalov, MOBS team, Paris,
% 15/09/2021
% github.com/bryzgalovdm

%% Parameters
type_sleep = 'NREM';
sav = false;

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help EV_Behav_corr_DB_new'' for details).');
            end
    end
end


%% Calculate behavioral results
% Manage experiment
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';
end

% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);

%% Calculate Ripples
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        try
            Rip{cnt} = load([Dir.path{imouse}{isession} 'SWR.mat'], 'ripples');
        catch
            Rip{cnt} = load([Dir.path{imouse}{isession} 'Ripples.mat'], 'ripples');
        end
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'], 'behavResources', 'SessionEpoch', 'FreezeAccEpoch', 'ZoneEpoch');
        try
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Wake');
        catch
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Wake'); % REM and Sleep are not used
        end
        cnt=cnt+1;
    end
end

% Prepare intervalSets for ripples
for i=1:length(b)
    if strcmp(expe, 'PAG')
        [~, ~, CondEpoch{i}, ~, ~] = ReturnMnemozyneEpochs(b{i}.SessionEpoch);
        CondFreezeEpoch{i} = and(CondEpoch{i}, b{i}.FreezeAccEpoch);
    else
        [~, ~, CondEpoch{i}, ~, ~] = ReturnMnemozyneEpochs(b{i}.SessionEpoch);
    end
    try
        PreNREMEpoch{i} = and(RestrictToTime(b{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleep{i}.SWSEpoch);
        PostNREMEpoch{i} = and(RestrictToTime(b{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleep{i}.SWSEpoch);
    catch
        PreNREMEpoch{i} = and(b{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch);
        PostNREMEpoch{i} = and(b{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch);
    end
end

% Restrict sleeps to first 30 min
for i = 1:length(b)
    try
        EarlyNREM_pre{i} = RestrictToTime(PreNREMEpoch{i}, 30*60*1e4);
        EarlyNREM_post{i} = RestrictToTime(PostNREMEpoch{i}, 30*60*1e4);
    catch
        EarlyNREM_pre{i} = PreNREMEpoch{i};
        EarlyNREM_post{i} = PostNREMEpoch{i};
    end
end

% Extract ripples for presleep and postsleep (during detected sleep only)
for i = 1:length(b)
    ripplesPeak=ts(Rip{i}.ripples(:,2)*1e4);
    PreRipples{i}=Restrict(ripplesPeak,PreNREMEpoch{i});
    PostRipples{i}=Restrict(ripplesPeak,PostNREMEpoch{i});
    if strcmp(expe, 'PAG')
        CondRipples{i}=Restrict(ripplesPeak,CondFreezeEpoch{i});
    else
        CondRipples{i}=Restrict(ripplesPeak,CondEpoch{i});
    end
    
    
    PreR_density(i) = length(Range(PreRipples{i})) / tot_length(PreNREMEpoch{i}, 's');
    PostR_density(i) = length(Range(PostRipples{i})) / tot_length(PreNREMEpoch{i}, 's');
    if strcmp(expe, 'PAG')
        CondR_density(i) = length(Range(CondRipples{i})) / tot_length(CondFreezeEpoch{i}, 's');
    else
        CondR_density(i) = length(Range(CondRipples{i})) / tot_length(CondEpoch{i}, 's');
    end
    
    EarlyPreR_density(i) = length(Restrict(ripplesPeak,EarlyNREM_pre{i})) / tot_length(EarlyNREM_pre{i}, 's');
    EarlyPostR_density(i) = length(Restrict(ripplesPeak,EarlyNREM_post{i})) / tot_length(EarlyNREM_post{i}, 's');
    
end

%% Freezing
% Prepare epochs
for isession = 1:length(b)
    [~, ~, ~, CondEpoch{isession}, ~, TestPostEpoch{isession}] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    FreezeEpoch{isession} = and(sleep{isession}.Wake, and(b{isession}.FreezeAccEpoch, CondEpoch{isession}));
    FreezeEpochShock{isession} = and(FreezeEpoch{isession}, b{isession}.ZoneEpoch.Shock);
    FreezeEpochSafe{isession} = and(FreezeEpoch{isession}, b{isession}.ZoneEpoch.NoShock);
    
    FreezeEpochPost{isession} = and(sleep{isession}.Wake, and(b{isession}.FreezeAccEpoch, TestPostEpoch{isession}));
    FreezeEpochPostShock{isession} = and(FreezeEpochPost{isession}, b{isession}.ZoneEpoch.Shock);
    FreezeEpochPostSafe{isession} = and(FreezeEpochPost{isession}, b{isession}.ZoneEpoch.NoShock);
end

% Calculate percentage of freezing
for isession = 1:length(b)
    OverallFreeze(isession) = tot_length(FreezeEpoch{isession})/tot_length(CondEpoch{isession})*100;
    ShockFreeze(isession) = tot_length(FreezeEpochShock{isession})/tot_length(CondEpoch{isession})*100;
    SafeFreeze(isession) = tot_length(FreezeEpochSafe{isession})/tot_length(CondEpoch{isession})*100;
    
    OverallFreezePost(isession) = tot_length(FreezeEpochPost{isession})/tot_length(TestPostEpoch{isession})*100;
    ShockFreezePost(isession) = tot_length(FreezeEpochPostShock{isession})/tot_length(TestPostEpoch{isession})*100;
    SafeFreezePost(isession) = tot_length(FreezeEpochPostSafe{isession})/tot_length(TestPostEpoch{isession})*100;
end

% % Find indices of PreTests and PostTest session in the structure
% id_Pre = cell(1,length(a));
% id_Post = cell(1,length(a));
% 
% for i=1:length(a)
%     id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
%     id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
% end
% 
% % Calculate average occupancy
% occup_shock = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post
% for i=1:length(a)
%     for k=1:length(id_Pre{i})
%         temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k)));
%         occup_shock(i,k,1) = temp(1);
%         temp = CalculateZoneOccupancy(a{i}.behavResources(id_Post{i}(k)));
%         occup_shock(i,k,2) = temp(1);
%     end
% end
% 
% occup_shock_mean = nan(length(Dir.path), 2);
% for izone = 1:2 % 1 codes for preTest, 2 for postTest
%     occup_shock_mean(:,izone) = mean(squeeze(occup_shock(:,:,izone)),2);
% end
% 
% % Prepare the 'first enter to shock zone' array
% EntryTime_shock = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post
% for i = 1:length(a)
%     for k=1:length(id_Pre{i})
%         temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Pre{i}(k)), 240);
%         EntryTime_shock(i,k,1) = temp(1);
%         temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Post{i}(k)), 240);
%         EntryTime_shock(i,k,2) = temp(1);
%     end  
% end
%     
% EntryTime_shock_mean = nan(length(Dir.path), 2);
% for izone = 1:2 % 1 codes for preTest, 2 for postTest
%     EntryTime_shock_mean(:,izone) = mean(squeeze(EntryTime_shock(:,:,izone)),2);
% end
% 
% Diff_occup = diff(occup_shock_mean,1,2)*100;
% Diff_entrytime = diff(EntryTime_shock_mean,1,2);

%% Calculate freezing

%% Choose epochs
% num_epoch_wake = strcmp(states_wake, EV_epoch);
% num_epoch_sleep = strcmp(states_sleep, type_sleep);
% EV_4corr = nonzeros(EV{num_epoch_sleep}{num_epoch_wake});
% EV_4corr(isnan(EV_4corr)) = [];

%% Calculate correlations
% Occupancy
[rho, p] = corr(PostR_density', OverallFreeze', 'Type', 'Spearman'); 
r(1) = rho;
pvalue(1) = p;

% Entries
[rho, p] = corr(EarlyPostR_density', ShockFreeze', 'Type', 'Spearman'); 
r(2) = rho;
pvalue(2) = p;

%% Plot figures
f= figure('units', 'normalized', 'outerposition', [0 0 1 0.65]);
subplot(121)
scatter(PostR_density', OverallFreeze', 150, 'k', 'filled')
t1 = text(.4, .9, ['rho = ' num2str(round(r(1), 3))], 'sc', 'FontSize', 13);
t2 = text(.4, .8, ['p = ' num2str(round(pvalue(1), 3))], 'sc', 'FontSize', 13);
l1 = lsline;
l1.Color = 'k';
l1.LineWidth = 1.5;
xlabel('Ripples density in PostSleep (rip/s)');
ylabel('Overall freezing % during Cond');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(gca, 'LineWidth', 1);
makepretty_DB

subplot(122)
scatter(EarlyPostR_density', SafeFreeze', 150, 'k', 'filled')
t1 = text(.4, .9, ['rho = ' num2str(round(r(2), 3))], 'sc', 'FontSize', 13);
t2 = text(.4, .8, ['p = ' num2str(round(pvalue(2), 3))], 'sc', 'FontSize', 13);
l1 = lsline;
l1.Color = 'k';
l1.LineWidth = 1.5;
xlabel('Ripples density in PostSleep in the first 30 min (rip/s)');
ylabel('SafeZone freezing % during Cond');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(gca, 'LineWidth', 1);
makepretty_DB

if sav
    foldertosave = ChooseFolderForFigures_DB('ReactReplay');
    saveas(f,[foldertosave filesep 'EV' filesep 'Correlations' filesep 'CorrRip_Freeze.fig']);
    saveFigure(f,'CorrRip_Freeze', [foldertosave filesep 'EV' filesep 'Correlations']);
end

end