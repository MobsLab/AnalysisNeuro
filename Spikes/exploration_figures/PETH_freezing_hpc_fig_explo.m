clear all

%% Parameters
IsSave=0;
mice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168 1182];
wi = 2; % in sec
binsize = 0.05; % in sec
nbins = 81;
speed_thresh = 5; % in cm/s

%% Get data
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

% Allocate arrays
b = cell(length(Dir.path),1); % behavior
s = cell(length(Dir.path),1); % spikes
CondEpoch = cell(length(Dir.path),1);
UMazeEpoch = cell(length(Dir.path),1);
FreezeEpoch = cell(length(Dir.path),1);
Neurons = cell(length(Dir.path),1);
PlaceCells = nan(5e4,1);
SpInfo = nan(5e4,1);
FR = nan(5e4,1);
NeuronClass = nan(5e4,1);

for imouse = 1:length(Dir.path)
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch', 'FreezeAccEpoch', 'MovAcctsd', ...
        'CleanAlignedXtsd', 'CleanAlignedYtsd', 'CleanVtsd');
    s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells');
end

%% Organize epochs and data
% Epochs
for imouse = 1:length(Dir.path)
    [~, ~, CondEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch);
    [~, UMazeEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch,...
        'Speed', b{i}.CleanVtsd, 'SpeedThresh', speed_thresh);
    
    FreezeEpoch{imouse} = and(b{imouse}.FreezeAccEpoch,CondEpoch{imouse});
    FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
    FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4);
end
% Data
for imouse = 1:length(Dir.path)
    Neurons{imouse} = cell(length(s{imouse}.S),1);
    Q = MakeQfromS(s{imouse}.S, binsize*1e4);
    time = Restrict(Q, CondEpoch{imouse});
    tempdat = zscore(full(Data(Restrict(Q, CondEpoch{imouse}))));
    for icell = 1:length(s{imouse}.S)
        Neurons{imouse}{icell} = [Range(time)/1e4 tempdat(:,icell)];
    end
end
    
%% Calculate PETH
cnt=1;
cnt_pc=1;
for imouse = 1:length(Dir.path)
    for icell = 1:length(Neurons{imouse})
        [r,i] = Sync(Neurons{imouse}{icell}, Start(FreezeEpoch{imouse})/1e4 ,'durations',[-wi wi]);
        T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
        OnResp(cnt,:) = mean(T);
        
        [r,i] = Sync(Neurons{imouse}{icell}, End(FreezeEpoch{imouse})/1e4 ,'durations',[-wi wi]);
        T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
        OffResp(cnt,:) = mean(T);
        
        % Calculate auxiliary info
        if IsCalcAdditional
            % Spatial Info
            [~, ~, stats] = PlaceField_DB(s{imouse}.S{icell},...
                b{imouse}.CleanAlignedXtsd, b{imouse}.CleanAlignedYtsd, 'Epoch', UMazeEpoch{imouse},...
                'PlotResults',0, 'PlotPoisson',0);
            SpInfo(cnt) = stats.spatialInfo;
            % Place Cell or not
            if sum(s{imouse}.PlaceCells.idx == icell) > 0
                PlaceCells(cnt) = true;
            else
                PlaceCells(cnt) = false;
            end
            % Firing rate
            FR(cnt) = s{imouse}.BasicNeuronInfo.firingrate(icell);
            % Neuron class
            NeuronClass(cnt) = s{imouse}.BasicNeuronInfo.neuroclass(icell);
        end
        
        cnt = cnt + 1;
        
    end
    
end

if IsCalcAdditional
    SpInfo(isnan(SpInfo)) = [];
    PlaceCells(isnan(PlaceCells)) = [];
    FR(isnan(FR)) = [];
    NeuronClass(isnan(NeuronClass)) = [];
end

%% Calculate accelerometer trace
Accelero = nan(length(Dir.path), 201*2);
for imouse = 1:length(Dir.path)
    M_st = PlotRipRaw(b{imouse}.MovAcctsd,Start(FreezeEpoch{imouse},'s'),2000,1,0,0);
    M_end = PlotRipRaw(b{imouse}.MovAcctsd,End(FreezeEpoch{imouse},'s'),2000,1,0,0);
    
    time = [(M_st(:,1))' (M_end(:,1)+4)'];
    Accelero(imouse,:) = [(M_st(:,2))' (M_end(:,2))'];
end
    

%% Plot figure
f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.4 0.8]);
DatNormZ = [OnResp,OffResp];
SustVal = nanmean(DatNormZ(:,60:100),2);
% SustVal = nanmean(DatNormZ(:,121:141),2);
% SustVal = nanmean(DatNormZ(:,40:120),2)./abs(nanmean(DatNormZ(:,121:141),2));
UseForTresh = SustVal;
[a,ind_sort] = sort(UseForTresh);
subplot(4,1,1:3)
imagesc(-2:0.05:6,1:size(DatNormZ,1),DatNormZ(ind_sort,:))
caxis([-0.6 0.6])
line([0 0],ylim,'color','k')
line([4 4],ylim,'color','k')
line(xlim, [100 100], 'Color',[.9856, .7372, .2537], 'LineStyle', '--');
line(xlim, [600 600], 'Color', 'm', 'LineStyle', '--')
makepretty
ylabel('Neuron #')
xlabel('Time to freeze on (s)')
% Accelero meter
subplot(4,1,4)
if size(Accelero,1) > 1
    shadedErrorBar(time,mean(Accelero, 1), std(Accelero, 1), 'k')
else
    plot(time,Accelero,'color','k')
end
makepretty
ylabel('Accelero')
xlabel('Time to freeze on (s)')

% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1,[foldertosave '/Freezing/PETHonFreezing_hpc.fig']);
    saveFigure(f1, 'PETHonFreezing_hpc', [foldertosave '/Freezing/']);
end


%% Find out how to divide in groups
A = mean(DatNormZ(ind_sort,:));
f2 = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.6]);
plot(mean(DatNormZ(ind_sort,60:100),2))
ylim([-0.6 0.6])
l1 = line(xlim, [mean(A) mean(A)], 'Color', 'b', 'LineWidth', 1);
l2 = line(xlim, [-0.1 -0.1], 'Color', 'r', 'LineWidth', 2);
line([100 100], ylim, 'Color', 'r', 'LineWidth', 2);
l3 = line(xlim, [0.08 0.08], 'Color', 'm', 'LineWidth', 2);
line([600 600], ylim, 'Color', 'm', 'LineWidth', 2);
xlabel('Neuron number')
ylabel('Firing at freezing zscored')
legend([l1 l2 l3], 'Mean', 'OFF-Freezing', 'ON-Freezing', 'Location', 'NorthWest');
makepretty

%% Groups mean
group_data = DatNormZ(ind_sort,:);
group_data = {group_data(1:100,:), group_data(101:599,:), group_data(600:end,:)};
titles = {'OFF on freezing', 'Not affected', 'ON on freezing'};
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for igroup = 1:length(group_data)
    subplot(3,1,igroup)
    plot(-2.05:0.05:6, mean(group_data{igroup}), 'Color', 'k');
    xlim([-2 6])
    ylim([-0.25 0.2])
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    ylabel('Firing zscored')
    xlabel('Time to freeze on (s)')
    title(titles{igroup});
    makepretty
end

if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f3,[foldertosave '/Freezing/SpikesOnFreezing_hpc_groups.fig']);
    saveFigure(f3, 'SpikesOnFreezing_hpc_groups', [foldertosave '/Freezing/']);
end

%% Try first 3 PCs 
[co, sc, ~,~,explained] = pca(DatNormZ');
f4 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for ipc = 1:3
    subplot(3,1,ipc)
    plot(-2.05:0.05:6, sc(:, ipc), 'Color', 'k');
    xlim([-2 6])
    ylim([-2 4.5])
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    ylabel(['PC#' num2str(ipc)])
    xlabel('Time to freeze on (s)')
    title(['PC#' num2str(ipc) ', ' num2str(round(explained(ipc),1)) '% of data explained']);
    makepretty
end

if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f4,[foldertosave '/Freezing/PCAOnFreezing_hpc_groups.fig']);
    saveFigure(f4, 'PCAOnFreezing_hpc_groups', [foldertosave '/Freezing/']);
end

%% To save
% DatNormZ, ind_sort, PlaceCells, SpInfo, FR, NeuronClass, Accelero