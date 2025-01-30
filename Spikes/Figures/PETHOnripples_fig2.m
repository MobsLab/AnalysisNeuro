%% Parameters
IsSave=0;
mice = [797 798 828 861 882 905 906 911 912 977 994];
wi = 0.5; % in sec
binsize = 0.01; % in sec
nbins = 100;
titles = {'PreSleep', 'Conditioning', 'PostSleep'};

%% Get data
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

% Allocate arrays
b = cell(length(Dir.path),1); % behavior
r = cell(length(Dir.path),1); % ripples
s = cell(length(Dir.path),1); % spikes
ss = cell(length(Dir.path),1); % sleep scoring

for imouse = 1:length(Dir.path)
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch');
    r{imouse} = load([Dir.path{imouse}{1} '/Ripples.mat'], 'ripples');
    s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells');
    if strcmp(Dir.name{imouse},'Mouse861') || strcmp(Dir.name{imouse},'Mouse906') % bad scoring for 861 and no scoring for 906
        ss{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_Accelero.mat'], 'REMEpoch', 'SWSEpoch', 'Wake', 'Sleep');
    else
        ss{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_OBGamma.mat'], 'REMEpoch', 'SWSEpoch', 'Wake', 'Sleep');
    end
end

%% Allocate arrays
% Neurons
Neurons_SZ = cell(length(Dir.path), 1);
Neurons_PC_O = cell(length(Dir.path), 1);
Neurons_OC = cell(length(Dir.path), 1);
% Epochs
ripIS = cell(length(Dir.path), 1);
noripIS = cell(length(Dir.path), 1);
RipPreIS = cell(length(Dir.path), 1);
RipPostIS = cell(length(Dir.path), 1);
AllRipples = cell(length(Dir.path), 1);
NoRipples = cell(length(Dir.path), 1);

%% Get cells and epochs
for imouse = 1:length(Dir.path)
    
    % Sort neuros into place cells and other tentative pyr cells
    if ~isempty(s{imouse}.PlaceCells.idx)
        % Place cells shock zone
        cells_PC_SZ = s{imouse}.S(s{imouse}.PlaceCells.SZ);
        if length(cells_PC_SZ) > 0
            Q = MakeQfromS(cells_PC_SZ, binsize*1e4);
            tempdat = zscore(full(Data(Q)));
            for icell = 1:length(cells_PC_SZ)
                Neurons_SZ{imouse}{icell} = [Range(Q)/1e4 tempdat(:,icell)];
            end
        end
        
        % Place cells others
        cells_PC_O = s{imouse}.S(setdiff(s{imouse}.PlaceCells.idx, s{imouse}.PlaceCells.SZ));
        if length(cells_PC_O) > 0
            Q = MakeQfromS(cells_PC_O, binsize*1e4);
            tempdat = zscore(full(Data(Q)));
            for icell = 1:length(cells_PC_O)
                Neurons_PC_O{imouse}{icell} = [Range(Q)/1e4 tempdat(:,icell)];
            end
        end
        
        % Other pyramidal cells
        cells_OC = s{imouse}.S(intersect(setdiff(s{imouse}.BasicNeuronInfo.idx_SUA,s{imouse}.PlaceCells.idx),...
            find(s{imouse}.BasicNeuronInfo.neuroclass>0)));
        if length(cells_OC) > 0
            Q = MakeQfromS(cells_OC, binsize*1e4);
            tempdat = zscore(full(Data(Q)));
            for icell = 1:length(cells_OC)
                Neurons_OC{imouse}{icell} = [Range(Q)/1e4 tempdat(:,icell)];
            end
        end
        
    else
        Neurons_SZ{imouse} = [];
        Neurons_PC_O{imouse} = [];
        Neurons_OC{imouse} = [];
    end
    
    % Create the most important epochs
    % Ripples
    ripIS{imouse}=intervalSet(r{imouse}.ripples(:,1)*1E4, r{imouse}.ripples(:,3)*1E4);
    noripIS{imouse}=intervalSet((r{imouse}.ripples(:,1)+0.55)*1E4, (r{imouse}.ripples(:,3)+0.55)*1E4); % 550 ms forward
    % CondEpoch
    CondEpoch = or(b{imouse}.SessionEpoch.Cond1, b{imouse}.SessionEpoch.Cond2);
    CondEpoch = or(CondEpoch, b{imouse}.SessionEpoch.Cond3);
    CondEpoch = or(CondEpoch, b{imouse}.SessionEpoch.Cond4);
    
    %%% Create epochs
    % PreSleep
    SWSPreIS = and(ss{imouse}.SWSEpoch,b{imouse}.SessionEpoch.PreSleep);
    RipPreIS{imouse} = and(ripIS{imouse}, SWSPreIS);
    AllRipples{imouse}{1} = RipPreIS{imouse};
    NoRipples{imouse}{1} = and(noripIS{imouse}, SWSPreIS);
    % Cond
    RipCondIS = and(ripIS{imouse},CondEpoch);
    AllRipples{imouse}{2} = RipCondIS;
    NoRipples{imouse}{2} = and(noripIS{imouse}, CondEpoch);
    % PostSleep
    SWSPost = and(ss{imouse}.SWSEpoch,b{imouse}.SessionEpoch.PostSleep);
    RipPostIS{imouse} = and(ripIS{imouse},SWSPost);
    AllRipples{imouse}{3} = RipPostIS{imouse};
    NoRipples{imouse}{3} = and(noripIS{imouse}, SWSPost);
end

%% Plot PEST
% Allocate
SZ = cell(length(titles),1); % peth of all shock cells
OPC = cell(length(titles),1); % peth of other place cells
OCC = cell(length(titles),1); % peth of all pyramidal non-place cells
SZ_norip = cell(length(titles),1); % peth of all shock cells
OPC_norip = cell(length(titles),1); % peth of other place cells
OCC_norip = cell(length(titles),1); % peth of all pyramidal non-place cells

% Shock zone
cnt=0;
for imouse = 1:length(Dir.path)
    for num = 1:length(Neurons_SZ{imouse})
        cnt=cnt+1;
        for isess = 1:length(titles)
            % Ripples
            [r,i] = Sync(Neurons_SZ{imouse}{num}, Start(AllRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            SZ{isess}(cnt,:) = mean(T);
            % No ripples
            [r,i] = Sync(Neurons_SZ{imouse}{num}, Start(NoRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            SZ_norip{isess}(cnt,:) = mean(T);
        end
    end
end
% Other place cells
cnt=0;
for imouse = 1:length(Dir.path)
    for num = 1:length(Neurons_PC_O{imouse})
        cnt=cnt+1;
        for isess = 1:length(titles)
            % Ripples
            [r,i] = Sync(Neurons_PC_O{imouse}{num}, Start(AllRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OPC{isess}(cnt,:) = mean(T);
            % No ripples
            [r,i] = Sync(Neurons_PC_O{imouse}{num}, Start(NoRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OPC_norip{isess}(cnt,:) = mean(T);
        end
    end
end
% Other cells
cnt=0;
for imouse = 1:length(Dir.path)
    for num = 1:length(Neurons_OC{imouse})
        cnt=cnt+1;
        for isess = 1:length(titles)
            % Ripples
            [r,i] = Sync(Neurons_OC{imouse}{num}, Start(AllRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OCC{isess}(cnt,:) = mean(T);
            % No ripples
            [r,i] = Sync(Neurons_OC{imouse}{num}, Start(NoRipples{imouse}{isess})/1e4 ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OCC_norip{isess}(cnt,:) = mean(T);
        end
    end
end

%% Figure
numShock = size(SZ{1},1);
num_OPC = size(OPC{1},1);
num_OCC = size(OCC{1},1);
disp(['Number of shock zone place cells: ' num2str(numShock)]);
disp(['Number of other place cells: ' num2str(num_OPC)]);
disp(['Number of other neurons: ' num2str(num_OCC)]);

border_time = [-wi*1e3 wi*1e3 wi*1e3 -wi*1e3 -wi*1e3];
border_SZ = [0 0 numShock numShock 0];
border_OPC = [numShock+1 numShock+1 num_OPC+numShock num_OPC+numShock numShock+1];
border_OCC = [numShock+num_OPC+2 numShock+num_OPC+2....
    numShock+num_OPC+num_OCC numShock+num_OPC+num_OCC numShock+num_OPC+2];

% Sort conditioning {2}
SustVal = nanmean(SZ{2}(:,40:65),2);
[~,ind1] = sort(SustVal);
SZ_toplot{2} = SZ{2}(ind1,:);
SZ_norip_toplot{2} = SZ_norip{2}(ind1,:);
% Other place cells
SustVal = nanmean(OPC{2}(:,40:65),2);
[~,ind2] = sort(SustVal);
OPC_toplot{2} = OPC{2}(ind2,:);
OPC_norip_toplot{2} = OPC_norip{2}(ind2,:);
% Other cells
SustVal = nanmean(OCC{2}(:,40:65),2);
[~,ind3] = sort(SustVal);
OCC_toplot{2} = OCC{2}(ind3,:);
OCC_norip_toplot{2} = OCC_norip{2}(ind3,:);

% Sort others based on conditioninig
for isess = [1 3]
    % Shock
    SZ_toplot{isess} = SZ{isess}(ind1,:);
    SZ_norip_toplot{isess} = SZ_norip{isess}(ind1,:);
    % Other place cells
    OPC_toplot{isess} = OPC{isess}(ind2,:);
    OPC_norip_toplot{isess} = OPC_norip{isess}(ind2,:);
    % Other cells
    OCC_toplot{isess} = OCC{isess}(ind3,:);
    OCC_norip_toplot{isess} = OCC_norip{isess}(ind3,:);
end

% Figure
f1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
ax = arrayfun(@(i) subplot(2,3,i, 'NextPlot', 'add', 'Box', 'off'), 1:6);
for isess = 1:length(titles)
    axes(ax(isess));
    DatNormZ = [SZ_toplot{isess}; OPC_toplot{isess}; OCC_toplot{isess}];
    imagesc(-wi*1e3:binsize*1e3:wi*1e3,1:size(DatNormZ,1),DatNormZ);
    caxis([-0.3 0.3])
    axis ij
    xlim([-wi*1e3 wi*1e3])
    ylim([0 length(DatNormZ)])
    ylabel('Neuron #')
    xlabel('Time around a ripple (ms)')
    title(titles{isess})
    % Delineate groups
    line([0 0], ylim, 'Color', 'k')
    plot(border_time, border_SZ, 'Color', 'r');
    plot(border_time, border_OPC, 'Color', 'b');
    plot(border_time, border_OCC, 'Color', 'g');
    makepretty
    
    axes(ax(isess+length(titles)))
    DatNormZ = [SZ_norip_toplot{isess}; OPC_norip_toplot{isess}; OCC_norip_toplot{isess}];
    imagesc(-wi*1e3:binsize*1e3:wi*1e3,1:size(DatNormZ,1),DatNormZ);
    caxis([-0.3 0.3])
    axis ij
    xlim([-wi*1e3 wi*1e3])
    ylim([0 length(DatNormZ)])
    ylabel('Neuron #')
    xlabel('Time around a no-ripple timepoint (ms)')
    title(titles{isess})
    % Delineate groups
    line([0 0], ylim, 'Color', 'k')
    plot(border_time, border_SZ, 'Color', 'r');
    plot(border_time, border_OPC, 'Color', 'b');
    plot(border_time, border_OCC, 'Color', 'g');
    makepretty
end

%% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1, [foldertosave '/PETHOnRipples.fig']);
    saveFigure(f1, 'PETHOnRipples', foldertosave);
end