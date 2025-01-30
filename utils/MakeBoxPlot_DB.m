function [boxhandle, pointshandle] = MakeBoxPlot_DB(A,Cols,X,Legends,ShowPoints, varargin)
% INPUT:
%
% - A                   data in a cell format 
%                       A = {Data1,Data2,Data3}
% - Cols (optional)     The color for each set of data
%                       Cols = {[1 0 0],[0 0 1],[0 1 0]}
%                       If you don't care about colors just put {} and everything wille be grey
% - X (optional)        The position to plot your data 
%                       X = [1,2,3]
%                       If you don't care just put []
% - Legends (optional)  The identity of your datasets for xlabels 
%                       Legends = {'bla' 'bla' 'bla'}
%                       If you don't care about colors just put {} 
% - ShowPoints          Shows data points.
%                       = 0 for no points; = 1 for points
% 
% Modified from MakeSpreadAndBoxPlot_SB.m

% Parameters
ConnectDots = false;

%% Optional parameters handling
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'connectdots'
            ConnectDots = varargin{i+1};
            if ~ShowPoints
                error('ConnectDots will work only with ShowPonts = true');
            else
                if ConnectDots ~= 0 && ConnectDots ~= 1 
                    error('Incorrect value for property ''IsSaveFig'' (type ''help PETH_freezing_hpc_fig'' for details).');
                end
            end
    end
end

%% Handling of arguments
if isempty(Cols)
    for i = 1:length(A)
        Cols{i} = [0.6 0.6 0.6];
    end
end

if isempty(X)
    for i = 1:length(A)
        X(i) = i*2;
    end
end

if isempty(Legends)
    for i = 1:length(A)
        Legends{i} = num2str(i);
    end
end

for i = 1:length(A)
    lenData(i) = length(A{i});
end
if ConnectDots
    if sum(diff(lenData)) > 0
        ConnectDots = false;
        warning('The number of points for each box plots is different. The dots will not be connected')
    end
    if sum(lenData > 20)
        ConnectDots = false;
        warning('The number of points for at least one box plot is higher than 20. The dots will not be connected')
    end
end


%% Plot
boxhandle = cell(length(A), 1);
if ShowPoints
    pointshandle = cell(length(A),1);
end
for k = 1 : length(A)
    if sum(isnan(A{k}))<length(A{k})
        boxhandle{k}=iosr.statistics.boxPlot(X(k),A{k}(:),'boxColor', Cols{k}, 'lineColor',[0.95 0.95 0.95],...
            'medianColor','k','boxWidth','auto', 'LineColor', 'k', 'LineWidth', 3, 'showOutliers',false);
        boxhandle{k}.handles.upperWhiskers.Visible='on';boxhandle{k}.handles.upperWhiskerTips.Visible='on';
        boxhandle{k}.handles.lowerWhiskers.Visible='on';boxhandle{k}.handles.lowerWhiskerTips.Visible='on';
        boxhandle{k}.handles.medianLines.LineWidth = 5;
        hold on
        if ShowPoints
            pointshandle{k}=plotSpread(A{k}(:),'distributionColors','k','xValues',X(k),'spreadWidth',0.8); hold on;
            set(pointshandle{k}{1},'MarkerSize',25)
            pointshandle{k}=plotSpread(A{k}(:),'distributionColors',Cols{k}*0.4,'xValues',X(k),'spreadWidth',0.8); hold on;
            set(pointshandle{k}{1},'MarkerSize',15)
            if ConnectDots
                temp = get(get(pointshandle{k}{3}, 'Children'), 'xdata');
                xdata{k} = temp{1};
                temp = get(get(pointshandle{k}{3}, 'Children'), 'ydata');
                ydata{k} = temp{1};
            end
        end
    end
    MinA(k) = min(A{k});
    MaxA(k) = max(A{k});
end


xlim([min(X)-1 max(X)+1])
rg = max(MaxA)-min(MinA);

if ShowPoints 
    ylim([min(MinA)-rg*0.1 max(MaxA)+rg*0.1])
end

if exist('Legends')
    set(gca,'XTick',X,'XTickLabel',Legends)
    box off
    set(gca,'FontSize',10,'Linewidth',1)
end

%% Connect dots
if ConnectDots
    for ibox = 1:length(xdata)-1
        for idot = 1:length(xdata{ibox})
            plot([xdata{ibox}(idot) xdata{ibox+1}(idot)], [ydata{ibox}(idot) ydata{ibox+1}(idot)], 'k', 'LineWidth', 2);
        end
    end
end
