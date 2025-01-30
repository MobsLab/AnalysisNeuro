% 
% Detail:   For validation of the encoding. 
% 
% Need:     Used with the inferring.mat variable created with the decoder
%
% By Samuel Laventure - 2020-07
% 

clear all
load('nnBehavior.mat','behavior')
load('inferring.mat');
% treshall = [.05 .005 -.001 -.0095];
treshall = [0.005 0.001 -.001];
szside = 'right';
nopos = 1;  % if there are no position in the var pos use this (nopos = 1)/for complete session only
sesslen = times(end);
testst_id = find(times>=behavior.testEpochs(1),1,'first');

if nopos
    st = find(times>=behavior.trainEpochs(1),1,'first');
    posx = [zeros(1,st), interp1(1:length(behavior.positions(:,1)),behavior.positions(:,1)',linspace(1,length(behavior.positions(:,1)),length(inferring)-st))];
    posy = [zeros(1,st), interp1(1:length(behavior.positions(:,2)),behavior.positions(:,2)',linspace(1,length(behavior.positions(:,2)),length(inferring)-st))];
    
    x = posx/max([posx(1:testst_id-1) posy(1:testst_id-1)]);
    y = posy/max([posx(1:testst_id-1) posy(1:testst_id-1)]);
    pos = [x;y]';
end
%%
for itresh=1:length(treshall)
    tresh=treshall(itresh);
    infer_id = find(inferring(:,3)<tresh);
    decodedposX = inferring(infer_id,1);
    decodedposY = inferring(infer_id,2);
    posx = pos(infer_id,1);
    posy = pos(infer_id,2);
    predx = decodedposX;
    predy = decodedposY;
    % eleminate non-value
    zerx = find(decodedposX==0);
    zery = find(decodedposY==0);
    decodedposX(zerx)=[];
    decodedposY(zery)=[];
    posx(zerx)=[];
    posy(zery)=[];
    zerx = find(posx==0);
    zery = find(posy==0);
    decodedposX(zerx)=[];
    decodedposY(zery)=[];
    posx(zerx)=[];
    posy(zery)=[];

    maxt = length(posx);
    tps=1:length(posx);
    
    % Set the boundaries
    switch szside
        case 'right'
            minx=.8;
            maxx=1;
            miny=0;
            maxy=.55;
            ShockZone = [.78 .38 .26 .21];

        case 'left'
            minx=.35;
            maxx=.6;
            miny=.38;
            maxy=.58;
            ShockZone = [.36 .38 .26 .21];
    end
    maze = [0.36 0.38; 0.36 .96; 1.04 .96; 1.04 0.38; 0.78 0.38; 0.78 0.80; 0.62 0.80; 0.62 0.38; 0.36 0.38];
    
    % restrict to shock zone position
    for i=1:length(decodedposX)
        idx=find(decodedposX>minx & decodedposX<maxx);
        idy=find(decodedposY>miny & decodedposY<maxy);
        [val val_id] = intersect(idx,idy);
    end
    for i=1:size(inferring,1)
        idx=find(inferring(infer_id,1)>minx & inferring(infer_id,1)<maxx);
        idy=find(inferring(infer_id,2)>miny & inferring(infer_id,2)<maxy);
        [inZone val_id] = intersect(idx,idy);
    end
    
    % Figure: x and y linear trajectories real (red) + predicted (blue)
    if nopos   
        figure('Color',[1 1 1],'rend','painters','pos',[100+itresh*10 20+itresh*10 1600 800])
            subplot(2,2,1:2)
                hold on
                rectangle('Position',[0 minx length(x)  maxx-minx],'FaceColor',[.9 .9 .9])
                hold on
                plot(x,'r','LineWidth',1)
                hold on
                scatter(infer_id,inferring(infer_id,1),15,'filled')
                scatter(infer_id(inZone),inferring(infer_id(inZone),1),15,[0 .6 .2],'filled')
                xlim([1 length(x)])
                ylim([.25 1])
                title(['Position X (RED=real; GREEN=pred IN shk zone; BLUE=pred outside shk zone): threshold @' num2str(tresh)])
            subplot(2,2,3:4)
                hold on
                rectangle('Position',[0 miny length(x)  maxy-miny],'FaceColor',[.9 .9 .9])
                hold on
                plot(y,'r','LineWidth',1)
                hold on
                scatter(infer_id,inferring(infer_id,2),15,'filled')
                scatter(infer_id(inZone),inferring(infer_id(inZone),2),15,[0 .6 .2],'filled')
                xlim([1 length(x)])
                ylim([.25 1])
                title(['Position Y (RED=real; GREEN=pred IN shk zone; BLUE=pred outside shk zone): threshold @' num2str(tresh)])
    end 
    %%
end

    %% real vs prediction relation in space

    figure('Color',[1 1 1],'rend','painters','pos',[400+itresh*10 400+itresh*10 900 600])
        plot(pos(:,1),pos(:,2),'k','linewidth',1) 
        hold on
        plot(posx(1:maxt),posy(1:maxt),'ko','markerfacecolor','k')
        line([posx(1:maxt) decodedposX(1:maxt)]', [posy(1:maxt) decodedposY(1:maxt)]','color',[0.7 0.7 0.7]) 
        hold on
        plot(decodedposX(1:maxt), decodedposY(1:maxt),'k.')
        hold on
        scatter(posx,posy,15,tps,'filled')
        colorbar
        % plot(posx(1:maxt),posy(1:maxt),'k','linewidth',1) 
        hold on 
        plot(maze(:,1),maze(:,2),'k','LineWidth',2)
        rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
        title(['Thresh = ' num2str(tresh)])
        xlim([.3 1.2])
        ylim([.3 1.1])


    %% Location of stim vs real position
    figure('Color',[1 1 1],'rend','painters','pos',[100+itresh*10 20+itresh*10 1600 600])
            subplot(121)
                scatter(posx(val),posy(val),50,tps(val)','filled');
                hold on
                plot(maze(:,1),maze(:,2),'k','LineWidth',2)
                rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
                hold off
                xlim([0.335 1.045]);
                ylim([0.36 1])
                set(gca,'visible','off')
                set(gca, 'XTickLabel', {}, 'YTickLabel', {});
                title('Real position', 'FontSize', 16)
            subplot(122)
                scatter(decodedposX(val),decodedposY(val),50,tps(val)','filled');
                colorbar
                hold on
                plot(maze(:,1),maze(:,2),'k','LineWidth',2)
                rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
                hold off
                xlim([0.335 1.045])
                ylim([0.36 1])
                set(gca,'visible','off')
                set(gca, 'XTickLabel', {}, 'YTickLabel', {});
                title('Online inference', 'FontSize', 16)
                annotation('textbox', [0.4 0.87 0.1 0.1], 'String', ['Thresh = ' num2str(tresh)], 'FontWeight', 'bold', 'FontSize', 12,...
                    'EdgeColor', 'none')

end
            
%% color coded trajectories and prediction

%FigureRealVSInferred(pwd, tresh)
        
        
%% real vs decoded shock zone        


