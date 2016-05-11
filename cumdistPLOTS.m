%% plot the cummulative distribution of the nearest neighbors
distNeighbor = [];
%%
try
    TPN; 
catch ME    
    TPN = GetMyDir;
end
load([TPN 'Grouped.mat'])

id = size(distNeighbor,1);
distNeighbor(id+1,:).NN = Grouped.NNDist;
distNeighbor(id+1,:).NoN = Grouped.NoNDist;
distNeighbor(id+1,:).Name = TPN;
clear Grouped; clear TPN;
%%

whitebg 'black' % though when saving black background is not saved
hold on
nnPSD = [];
nonPSD = [];
nnGRGM = [];

for i = 1:size(distNeighbor,1)/2
    [nnHand, psdNNstats] = cdfplot(distNeighbor(i*2-1).NN);
    [nonHand, psdNoNstats] = cdfplot(distNeighbor(i*2-1).NoN);
    set(nnHand, 'Color', 'm','LineWidth',.5,'LineStyle','-.')
    set(nonHand, 'Color', [0.7 0.7 0.7],'LineWidth',.5,'LineStyle','-.')
    nnPSD = cat(2,nnPSD, distNeighbor(i*2-1).NN);
    nonPSD = cat(2,nonPSD, distNeighbor(i*2-1).NoN);
    [nnHand, grgmNNstats] = cdfplot(distNeighbor(i*2).NN);
    set(nnHand, 'Color', 'y','LineWidth',.1,'LineStyle','-.')
    nnGRGM = cat(2,nnGRGM,distNeighbor(i*2).NN);
    psdNNstats
    psdNoNstats
    grgmNNstats
end

[nnHand, psdNNstats] = cdfplot(nnPSD);
[nonHand, psdNoNstats] = cdfplot(nonPSD);
set(nnHand, 'Color', 'm','LineWidth',3,'LineStyle','-')
set(nonHand, 'Color', [0.7 0.7 0.7],'LineWidth',3,'LineStyle','-')
[nnHand, grgmNNstats] = cdfplot(nnGRGM);
set(nnHand, 'Color', 'y','LineWidth',3,'LineStyle','-')

xlim([0 5])
grid off
box off
title 'Nearest Neighbor Cumulative Distribution'
xlabel 'microns'

%%
TPN = GetMyDir;
saveas(gcf,[TPN 'NNcdf'],'epsc')
save([TPN 'P12distNeighbor.mat'],'distNeighbor');