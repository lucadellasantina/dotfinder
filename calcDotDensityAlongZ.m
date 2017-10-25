function calcDotDensityAlongZ(TPN, Dots, SG)
% Accumulate passing dots coordinates (xyz) into dPosPassF
passingIDs = SG.passI;
PassDotIDs = find(passingIDs==1);
dotDensity.dPosPassI = Dots.Pos(PassDotIDs,:); %(dotPassingID,:); % create directory of passing dots positions
dotDensity.ImSize = Dots.ImSize;

%dotDensity.zStart = 33; % use custom start position in the volume depth
%dotDensity.zEnd = 176;  % use custom start position in the volume depth
dotDensity.zStart = 1; % Use this to analyze entire volume
dotDensity.zEnd = dotDensity.dZsize; % Use this to analyze entire volume

dotDensity.binSize = (dotDensity.zEnd - dotDensity.zStart)/100; %1 Percent size of binning

tmpDensity=[];
for i=dotDensity.zStart:dotDensity.binSize:dotDensity.zEnd -1;
   tmpF = [];
   tmpF = find(dotDensity.dPosPassI(:,3) >= (i-dotDensity.binSize) & dotDensity.dPosPassI(:,3)<(i+dotDensity.binSize));
   tmpDensity(end+1) = numel(tmpF);
end

dotDensity.density = tmpDensity;

save([TPN 'dotDensity'], 'dotDensity') %fixed to save only Settings (9/2/09 HO)
clear tmp* i passingIDs PassDotIDs ans Dots SG;

%%
% Plot dot density distribution as a function of Volume depth.
tmpH = figure('Name', 'Puncta density');
set(tmpH, 'Position', [100 200 1200 500]);
set(gcf, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);
set(gcf, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 12);
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

hold on;
tmpY = dotDensity.density;
tmpX = [1:100];
plot(tmpY, 'k', 'MarkerSize', 8);

box off;
set(gca, 'color', 'none',  'TickDir','out');
ylabel('Number of puncta');
xlabel('Volume depth percentage');

clear tmp* plot_variance;
end
