%% Plot and compare dot density for selected experiments
tmpXYres=0.098;
tmpZres=0.3;

tmpVars=who;
[tmpSel] = listdlg('PromptString','Choose Control', 'ListString', tmpVars);

tmpPuncta = [];
tmpDensity = [];
for i=1:numel(tmpSel)
    tmpDotDensity = evalin('base',char(tmpVars(tmpSel(i))));
    tmpPuncta(i, :) = tmpDotDensity.density(:);
    tmpDensity(i) = sum(tmpDotDensity.density)/(50.11*50.11*tmpZres*(tmpDotDensity.zEnd-tmpDotDensity.zStart));
end
tmpPunctaCtrl = tmpPuncta;
tmpDensityCtrl = tmpDensity;
tmpAvgPunctaCtrl = mean(tmpPuncta,1);

% Plot average with SEM.
tmpH = figure('Name', 'Puncta density');
set(tmpH, 'Position', [100 200 1200 500]);

set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontSize', 12);
set(gcf, 'DefaultAxesFontName', 'Arial')
set(gcf, 'DefaultAxesFontSize', 12)
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

%Plot control density
subplot('position', [0.06 0.12 0.7 0.8]);
hold on;
tmpY = tmpAvgPunctaCtrl;
tmpX = [1:100];
margin = std(tmpPuncta)/sqrt(size(tmpPuncta,1));
plot_variance(tmpX,tmpY-margin,tmpY+margin,[0.5 0.5 0.5]);
plot(tmpY, 'k', 'MarkerSize', 8);

box off;
set(gca, 'color', 'none',  'TickDir','out');
xlabel('Volume depth percentage');
ylabel('Number of puncta in a 50x50 µm^2 IPL stack');


tmpVars=who;
[tmpSel] = listdlg('PromptString','Choose 7 days Laser', 'ListString', tmpVars);

tmpPuncta = [];
tmpDensity = [];
for i=1:numel(tmpSel)
    tmpDotDensity = evalin('base',char(tmpVars(tmpSel(i))));
    tmpPuncta(i, :) = tmpDotDensity.density(:);
    tmpDensity(i) = sum(tmpDotDensity.density)/(50.11*50.11*tmpZres*(tmpDotDensity.zEnd-tmpDotDensity.zStart));
end
tmpPunctaTreat1 = tmpPuncta;
tmpDensityTreat1 = tmpDensity;
tmpAvgPunctaTreat1 = mean(tmpPuncta,1);

%Plot treated 1 density
hold on;
tmpY = tmpAvgPunctaTreat1;
tmpX = [1:100];
margin = std(tmpPuncta)/sqrt(size(tmpPuncta,1));
plot_variance(tmpX,tmpY-margin,tmpY+margin,[0.5 0 0]);
plot(tmpY, 'r', 'MarkerSize', 8);

tmpVars=who;
[tmpSel] = listdlg('PromptString','Choose 7 days Laser', 'ListString', tmpVars);

tmpPuncta = [];
tmpDensity = [];
for i=1:numel(tmpSel)
    tmpDotDensity = evalin('base',char(tmpVars(tmpSel(i))));
    tmpPuncta(i, :) = tmpDotDensity.density(:);
    tmpDensity(i) = sum(tmpDotDensity.density)/(50.11*50.11*tmpZres*(tmpDotDensity.zEnd-tmpDotDensity.zStart));
end
tmpPunctaTreat2 = tmpPuncta;
tmpDensityTreat2 = tmpDensity;
tmpAvgPunctaTreat2 = mean(tmpPuncta,1);

%Plot treated 2 density
hold on;
tmpY = tmpAvgPunctaTreat2;
tmpX = [1:100];
margin = std(tmpPuncta)/sqrt(size(tmpPuncta,1));
plot_variance(tmpX,tmpY-margin,tmpY+margin,[0 0.5 0]);
plot(tmpY, 'g', 'MarkerSize', 8);

% Average values and statistical comparisons
disp('-------------------------------------');
disp('Average density (dots/micron^3, SEM)   ');
disp('-------------------------------------');
disp(['Control: ' num2str(mean(tmpDensityCtrl)) setstr(177) num2str(mean(tmpDensityCtrl)/sqrt(numel(tmpDensityCtrl))) ' n=' num2str(numel(tmpDensityCtrl))]);
disp(['Treat1: ' num2str(mean(tmpDensityTreat1)) setstr(177) num2str(mean(tmpDensityTreat1)/sqrt(numel(tmpDensityTreat1))) ' n=' num2str(numel(tmpDensityTreat1))]);
disp(['Treat2: ' num2str(mean(tmpDensityTreat2)) setstr(177) num2str(mean(tmpDensityTreat2)/sqrt(numel(tmpDensityTreat2))) ' n=' num2str(numel(tmpDensityTreat2))]);
disp(' ');
disp('-------------------------------------');
disp('Statistical comparison (ranksum test)');
disp('-------------------------------------');
tmpPc1 = ranksum(tmpDensityTreat1, tmpDensityCtrl);
disp(sprintf('Avg density treat 1 vs ctrl.... p = %3.5f', tmpPc1));
tmpPc2 = ranksum(tmpDensityTreat2, tmpDensityCtrl);
disp(sprintf('Avg density treat 2 vs ctrl.... p = %3.5f', tmpPc1));
tmpPc12 = ranksum(tmpDensityTreat2, tmpDensityTreat1);
disp(sprintf('Avg density treat 2 vs treat 1.... p = %3.5f', tmpPc12));


%Plot average density
subplot('position', [0.8 0.12 0.18 0.8]);
hold on;
tmpH = bar([1:3], [mean(tmpDensityCtrl), mean(tmpDensityTreat1), mean(tmpDensityTreat2)], 'k');
set(tmpH,'FaceColor','none')

errorbar(1, mean(tmpDensityCtrl), std(tmpDensityCtrl)/sqrt(numel(tmpDensityCtrl)), 'k');
errorbar(2, mean(tmpDensityTreat1), std(tmpDensityTreat1)/sqrt(numel(tmpDensityTreat1)), 'r');
errorbar(3, mean(tmpDensityTreat2), std(tmpDensityTreat2)/sqrt(numel(tmpDensityTreat2)), 'g');

for i=1:numel(tmpDensityCtrl);
    plot(1,tmpDensityCtrl(i),'--ok','MarkerSize', 8);
end
for i=1:numel(tmpDensityTreat1);
    plot(2,tmpDensityTreat1(i),'--or','MarkerSize', 8);
end
for i=1:numel(tmpDensityTreat2);
    plot(3,tmpDensityTreat2(i),'--og','MarkerSize', 8);
end

box off;
title('Average density (puncta/µm^3)');
set(gca, 'color', 'none',  'TickDir','out');
xlim([0.5 3.5]);
ylim([0 1]);
set(gca,'XTick',[1 2 3]);
set(gca,'xtickLabel',{sprintf('Saline (%d)', numel(tmpDensityCtrl)),...
                      sprintf('Beads (%d)',   numel(tmpDensityTreat1)),...
                      sprintf('Beads (%d)',   numel(tmpDensityTreat2))});
text(1.75, 1.4, sprintf('p ctrl=\n%3.4f', tmpPc1));
text(2.75, 1.4, sprintf('p ctrl=\n%3.4f', tmpPc2));
text(2.75, 1.2, sprintf('p vs 1=\n%3.4f', tmpPc12));


clear tmp* i margin plot_variance ans;
end