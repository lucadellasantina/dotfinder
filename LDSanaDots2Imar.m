%% comments and errors
% 20101008 AB started function for adding ML passing and non passing spots
% to Imaris in order to see 3D 
%
%% 
% function[] = HOAnaDots2Imar(TPN)
% send dots found with dotfinder to Imaris
% imaris spots will be displayed as passing and nonpassing dots. 
% can also choose the statistics (dot properties) sent into imaris for both
% passing and non passing dots
%% get directory
if ~exist('TPN')
    TPN = GetMyDir;
end

%% connect to Imaris Com interface
try vImarisApplication.mVisible 
catch 
    vImarisApplication = actxserver('Imaris.Application');
    vImarisApplication.mVisible = true;
    pause(2) % imaris startup is slow sometimes,, and calling a file will cause a crash if imaris isnt running
end

%I changed the following part for my folder structure, but loading my imx
%file takes too long time, and it is not necessary to load filament, so
%load cell fill channel and MASKED PSD channel (PostMask under find folder)
%from Imaris, adjust the scale, then move onto the following part of the 
%code. HO 6/16/2011


% open correct file 
% if isempty(vImarisApplication.GetCurrentFileName)
% %     TPNfolder = TPN(1:end-1); %ab  remove '\' from TPN to get parent directory 
% %     [pathSTR,~,~] = fileparts(TPNfolder);%ab   get path to parent folder of matlab for Fiji folder access 
% %     TPNi = [pathSTR filesep 'ImarSkel' filesep]; %ab   filesep allows for cross platform accessability .. now TPNi accesses Fiji
%     TPNi = [TPN 'ImarSkel\'];
%     dTPNi = dir([TPNi '*.imx']);
%     vImarisApplication.FileOpen([TPNi dTPNi(1).name], ['LoadDataSet="eDataSetYes"']);
% end
'load the correct Imaris dataset'

%% Aquire Matlab Dot information  
if ~exist('TPN')
    TPN=GetMyDir;  % get directory of matlab dots
end
load([TPN 'Settings.mat']);
xyum = Settings.ImInfo.xyum; 
zum = Settings.ImInfo.zum;
load([TPN 'find' filesep 'SG.mat']);
load([TPN 'Dots.mat']);

% get the passing IDS
Iprompt = {'SG.passF:','SG.passI:'};
Idlg_title = 'which Grouping do you want:';
Inum_lines = 1;
if isfield(SG, 'passI')
    Idef = {'0','1'};
else
    Idef = {'1','0'};
end
Answer = inputdlg(Iprompt,Idlg_title,Inum_lines ,Idef);
if isempty(Answer), return, end

if str2double(cell2mat(Answer(1)))==1
    passingIDs = SG.passF';
else
    passingIDs = SG.passI';
end
%load([TPN 'Grouped.mat']);
%%
%HO 9/10/2011 changing passingIDs from equivalent to passF or passI (0 or 1
%in each element) to the dot IDs (1 to total number of dots)
PassDotIDs = find(passingIDs==1);
NoPassDotIDs = find(passingIDs==0);

dPosPassF = Dots.Pos(PassDotIDs,:); %(dotPassingID,:); % create directory of passing dots positions
dPosPassF(:,1:2)=(dPosPassF(:,1:2)-0.5)*xyum; %convert dots into actual values(um)
dPosPassF(:,3)=(dPosPassF(:,3)-0.5)*zum;
SPosPassF=[dPosPassF(:,2),dPosPassF(:,1),dPosPassF(:,3)]; % transpose x and y to convert from Matlab to Imaris
xyzVolConv = xyum^2*zum;
dVolpassF = Dots.Vol(PassDotIDs).*xyzVolConv;
dRadiusPassF = (dVolpassF.*3/(4*pi)).^(1/3);

%non passing dots  % 20101103 AB still going to use nonpassing dots from
%Dots as non passing dots (trying to encorporate grouped dots)
dPosNoPassF = Dots.Pos(NoPassDotIDs,:); %(dotPassingID,:); % create directory of passing dots positions
dPosNoPassF(:,1:2)=dPosNoPassF(:,1:2)*xyum; %convert dots into actual values(um)
dPosNoPassF(:,3)=dPosNoPassF(:,3)*zum;
SPosNoPassF=[dPosNoPassF(:,2),dPosNoPassF(:,1),dPosNoPassF(:,3)]; % transpose x and y to convert from Matlab to Imaris
xyzVolConv = xyum^2*zum;
dVolNopassF = Dots.Vol(NoPassDotIDs).*xyzVolConv;
dRadiusNoPassF = (dVolNopassF.*3/(4*pi)).^(1/3); 


%% create spots from matlab data
%passing dots
vSpotsAPosXYZ = SPosPassF;
vSpotsARadius = dRadiusPassF;
vSpotsAPosT = zeros(1,length(dPosPassF));
%non passing dots
vSpotsBPosXYZ = SPosNoPassF;
vSpotsBRadius = dRadiusNoPassF;
vSpotsBPosT = zeros(1,length(dPosNoPassF));


%% add new spots pass / non passing
% add passing spots
vSpotsA = vImarisApplication.mFactory.CreateSpots;
vSpotsA.Set(vSpotsAPosXYZ, vSpotsAPosT, vSpotsARadius);
vSpotsA.mName = sprintf('passing');
vSpotsA.SetColor(0.0, 1.0, 0.0, 0.0);
vImarisApplication.mSurpassScene.AddChild(vSpotsA);
%add non passing spots
vSpotsB = vImarisApplication.mFactory.CreateSpots;
vSpotsB.Set(vSpotsBPosXYZ, vSpotsBPosT, vSpotsBRadius);
vSpotsB.mName = sprintf('nonpassing');
vSpotsB.SetColor(0.0, 0.0, 1.0, 0.0);
vImarisApplication.mSurpassScene.AddChild(vSpotsB);
vSpotsB.mVisible = 0;

%% get spots statistics and create new from RunCell data
[aNames,aValues,aUnits,aFactors,aFactorNames,aIds]=vSpotsA.GetStatistics;
clear aNames; clear aValues; clear aUnits; clear aFactors; clear aIds;
%aFactorNames = {'Category';'Channel';'Collection';'Time'}; %necessary element for adding statistics

%choose the passing dots stats AB 20101105
%dotFN = fieldnames(Grouped);
dotFN = fieldnames(Dots);
[vDotsStats,vOk] = listdlg('ListString', dotFN,...
        'SelectionMode','multiple', ...
        'ListSize',[300 300], 'Name','DotsStats', ...
        'PromptString',{'Please select passing dot stats:'});
        if vOk<1, return, end
dotStatsNames = dotFN(vDotsStats);

for i = 1:length(dotStatsNames);
    for j = 1:length(vSpotsAPosXYZ)
        aNames{j,1} = strcat('RC_' ,dotStatsNames{i});
        %aValues(j,1) = single(Grouped.(dotStatsNames{i})(j));
        aValues(j,1) = single(Dots.(dotStatsNames{i})(PassDotIDs(j)));
        aUnits{j,1} = 'arb';
        aFactors{1,j} = 'Spots';
        aFactors{2,j} = '';
        aFactors{3,j} = '';
        aFactors{4,j} = '1';
        aIds(j,1) = int32(j-1);
    end
    vSpotsA.AddStatistics(aNames,aValues,aUnits,aFactors,aFactorNames,aIds);
    clear aNames; clear aValues; clear aUnits; clear aFactors; clear aIds;
end
%% non passing dots
% [aNames,aValues,aUnits,aFactors,aFactorNames,aIds]=vSpotsB.GetStatistics;
% clear aNames; clear aValues; clear aUnits; clear aFactors; clear aIds;
% %choose the dots stats AB 20101105
% dotBFN = fieldnames(Dots);
% [vDotsBStats,vOk] = listdlg('ListString', dotFN,...
%         'SelectionMode','multiple', ...
%         'ListSize',[300 300], 'Name','DotsStats', ...
%         'PromptString',{'Please select nonpassing dot stats:'});
%         if vOk<1, return, end
% dotStatsBNames = dotBFN(vDotsBStats);
% 
% 
% dotBID = single(find(passingIDs == 0));
% for i = 1:length(dotStatsNames);
%     for j = 1:length(vSpotsBPosXYZ);
%         aNames{j,1} = strcat('RC_' ,dotStatsBNames{i});
%         aValues(j,1) = single(Dots.(dotStatsBNames{i})(dotBID(j)));
%         aUnits{j,1} = 'arb';
%         aFactors{1,j} = 'Spots';
%         aFactors{2,j} = '';
%         aFactors{3,j} = '';
%         aFactors{4,j} = '1';
%         aIds(j,1) = int32(j-1);
%     end
%     vSpotsB.AddStatistics(aNames,aValues,aUnits,aFactors,aFactorNames,aIds);
%     clear aNames; clear aValues; clear aUnits; clear aFactors; clear aIds;
% end
'transfer of dots complete'
%% clear redundancies
clear vSpotsBPosT; clear vSpotsBPosXYZ; clear vSpotsBRadius; clear vSpotsA;
clear vSpotsAPosT; clear vSpotsAPosXYZ; clear vSpotsARadius; clear vSpotsB;
clear vDotsStats; clear passingIDs; 
clear dPosNoPassF; clear dPosPassF; clear dRadiusNoPassF; clear dRadiusPassF; clear dVolNopassF; clear dVolpassF; 
clear SPosNoPassF; clear SPosPassF;
clear Dots;
clear Grouped;
clear Settings;
clear SG;
clear aFactorNames;
clear dotStatsNames;

%% nothing 

