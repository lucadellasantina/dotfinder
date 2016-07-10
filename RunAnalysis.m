%% Get settings for the analysis

TPN = GetMyDir; %% Retrieves path for Folder that contains image folder (I)
save([TPN 'TPN.mat'], 'TPN');

if exist([TPN 'Settings.mat']); load([TPN 'Settings.mat']); end

%get cell info HO 6/25/2010
v.Animal = 'CD1';
v.CellName = 'Cell_';
v.CellType = 'OFFs';%'Type 1 BC';%'Flag N AC';
v.Age = 45;
v.Prep = 'Whole-mount GCL';
v.Bullet1 = 'CMV::CFP';
v.Bullet2 = 'CMV::PSD95-YFP';
v.Bullet3 = 'None';
v.Immuno1 = 'CtBP2Alexa568';
v.Immuno2 = 'Syt2Dylight649';
v.Immuno3 = 'None';
v.Immuno4 = 'None';
% v.RetinaUpperOrientation = 'NA';%270; %this will give standard to the cell orientation
% v.RetinaNasalOrientation = 'NA';%270;%'NA'; %this will give standard to the cell orientation
% v.CellOrientation = 'NA';%90;
% v.CellEccentricity = 0.3; %0 at the optic disc, 1 at the edge of the retina

v = LDSgetVars(v,'Enter cell information');
CellInfo = v;
save([TPN 'CellInfo.mat'], 'CellInfo')


%HO added puncta volume to voxel number conversion because images with
%different xyz resolution will be analyzed. 6/8/2010
Answer = inputdlg({sprintf(['xy resolution in um :\n']), sprintf(['z resolution in um :\n']), ...
  sprintf(['xy diameter of the max dot in um (normally 1) :\n']), ... %most largest CtBP2 puncta, which is elipsoid of 1um diameter for xy and 2um diamter for z, this will correspond to ~330 voxels with 0.103um xy 0.3um z step voxel dimention. HO 1/25/2010
  sprintf(['z diameter of the max dot in um (normally 2) :\n']), ...
  sprintf(['xy diameter of the min dot in um (normally 0.25) :\n']), ...
  sprintf(['z diameter of the min dot in um (normally 0.5) :\n'])}, ...
  'Volume to Voxel Number Conversion', 1, {'0.103','0.3','1','2','0.25','0.5'});
if isempty(Answer), return, end

MaxDotSize = (4/3*pi*(str2double(Answer(3))/2)*(str2double(Answer(3))/2)*(str2double(Answer(4))/2)) / (str2double(Answer(1))*str2double(Answer(1))*str2double(Answer(2)));
MinDotSize = (4/3*pi*(str2double(Answer(5))/2)*(str2double(Answer(5))/2)*(str2double(Answer(6))/2)) / (str2double(Answer(1))*str2double(Answer(1))*str2double(Answer(2)));
BlockBuffer = round(1.5/str2double(Answer(1)));

clear v
v.blockSize = 90;
v.blockBuffer=BlockBuffer; %changing this between 15 to 50 didn't make much difference, 10 could reduce the number of dots found, use 15 for safety. HO 1/5/2010 Then switched from 15 to BlockBuffer to process images with different resolution. HO 6/8/2010
v.thresholdStep = 2;
v.maxDotSize = MaxDotSize; %max dot size for single-peak dot DURING ITERATIVE THRESHOLDING, NOT FINAL, switched from 300 to MaxDotSize to process images with different resolution HO 6/8/2010
v.minDotSize=3; %min dot size DURING ITERATIVE THRESHOLDING, NOT FINAL.
v.MultiPeakDotSizeCorrectionFactor = 0; %added by HO 2/8/2011, maxDotSize*MultiPeakDotSizeCorrectionFactor will be added for each additional dot joined to the previous dot, see dotfinder. With my PSD95CFP dots, super multipeak dots are rare, so put 0 for this factor.
v.itMin = 2; % added by HO 2/9/2011 minimum iterative threshold allowed to be analyzed as voxels belonging to any dot...filter to remove value '1' pass thresholds. value '2' is also mostly noise for PSD95 dots, so 3 is the good starting point HO 2/9/2011
v.peakCutoffLowerBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
v.peakCutoffUpperBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
v.minFinalDotITMax = 3; % minimum ITMax allowed as FINAL dots. Any found dot whose ITMax is below this threshold value is removed from the registration into Dots. 5 will be the best for PSD95. HO 1/5/2010
v.minFinalDotSize = MinDotSize; %minimum dot size allowed as FINAL dots, switched from 3 to MinDotSize to process images with different resolution HO 6/8/2010

v = getVars(v,'Define Dotfinder Variables for Post dot');
Settings.dotfinder = v;
save([TPN 'Settings.mat'], 'Settings')

clear BlockBuffer MaxDotSize MinDotSize CellInfo v                   % Cleanings


%% Read images to be used for analysis

LDSanaRead(TPN);
%CellFillFlag = input('Do you have cell fill channel? Type 1 for yes, 0 for no.\n');
%
%if CellFillFlag;
%    
%     if ~exist('D')
%         disp('Select mask file');
%         D = LDStif2mat;
%         D = D / max(max(max(D))); % Normalize mask value, 0 /1 (useful for Fiji created mask that contains 0 / 255 values)
%         save([TPN 'D.mat'],'D');
%     end
%     disp('Select cell filling file');
%     Dend = LDStif2mat;
%     save([TPN 'Dend.mat'],'Dend');
%else;
    % Create dummy D for dotFinderInMask and PCA analysis
%    D = ones(size(Post), 'uint8');
%    save([TPN 'D.mat'],'D');
%end

%disp('Select dots image');
%Post = LDStif2mat;
%save([TPN 'Post.mat'],'Post');

%[ImInfo.yNumVox, ImInfo.xNumVox, ImInfo.zNumVox] = size(Post);
%ImInfo.DenCh = 0;
%ImInfo.PostCh = 1;
%ImInfo.ColoCh = 0;
%ImInfo.xyum = str2double(Answer(1));
%ImInfo.zum = str2double(Answer(2));

%if exist([TPN 'Settings.mat'])
%    load([TPN 'Settings.mat']);
%end
%Settings.ImInfo = ImInfo;
%save([TPN 'Settings.mat'], 'Settings') %fixed to save only Settings (9/2/09 HO)

clear Answer;
%% Dot finding process

LDSdotFinderInMaskWS(TPN) %in case you have a mask to restrict search into


%% make finer Skel, also calculate path distance of skels from soma 10/18/2011 HO
load([TPN 'Skel.mat'])
load([TPN 'Settings.mat'])
Skel = LDSSkelPathLengthCalculator(Skel);
save([TPN 'Skel.mat'],'Skel')
SkelFiner = LDSSkelFinerGenerator(Skel, Settings.ImInfo .xyum);
clear Skel;
Skel = SkelFiner;
save([TPN 'SkelFiner.mat'], 'Skel');
Skel = LDSSkelPathLengthCalculator(Skel);
save([TPN 'SkelFiner.mat'], 'Skel');
clear Skel SkelFiner Settings;


%% Post processing to compute additional information about dots

% Ratioing (create fields for ratio between the two channels)
LDSanaRa(TPN);

% Rounding (create fields about sphericity of dots)
LDSanaRd(TPN)

% Calculate position of cell body
load([TPN 'Settings.mat']);
LDSanaCB(TPN)
%% Post processing to remove error dots

% Running SG once, use only edgedot + singlez
LDSanaSGPCA2(TPN)

% In imaris import ItMax and ratio from SG.PassF and filter dots for ratio, find good
% values to use as threshold in SG.
open LDSanaDots2Imar 

% Run again first pass, now also threshold dots with the ratio value picked up from imaris
LDSanaSGPCA2(TPN)

% Now manually remove or add dots (Shify+click)
open LDSanaDots2Imar 

% Now do not run first pass (hit 0), this will only grab yes no from Imaris
LDSanaSGPCA2(TPN)
%%
% Group dots 
LDSanaGroup(TPN)
%HOanaGroupManual(TPN)

HOanaCBGrouped(TPN) %7/6/2010 HO
HOanaRdGrouped(TPN) %7/6/2010 HO

%register the path length of each dot from soma 10/18/2011 HO
%The function will also spit out the distance of dots to the closest Skel.
load([TPN 'Grouped.mat']);
load([TPN 'Settings.mat']);
load([TPN 'SkelFiner.mat']); %use finer version for accuracy
Grouped = LDSClosestSkelFinder(Settings, Grouped, Skel);
save([TPN 'Grouped.mat'], 'Grouped');
clear Skel Grouped;

%% GENERATE USE
% Luca: This function gets all the dots contained in Grouped.
% Make sure that it contains only the final selected dots

close all;
LDSanaMakeUseOnec(TPN)

LDSDotsDD(TPN); %input argument of yxum and zum were removed, instead load Settings to get those within the program.
%Although I don't like the way this program defines stratification-related
%parameters, just keep the format as it is.
%Can be modified in the future for bistratified GCs.

% Luca: select on the first plot (red) the peak and hit return

%% GENERATE CA FOR AREA
% Generate heatmap plots of puncta over a MIP of the skeleton
close all;
[CA] = LDSCAsampleUse(TPN);

%% DRAW OUTPUT
% Replot horizontally the maps
LDSCAsampleCollect(TPN);

%% GENERATE DEPTHDEV

%Not necessary to do this in each arbor for bistratified RGC.
%So the program was not modified for bistratified RGC. 10/15/2011 HO.
close all;
[DepthDev] = LDSStratNoCBmedian(TPN);

%% Plot dots distributions as concentric distance form soma (Josh)
% first plot is 2d distance, second plot is 3d distance according to haru
% This plots the average distribution of dots along the dendrites
% P/A = puncta over area (puncta/um^2)
% D/A = dendrite over area ? (um/um^2)
% P/D = linear density = puncta over dendrite (puncta/um)
[Grad] = LDSGradient(TPN);

%% Plot dots distribution along dendrites length (Haru)
LDSPathLengthStats(TPN);
