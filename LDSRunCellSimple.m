% This RunCell simple version is for processing just puncta channel only,
% no cell fill, no mask, no colo channel.
% HOanaRead was replaced by just HOtif2mat and save the mat file as Post to
% take advantage of already written programs. HOdotfinderInMaskWS was
% replaced by HOdotfinderWS, which simply removed the masking of Post by D.
% 
% HO 3/22/2011

TPN = GetMyDir; %% Retrieves path for Folder that contains image folder (I)
save([TPN 'TPN.mat'], 'TPN');

if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end

%get cell info HO 6/25/2010
v.Animal = 'Thy1YFP';
v.CellName = 'Cell003';
v.CellType = 'G10';%'Type 1 BC';%'Flag N AC';
v.Age = 90;
v.Prep = 'Whole-mount GCL';
v.Bullet1 = 'None';
v.Bullet2 = 'None';
v.Bullet3 = 'None';
v.Immuno1 = 'GluR5alexa568';%'Syt2DyLight649';%'CtBP2Alexa514';
v.Immuno2 = 'TO-PRO3';%'RibeyeAdomainAlexa568';
v.Immuno3 = 'None';
v.Immuno4 = 'None';
% v.RetinaUpperOrientation = 'NA';%270; %this will give standard to the cell orientation
% v.RetinaNasalOrientation = 'NA';%270;%'NA'; %this will give standard to the cell orientation
% v.CellOrientation = 'NA';%90;
% v.CellEccentricity = 0.3; %0 at the optic disc, 1 at the edge of the retina

v = HOgetVars(v,'Enter cell information');
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
%v.percentBackground = 0.90; %this is not used HO
%v.punctaThreshold = 1; %this is not used HO
v.itMin = 2; % added by HO 2/9/2011 minimum iterative threshold allowed to be analyzed as voxels belonging to any dot...filter to remove value '1' pass thresholds. value '2' is also mostly noise for PSD95 dots, so 3 is the good starting point HO 2/9/2011
v.peakCutoffLowerBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
v.peakCutoffUpperBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
v.minFinalDotITMax = 3; % minimum ITMax allowed as FINAL dots. Any found dot whose ITMax is below this threshold value is removed from the registration into Dots. 5 will be the best for PSD95. HO 1/5/2010
v.minFinalDotSize = MinDotSize; %minimum dot size allowed as FINAL dots, switched from 3 to MinDotSize to process images with different resolution HO 6/8/2010
%v.roundThreshold = 50; %this is not used HO

v = getVars(v,'Define Dotfinder Variables for Post dot');
Settings.dotfinder = v;
save([TPN 'Settings.mat'], 'Settings')


%% Read image, generate dummy D if necessary
%IF YOU HAVE SECOND CYTOSOLIC FILL CHANNEL, RUN HOANAREAD JUST NORMALLY
CellFillFlag = input('Do you have cell fill channel? Type 1 for yes, 0 for no.\n');

if CellFillFlag == 1;
    HOanaRead(TPN);
    
elseif CellFillFlag == 0;
    %IF YOU DON'T HAVE SECOND CYTOSOLIC FILL CHANNEL, RUN THE FOLLOWING PART
    %load puncta image tif file and convert it to mat file and save.
    Post = HOtif2mat;
    save([TPN 'Post.mat'],'Post');
    D = ones(size(Post), 'uint8'); %create dummy D for dotFinderInMask and PCA analysis
    save([TPN 'D.mat'],'D');
    [ImInfo.yNumVox, ImInfo.xNumVox, ImInfo.zNumVox] = size(Post);
    ImInfo.DenCh = 0;
    ImInfo.PostCh = 1;
    ImInfo.ColoCh = 0;
    ImInfo.xyum = str2double(Answer(1));
    ImInfo.zum = str2double(Answer(2));

    if exist([TPN 'Settings.mat'])
        load([TPN 'Settings.mat']);
    end
    Settings.ImInfo = ImInfo;
    save([TPN 'Settings.mat'], 'Settings') %fixed to save only Settings (9/2/09 HO)

else
    'You need to enter 0 or 1...'
end

%%
%'Finding Dots'
LDSdotFinderInMaskWS(TPN) %in case you have a mask to restrict search into
%LDSdotFinderWS(TPN) %in case you don't have second cytosolic fill channel
%%
%no anaMa and anaCB, but do anaRa if you have second cytosolic fill channel
%'Ratioing'
if CellFillFlag == 1;
    HOanaRa(TPN);
end

% Rounding (rounds dots)
HOanaRd(TPN)

% Running SG once (refine the dot counting removing single Z plane dots and
% PCA analysis
HOanaSGPCA(TPN)

% Group dots 
HOanaGroup(TPN)

%HOanaGroupManual(TPN)

