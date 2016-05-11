%% Run all Programs up to 1/19/09
%%Sets up sequence of programs to run to find puncta on labeled processes
%%Copyright Josh Morgan and Daniel Kerchensteiner 2005-2009

%% Comments and Error
% 
% 052709 Adams Version of Newer started 
% 091009 need to correct anaSG TCrit output of image to match Dots output to
%   imaris. Imaris Dots are correct, but TCrit is wrong.adam
% 110309 -74 'added v.itMaxMin to establish minimum iT threshold first
%   passing; see dotFinder'.adam
% 011010 HO added the saving of TPN, and the saving of the image name (the
% name of the cell etc) under Settings.ImageName.
% 021210 HO changed the cutOff criteria for estimating the contour of dots
% in dotfinding from just setting one cutOff parameter to setting both the 
% upper bound and the lower bound of cutOff. This will give more control
% and also be more intuitive. If upper bound is 0.4 for example, smallest
% ITMax dots (actually the smallest dot was assumed to be ITMax = 0 for
% simplicity although not existing) will take voxels whose iterations 
% counted at least >40% of max iteration of the dot, and cut the bottom 40%
% from the dot contour estimation. If the lower bound is set to 0.1,
% possible largest ITMax dots will keep top 90% of voxels as the contour of
% the dots. See dotFinderInMask for the change.
% 060810 HO since I have to process images with different xyz voxel size
% (for BCs in Grm6tdTom mice), I changed blockbuffer, maxDotSize,
% minDotSize and minFinalDotSize from simply giving some nubmer of voxels
% to the number of voxels converted from desired volume in um3.
% 062510 HO instead of saving the image name under Settings, I generated
% CellInfo.mat that holds all the information about the cell so that later
% in analysis you can use these informations to sort data from many cells.
% 020811 HO added itMin in addition to itMaxMin for dotfinder parameter.
% Now itMin sets minimum iteration threshold to be considered as non-noise
% voxels. Any voxel below this value is removed as potentially a part of
% dots before taken to watershed. Then, itMaxMin removes dots whose ITMax
% is lower than itMaxMin from the registration into Dots.
% 061711 HO added Imaris yes no, which Adam had created.

%%
TPN = GetMyDir; %% Retrieves path for Folder that contains image folder (I)
save([TPN 'TPN.mat'], 'TPN');

%1/5/2010 HO
%Before doing dotfinder, go to Imaris, generate Filament, save the filament
%information with Filament2ML, generate and save masks with FilamentMask.
%This will eliminate running anaImar and anaSk.

%generate tight mask for dendrite stratification analysis and monte carlo analysis (no extra xy and z radius)
load([TPN 'D.mat']);
load([TPN 'Mask.mat']);
save([TPN 'DDend.mat'],'D');
if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end
Settings.MaskDend = Mask;
save([TPN 'Settings.mat'], 'Settings');
clear Settings D Mask;

%then if you need to generate larger mask for colo (CtBP2 dots with respect
%to RGC dendrites and PSD95 dots), go back to Imaris and generate larger
%mask (2um extra xy radius plus 0.5um more extra z radius worked well for
%CtBP2 dots), and run the following lines.
load([TPN 'D.mat']);
load([TPN 'Mask.mat']);
save([TPN 'DColo.mat'],'D');
if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end
Settings.MaskColo = Mask;
save([TPN 'Settings.mat'], 'Settings');
clear Settings D Mask;

%finally go back to Imaris and generate mask for post dots (PSD95 dots,
%1.5um xy redius plus 0.5um more extra z radius worked well for PSD95
%dots). No need to change D, but need to put Mask under Settings.
load([TPN 'Mask.mat']);
if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end
Settings.Mask = Mask;
save([TPN 'Settings.mat'], 'Settings');
clear Settings Mask;
%now you can delete Mask.mat

%get cell info HO 6/25/2010
v.Animal = 'aPax6CreMG6YSDTA1';%'MG6YSTENTC5';%'WT C57BL6';%'MG6YSDTA4';
v.CellName = '110917c4';%'Josh022dC';
v.CellType = 'G1 or 2 or 6 or 10 ON GC';%'G10 RGC';
v.Age = 23;
v.Prep = 'whole-mount';
v.Bullet1 = 'CMVPSD95CFP';
v.Bullet2 = 'CMVtdTomato';
v.Bullet3 = 'None';
v.Immuno1 = 'Syt2DyLight649';%'CtBP2DyLight649';%'CtBP2Alexa514';
v.Immuno2 = 'None';%'RibeyeAdomainAlexa568';
v.Immuno3 = 'None';

v = LDSgetVars(v,'Enter cell information');
CellInfo = v;
save([TPN 'CellInfo.mat'], 'CellInfo')
    
%% Find out what has been done so far
DoAll.doAll = 0;
DoAll.read = 1;
DoAll.imaris = 0;
DoAll.skeleton = 0;
DoAll.dotfinder = 1;  % Set up new status varible
DoAll.ratio = 1;
DoAll.mask = 1;
DoAll.round = 1;
DoAll.smartGuide = 1;
DoAll.group = 1;
DoAll.dotfinderColo = 1;  % added HO 1/15/2010
DoAll.roundColo = 1;  % added HO 6/8/2010
DoAll.smartGuideColo = 1; % added HO 1/15/2010
DoAll.groupColo = 1; % added HO 1/15/2010

if exist([TPN 'Status.mat'])
    load([TPN 'Status.mat']) %Retreive previous progress
else
    Status = DoAll;
end

Status = LDSgetVars(Status, 'Which programs would you like to run?');
if Status.doAll
    Status = DoAll;
end

save([TPN 'Status.mat'],'Status')
Settings.Status = Status;

%% Set up necessary Directories
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps

%% Get More Variables

if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
end

if ~isfield(Settings,'ImInfo')
    'Need to collect image info at least once'
    Status.read = 1;
end

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

if Status.dotfinder

    clear v
    v.blockSize = 90;
    v.blockBuffer=BlockBuffer; %changing this between 15 to 50 didn't make much difference, 10 could reduce the number of dots found, use 15 for safety. HO 1/5/2010 Then switched from 15 to BlockBuffer to process images with different resolution. HO 6/8/2010
    v.thresholdStep = 2;
    v.maxDotSize = MaxDotSize; %max dot size for single-peak dot DURING ITERATIVE THRESHOLDING, NOT FINAL, switched from 300 to MaxDotSize to process images with different resolution HO 6/8/2010
    v.minDotSize=3; %min dot size DURING ITERATIVE THRESHOLDING, NOT FINAL.
    v.MultiPeakDotSizeCorrectionFactor = 0; %added by HO 2/8/2011, maxDotSize*MultiPeakDotSizeCorrectionFactor will be added for each additional dot joined to the previous dot, see dotfinder. With my PSD95CFP dots, super multipeak dots are rare, so put 0 for this factor.
    %v.percentBackground = 0.90; %this is not used HO
    %v.punctaThreshold = 1; %this is not used HO
    v.itMin = 3; % added by HO 2/9/2011 minimum iterative threshold allowed to be analyzed as voxels belonging to any dot...filter to remove value '1' pass thresholds. value '2' is also mostly noise for PSD95 dots, so 3 is the good starting point HO 2/9/2011
    v.peakCutoffLowerBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
    v.peakCutoffUpperBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
    v.minFinalDotITMax = 5; % minimum ITMax allowed as FINAL dots. Any found dot whose ITMax is below this threshold value is removed from the registration into Dots. 5 will be the best for PSD95. HO 1/5/2010
    v.minFinalDotSize = MinDotSize; %minimum dot size allowed as FINAL dots, switched from 3 to MinDotSize to process images with different resolution HO 6/8/2010
    %v.roundThreshold = 50; %this is not used HO

    v = HOgetVars(v,'Define Dotfinder Variables for Post dot');
    Settings.dotfinder = v;
    save([TPN 'Settings.mat'], 'Settings')
end

if Status.dotfinderColo

    clear v
    v.blockSize = 90;
    v.blockBuffer=BlockBuffer; %changing this between 15 to 50 didn't make much difference, 10 could reduce the number of dots found, use 15 for safety. HO 1/5/2010 Then switched from 15 to BlockBuffer to process images with different resolution. HO 6/8/2010
    v.thresholdStep = 2;
    v.maxDotSize = MaxDotSize; %max dot size for single-peak dot DURING ITERATIVE THRESHOLDING, NOT FINAL, switched from 300 to MaxDotSize to process images with different resolution HO 6/8/2010
    v.minDotSize=3; %min dot size DURING ITERATIVE THRESHOLDING, NOT FINAL.
    v.MultiPeakDotSizeCorrectionFactor = 0.25; %added by HO 2/8/2011, maxDotSize*MultiPeakDotSizeCorrectionFactor will be added for each additional dot joined to the previous dot, see dotfinder. With my standard CtBP2 dotfinding setting, 0.35-1 picked up too much noise, 0.25 worked well.
    %v.percentBackground = 0.90; %this is not used HO
    %v.punctaThreshold = 1; %this is not used HO
    v.itMin = 3; % added by HO 2/9/2011 minimum iterative threshold allowed to be analyzed as voxels belonging to any dot...filter to remove value '1' pass thresholds. value '2' is also most likely noise, so 3 is the good starting point (I initially though 2 would be good for CtBP2 DL649, but this picks up lots of noise especially in DTA line, changing this to 3 made it much better 10/1/2011)
    v.peakCutoffLowerBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
    v.peakCutoffUpperBound = 0.2; %changed to the set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
    v.minFinalDotITMax = 4; % minimum ITMax allowed as FINAL dots. Any found dot whose ITMax is below this threshold value is removed from the registration into Dots. For one of my CtBP2 dots tested, ITMax of 5 caused loss of some dots, but 4 did not lose dots compared to 3, dots that were not detected with 4 were also not detected with 3, so must be lost for other reasons. HO 1/5/2010
    v.minFinalDotSize = MinDotSize; %minimum dot size allowed as FINAL dots, switched from 3 to MinDotSize to process images with different resolution HO 6/8/2010
    %v.roundThreshold = 50; %this is not used HO

    v = HOgetVars(v,'Define Dotfinder Variables for Colo dot');
    Settings.dotfinderColo = v;
    save([TPN 'Settings.mat'], 'Settings')
end

% if Status.skeleton
% 
%     clear v
%     v.minObjSize = 50;
%     v.minFillSize = 10;
%     v.maxSegLength = 5;
% 
%     v=getVars(v , 'Define Skeletonization Variables');
%     Settings.skeleton = v;
%     save([TPN 'Settings.mat'], 'Settings')
% end

clear Settings

%% Check out files
if Status.read
    LDSanaRead(TPN)
    Status.read=0
    save([TPN 'Status.mat'],'Status')
end

%% Run Dot Processing
% if Status.imaris
%     'Loading Imaris Skeleton'
%     HOanaImar(TPN)
%     Status.imaris = 0;
%     save([TPN 'Status.mat'],'Status')
% end
% if Status.skeleton
%     'Finding Skel'
%     anaSk(TPN)
%     Status.skeleton=0;
%     save([TPN 'Status.mat'],'Status')
% end

%% Run dot finding

if Status.dotfinder
    'Finding Dots'
    LDSdotFinderInMaskWS(TPN)
    %anaDF(TPN)
    Status.dotfinder=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.ratio
    'Ratioing'
    LDSanaRa(TPN)
    Status.ratio=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.mask
    'Masking'
    LDSanaMa(TPN)
    LDSanaCB(TPN)
    %     anaFSc(TPN, DPN) %check for shifts
    Status.mask=0;
    save([TPN 'Status.mat'],'Status')
end
if Status.round
    'Rounding'
    LDSanaRd(TPN)
    Status.round=0;
    save([TPN 'Status.mat'],'Status')            
end

'Running SG once'
LDSanaSGPCA(TPN)

%Use spectrum and change intensity 0 to black for both Imaris and Fiji.
%LDSSGPassDotIDColorful3DMapGenerator(TPN, 'Post') %6/17/2011 HO

%open LDSanaDots2Imar %6/17/2011 HO

%Use Imaris to see where is the good minimum threshold for ITMax and ratio.
%Run SG again using that threshold, play with PCA.
%Then come back to Imaris to do yes and no. Use the colorful DotID 3-D map
%in Imaris to find potential noise component in multi-peak dots. Use Spots
%function, LDSiPassingSpots2MLDotsVer (not the original grouped version) to
%save passI under SG. Repeat SG, but skip redoing, it will create ImarCrit,
%so check the result in 2-D in photoshop. If everything is fine, proceed to
%manual grouping.

%open LDSanaGroupManual

if Status.group
    LDSanaGroup(TPN)
    LDSanaMakeUseOnec(TPN)
else
    LDSanaMakeUseOne(TPN)
end
%LDSPunctaCrit

LDSanaCBGrouped(TPN) %7/6/2010 HO
LDSanaRdGrouped(TPN) %7/6/2010 HO

%LDSanaGrouped2Imar(TPN) %6/16/2011 HO

%Adam's way, useful when you want to map the result.

%open LDSRunAnalysis


%the following part is added for Colo dots 1/15/2010 HO
if Status.dotfinderColo
    'Finding Dots'
    LDSdotFinderInMaskWSColo(TPN)
    %anaDF(TPN)
    Status.dotfinderColo=0;
    save([TPN 'Status.mat'],'Status')
end

%HOdotmodulatorColo(TPN)

if Status.roundColo
    'Rounding'
    LDSanaRdColo(TPN)
    Status.roundColo=0;
    save([TPN 'Status.mat'],'Status')            
end

'Running SG once for Colo dots'
LDSanaSGPCAColo(TPN)

%check the result
%Use spectrum and change intensity 0 to black for both Imaris and Fiji.
%HOSGPassDotIDColorful3DMapGenerator(TPN, 'Colo') %9/8/2011 HO
%HOSGPassDotIDGray3DMapGenerator(TPN, 'Colo', 50) %9/8/2011 HO

%9/13/2011 HO
%choose difficult dots (could be pass or no pass) to display
%HOSGPassDimDotIDGray3DMapGenerator(TPN, 'Colo') %9/8/2011 HO
%do yes no in Fiji for the above difficult dots, save them as Yes3DSample
%or No3DSample tif file in find folder, and run the following to find
%minimum thresholds, and go back to HOanaSGPCAColo.
%HOSGColoMinThresholdFinder(TPN)


%open HOanaDotsColo2Imar %9/8/2011 HO

if Status.groupColo
    LDSanaGroupManualColo(TPN)
    %LDSanaMakeUseOnecColo(TPN) %this will be not necessary for presynaptic colocalization puncta
else
    %LDSanaMakeUseOneColo(TPN) %this will be not necessary for presynaptic colocalization puncta
end

LDSanaRdGroupedColo(TPN) %7/28/2010 HO

%open LDSColocAnalysis

%% Comment Error log

% 050409 aab added 'Rounding' function to Newer 'RunCell' including addition
    % of anaRd.mat from Josh-08\AnalysisFunctions to Newer\Analysis Functions