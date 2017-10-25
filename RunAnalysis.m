%% Semiautomated 3D object finder 
%
% *This scripts allows to analyze an image volume containing objects
% (i.e. synaptic punctate labeling) and segment each individual object.*
%
% The basic steps of this semi-automate approach involve:
%
% # Load image stacks of the volume to inspect and supplemental labelings
%
% <<LoadImages.png>>
%
% # A mask is optionally provided to limit the search volume
% # Candidate dots are detected automatically using an iterative 
%   thresholding method followed by watershed segmentation.
%
% <<IterativeThreshold.png>>
%
% # The user is requested to filter candidate objects according to one 
%   or more parameters, by setting a threshold (usually using ITmax)
% # The user refine the selection by 3D visual inspection using Imaris
%
% <<PunctaIdentification.png>>
%
% # Object distribution and density is calculated in the volume or across
%   the cell skeleton.
%
% The main search loop follows the logic described in:
%
%  Developmental patterning of glutamatergic synapses 
%  onto retinal ganglion cells.
%  Morgan JL, Schubert T, Wong RO. Neural Dev. (2008) 26;3:8.
%
% Originally written for the detection of postsynaptic PSD95 puncta on
% dendrites of gene-gun labelled retinal ganglion cells
%
% *Dependencies:*  
%
% * Imaris 7.2.3
% * Image Processing Toolbox
%

disp('---- ObjectFinder analysis ----');
% Get settings for the analysis
TPN = [pwd filesep]; % Instead reading TPN file, get current working folder
Settings.TPN = TPN;
save([TPN 'TPN.mat'], 'TPN');
save([TPN 'Settings.mat'], 'Settings');
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'find'])==0, mkdir([TPN 'find']); end %create directory to store steps

% Read images to be used for analysis
if exist([TPN 'Settings.mat'], 'file'), load([TPN 'Settings.mat']); end
tmpDir=[TPN 'I' filesep];
tmpFiles=dir([tmpDir '*.tif']);

% Get dimensions of the first image (assumes each channel stored as individual 3D .tif file coming from same acquisition)
tmpImInfo = imfinfo([tmpDir tmpFiles(1).name]);
zs = numel(tmpImInfo);
xs = tmpImInfo.Width;
ys = tmpImInfo.Height;

% Read all .tif images into a single Iraw matrix (X,Y,Z,Imge#)
Iraw=zeros(ys, xs, zs, numel(tmpFiles), 'uint8');
txtBar('Loading image stacks... ');
for i = 1:numel(tmpFiles)
    for j = 1:zs
        Iraw(:,:,j,i)=imread([tmpDir tmpFiles(i).name], j);
        txtBar( 100*(j+i*zs-zs)/(zs*numel(tmpFiles)) );
    end
end
txtBar('DONE');

Imax=squeeze(max(Iraw,[],3)); % Create a MIP of each image
cfigure(size(Imax,3)*10, 8); % Size panel according to # of images to display
for i = 1:size(Imax,3)    
    subplot(1,size(Imax,3),i)
    image(Imax(:,:,i)*(500/double(max(max(Imax(:,:,i))))))
    title(['Ch# ' num2str(i) ': ' tmpFiles(i).name]);
    set(gca,'box','off');
    set(gca,'YTickLabel',[],'XTickLabel',[]);
end
colormap gray(256);

% Ask user for image idendity settings
tmpPrompt = {'Neurites channel (0 = none):',...
             'Dots channel :',...
             'Mask channel (0 = use "D" mask):',...
             'xy resolution :',...
             'z resolution :',...
             '3D median filter post? (1=yes, 0=no)',...
             'Median filter kernel size: '};
tmpAns = inputdlg(tmpPrompt, 'Assign channels', 1, {'1', '3', '0',num2str(1/tmpImInfo(1).XResolution),'0.3','0','1'});

Settings.ImInfo.xNumVox = xs;
Settings.ImInfo.yNumVox = ys;
Settings.ImInfo.zNumVox = zs;
Settings.ImInfo.DenCh = str2double(tmpAns(1));
Settings.ImInfo.PostCh = str2double(tmpAns(2));
Settings.ImInfo.MaskCh = str2double(tmpAns(3));
Settings.ImInfo.xyum = str2double(tmpAns(4));
Settings.ImInfo.zum = str2double(tmpAns(5));
Settings.ImInfo.MedFilt = str2double(tmpAns(6));
Settings.ImInfo.MedFiltKern = str2double(tmpAns(7));

save([TPN 'Settings.mat'], 'Settings'); %fixed to save only Settings (9/2/09 HO)

% Write Channels into matlab files
if Settings.ImInfo.DenCh
    Dend = Iraw(:,:,:,Settings.ImInfo.DenCh);
    if Settings.ImInfo.MedFilt, Dend = medfilt3(Dend,[Settings.ImInfo.MedFiltKern,Settings.ImInfo.MedFiltKern,Settings.ImInfo.MedFiltKern]); end
    save([TPN 'Dend.mat'],'Dend');
end

if Settings.ImInfo.PostCh
    Post=Iraw(:,:,:,Settings.ImInfo.PostCh);
    if Settings.ImInfo.MedFilt, Post=medfilt3(Post,[Settings.ImInfo.MedFiltKern,Settings.ImInfo.MedFiltKern,Settings.ImInfo.MedFiltKern]); end
    save([TPN 'Post.mat'],'Post');
end

if Settings.ImInfo.MaskCh
    D = Iraw(:,:,:,Settings.ImInfo.MaskCh);
    D = D / max(max(max(D))); % Normalize mask max value to 1
    save([TPN 'D.mat'],'D');
    if isdir([TPN 'find'])==0, mkdir([TPN 'find']); end % Create directory to store steps
    saveastiff(Post, [TPN 'find' filesep 'PostMask.tif']); %save 3-D tiff image of the masked Post
elseif exist([TPN 'D.mat'], 'file')
    disp('Loading Mask from D.mat');
    load([TPN 'D.mat']);
else
    % Create dummy mask with all ones
    D = ones(size(Post), 'uint8');
    save([TPN 'D.mat'],'D');
end

close all;
clear i j Iraw Imax Is tmp* xs ys zs ans

% Ask user to input Experiment and ObjectFinder settings
tmpPrompt = {'Animal strain ', 'Age (Postnatal days)','Preparation (slice, whole-mount)',...
             'Labeling 1', 'Labeling 2', 'Labeling 3', 'Labeling 4',...
             'Orientation (NA,V,D,N,T)', 'Eccentricity (0=optic disc, 1=periphery)'};
tmpAns = inputdlg(tmpPrompt, 'Experiment data', 1, {'CD1', '45', 'Whole-mount','CtBP2','SyT2','','','NA', '0.5'});

Settings.VolumeInfo.Strain = tmpAns(1);
Settings.VolumeInfo.Age = str2double(tmpAns(2));
Settings.VolumeInfo.Prep = tmpAns(3);
Settings.VolumeInfo.Immuno1 = tmpAns(4);
Settings.VolumeInfo.Immuno2 = tmpAns(5);
Settings.VolumeInfo.Immuno3 = tmpAns(6);
Settings.VolumeInfo.Immuno4 = tmpAns(7);
Settings.VolumeInfo.Orientation = tmpAns(8);
Settings.VolumeInfo.Eccentricity = str2double(tmpAns(1));
save([TPN 'Settings.mat'], 'Settings')

% largest CtBP2 puncta, which is elipsoid of 1um diameter for xy and 
% 2um diamter for z, this will correspond to ~330 voxels with 
% 0.103um xy 0.3um z step voxel dimention. HO 1/25/2010

tmpPrompt = {'x-y pixel resolution (um) : ',...
             'z pixel resolution (um, default 0.3)',...
             'x-y diameter of the biggest dot (um, default 1)',...
             'z diameter of the biggest dot (um, default 2)',...
             'x-y diameter of the smallest dot (um, default 0.25)',...
             'z diameter of the smallest dot (um, normally 0.5)',...
             'Search block size (pixels, default 90) : ',...
             'Intensity thresholds stepping (default 2)',...
             'Multi-peak correction factor (default 0)',...
             'Minimum iteration threshold (default 2)'};
tmpAns = inputdlg(tmpPrompt, 'ObjectFinder settings', 1, {num2str(Settings.ImInfo.xyum),'0.3','1','2','0.25','0.5','90','2','0','2'});

MaxDotSize = (4/3*pi*(str2double(tmpAns(3))/2)*(str2double(tmpAns(3))/2)*(str2double(tmpAns(4))/2)) / (str2double(tmpAns(1))*str2double(tmpAns(1))*str2double(tmpAns(2)));
MinDotSize = (4/3*pi*(str2double(tmpAns(5))/2)*(str2double(tmpAns(5))/2)*(str2double(tmpAns(6))/2)) / (str2double(tmpAns(1))*str2double(tmpAns(1))*str2double(tmpAns(2)));
BlockBuffer = round(1.5/str2double(tmpAns(1)));

Settings.dotfinder.blockSize = str2double(tmpAns(7));
Settings.dotfinder.blockBuffer=BlockBuffer;  % Values between 15 to 50 didn't make much difference, 10 reduces number of dots found, use 15 for safety. HO 1/5/2010 Then switched from 15 to BlockBuffer to process images with different resolution. HO 6/8/2010
Settings.dotfinder.thresholdStep = str2double(tmpAns(8)); % step-size in the iterative search when looping through possible intensity values
Settings.dotfinder.maxDotSize = MaxDotSize; %max dot size exclusion criteria for single-peak dot DURING ITERATIVE THRESHOLDING, NOT FINAL.
Settings.dotfinder.minDotSize= 3; % min dot size exclusion criteria DURING ITERATIVE THRESHOLDING, NOT FINAL.
Settings.dotfinder.MultiPeakDotSizeCorrectionFactor = str2double(tmpAns(9)); % added by HO 2/8/2011, maxDotSize*MultiPeakDotSizeCorrectionFactor will be added for each additional dot joined to the previous dot, see dotfinder. With my PSD95CFP dots, super multipeak dots are rare, so put 0 for this factor.
Settings.dotfinder.itMin = str2double(tmpAns(10)); % added by HO 2/9/2011 minimum iterative threshold allowed to be analyzed as voxels belonging to any dot...filter to remove value '1' pass thresholds. value '2' is also mostly noise for PSD95 dots, so 3 is the good starting point HO 2/9/2011
Settings.dotfinder.peakCutoffLowerBound = 0.2; % set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
Settings.dotfinder.peakCutoffUpperBound = 0.2; % set threshold for all dots (0.2) after psychophysical testing with linescan and full 8-bit depth normalization HO 6/4/2010
Settings.dotfinder.minFinalDotITMax = 3; % minimum ITMax allowed as FINAL dots. Any found dot whose ITMax is below this threshold value is removed from the registration into Dots. 5 will be the best for PSD95. HO 1/5/2010
Settings.dotfinder.minFinalDotSize = MinDotSize; % Minimum dot size allowed for FINAL dots.

save([TPN 'Settings.mat'], 'Settings')
clear BlockBuffer MaxDotSize MinDotSize tmp*

% --- Find objects and calculate their properties ---

% Object finder (main loop)
Dots = findObjects(Post, D, Settings);
save([Settings.TPN 'Dots.mat'],'Dots');

% Create fields about sphericity of dots (Rounding)
Dots = fitSphere(Dots); 
save([TPN 'Dots.mat'],'Dots');

% Calculate distance to cell body position for each dot 
if exist([TPN 'Skel.mat'], 'file')
    Dots = distDotsToCellBody(TPN, Dots);
    save([TPN 'Dots.mat'],'Dots');
end

% Filter objects according to the following post-processing criteria (SG)
SGOptions.EdgeDotCut = 1;    % remove dots on edge of the expanded mask
SGOptions.SingleZDotCut = 1; % remove dots sitting on only one Z plane
SGOptions.xyStableDots = 0;
SGOptions.PCA = 0;
SGOptions.MinThreshold = 0;
filterObjects(TPN, SGOptions)

% Pass objects to Imaris, select ItMax as parameter and filter objects for it, 
% then export passing objects back to matlab values to use as threshold in SG.
exportObjectsToImaris(TPN)

% Grab passing sposts exported from Imaris
filterObjects(TPN)

% Group dots that are facing each other and recalculate properties
load([TPN 'Settings.mat']);
Grouped = groupFacingObjects(TPN, Dots, Settings);
Grouped = fitSphere(Grouped); 
save([TPN 'Grouped.mat'],'Grouped');

% Make finer Skeleton and calculate path distance of skels and dots from soma
if exist([TPN 'Skel.mat'], 'file')
    load([TPN 'Skel.mat'])
    load([TPN 'Settings.mat'])
    Skel = calcSkelPathLength(Skel);
    save([TPN 'Skel.mat'],'Skel')
    Skel = generateFinerSkel(Skel, Settings.ImInfo .xyum);
    Skel = calcSkelPathLength(Skel);
    save([TPN 'SkelFiner.mat'], 'Skel');
    Grouped = distDotsToCellBody(TPN, Grouped);
    Grouped = distDotsToSkel(Settings, Grouped, Skel);
    save([TPN 'Grouped.mat'],'Grouped');
end

% Generate Use.mat and calculate dendrite distribution
% TODO: Figure out what this function accumulates because used by all following
close all;
load([TPN 'SkelFiner.mat']);
anaMakeUseOnec(TPN, Dots, Grouped, Skel)
disp('Click on the maximum in the first plot (red) and hit return');
anaDotsDD(TPN);

% Generate heatmap plots (CA) of puncta over a MIP of the skeleton
close all;
[CA] = anaCAsampleUse(TPN);

% Plot dots distribution along dendrites length (Haru)
load([TPN 'SkelFiner.mat']);
calcPathLengthStats(TPN, Grouped, Skel);

disp('ObjectFinder analysis done!');
%% Change log
% _*Version 2.0*               created on 2017-10-24 by Luca Della Santina_
%
%  + Display more than 4 images if present in the I folder
%  + Automatically read x-y image resolution from tif files
%  + Added text progress bars to follow progress during processing steps
%  + Added debug=0/1 mode in subroutines to toggle text/graphic output
%  - Removed dependency from getVars() to get user input
%  - Removed redundances in user inputs (i.e. image resolution asked twice)
%  - Removed anaRa, anaRead, CAsampleCollect, StratNoCBmedian, Gradient
%  - Merged redundant scripts(anaCB/anaCBGrouped, anaRd/anaRdGrouped)
%  % Unified Imaris script under ObjectFinder_ names
%  % Imaris extensions provided visual confirmation dialog upon success
%  % Restructured main look into 4 distinct operations for maintainability
%  % Return Dots from objFinder and others instead of saving on disk
%  % Replaced GetMyDir dependency with pwd and work from current folder
%  % filterObjects(TPN) skips questions and just reprocess Imaris dots
%
% _*Version 1.5*   created on 2010-2011 by Hauhisa Okawa and Adam Bleckert_
%
%  % Improved speed of dot searching by resrtricting search within mask
%  % Improved speed of dot searching by working on Uint8 values
%  + Partial comments added to the original routines
%  + Puncta linear density, distances now calculated along dendritic path
%
% _*Version 1.0*                            created on 2008 by Josh Morgan_
%
% _*TODO*_
%
%  Find out why cannot re-run groupFacingObjects after analysis
%  Resolve minDotsize vs minFnalDotSize (current minDotSize fixed at 3px)
%  Investigate using multiple gmode to tune iterations for local contrast
%  Implement stepping=1 through gmode instead (current stepping of 2)
%  Remove dependency of HOgetVars from filterObjects.m
%  Check LDSCAsampleCollect and LDSCAsampleUse calculate density differently
%  When searching one signal (i.e. ctbp2 in IPL stack) prompt user several 
%  sample fields to choose an appropriate ITmax threshold
%
