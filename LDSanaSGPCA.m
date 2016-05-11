%% Comments and Error log (EL****)

% 041409 aab changed line 332 'if gP.reset, gP = deFault;' to 'if gP.reset, 
%   gP = default.gP;
% 042909 aab added readout for total puncta number and total passing puncta 
%   number line 504-505.
% 081809 aab  added Imaris YesNo optional input to define Yes and No puncta 
% 091009 had to turn off 'single' pass critereon (line 176)unique punctum 
%   finder. It was interfering with Imaris Dot surface integrations. Need to 
%   find what is wrong and fix.
% 092509 turned back on singles for normal test run
% 112009 -222 'added a filter for the TestCrit image so that noise in the
%   non arbor layers of Post is removed.' adam 
% 1/14/2010 HO added the option to remove dots facing outside the mask or
% the image.
% 1/15/2010 HO added 3D tif saving for TestCrit and SGCrit, plus Blank 3D
% tiff file to provide a background for yes and no when you do it in Fiji.
% Also, added forceCorrect for this 3D yes and no tiff files.
% 1/20/2010 HO added the skipping option for the first pass if you already
% did it and came back with yes no tif.
% 6/25/2010 HO All the lines with Dots.cut was commented out because we 
% don't do cut any more in anaMa.
% 7/9/2010 HO added 'if' loops in two places to skip loading Dend or DeltaF
% or DFOf if you don't have dendritic channel.
% 10/31/2010 Adam found a bug. I wrote 'if Settings.ImInfo.DendCh == 1' to
% enter enter/skip the 'if' loops that needs Dend channel, but DendCh is
% not always 1, so I changed them to 'if exist([TPN 'Dend.mat'])'.
% 1/4/2011 HO added 'if' loops before running EdgeDotCut, SingleZDotCut and
% PCA to skip them if you come back to re-do the analysis but not for all
% of them. (For example, just coming back to re-do PCA part.)
% 1/4/2011 HO modified the generation of TestCrit.tif to remove the
% background noise in Post channel by masking it with D.mat.
% 1/4/2011 HO added saving masked Post when YesNo3D was chosen so that yes
% no in 3D is easier if Post has lots of background.
% 1/29/2011 HO added forcing yes dots chosen in yes3Dforce.tif.
% 2/9/2011 HO removed Singles criterion because redundant voxels are
% now re-assigned to a signle dot at the end of dotfinder.
% 2/9/2011 HO moved the calculation of dP for all the dots inside if-end
% loop asking if PCA or MinThreshold will be run because there is no point
% of calculating it for every dot unless you use PCA or MinTreshold.
% 2/9/2011 HO removed the generation and the saving of Classify because
% most of them are duplicated in SG.
%
%AB 20101103 added simple definition of ratio to exclude dots outside
%dendrites
%
% 6/16/2011 HO added SGMinFieldOptions for choosing which criteria you want
% to use for minimum thresholding instead of setting non-chosen criteria to
% zero later because for some reason, some parameters (like deltaF) can be
% less than zero.

%%
function[] = LDSanaSGPCA(TPN)
%% Load variables and create directories

if ~exist('TPN'); TPN = GetMyDir; end
if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat'])
else
    Settings = [];
end

if ~exist([TPN 'find']), mkdir([TPN 'find']); end
if ~exist([TPN 'images']), mkdir([TPN 'images']); end

if exist([TPN 'find\SG.mat'])
    load([TPN 'find\SG.mat'])
else
    SG = [];
end

load([TPN 'Dots.mat'])


%% Skip the first pass if you come back to anaSG just for the latter smartGuide part 1/20/2010 HO
if exist([TPN 'find\SG.mat'])
    RedoFlag = input('First pass was already run. Redo? 1 = YES, 0 = No and go to smartGuide part.\n');
else
    RedoFlag = 1;
end

if RedoFlag == 1;
    %% Get user input, define variable Defaults
    
    SGOptions.EdgeDotCut = 1; %1 if you want to remove dots facing outside the mask or the image HO 1/14/2010
    SGOptions.SingleZDotCut = 1; %minimum number of voxels along z axis. Dots whose voxels are confined in a single z plane will be noise including speckling noise. This doesn't happen to PSD95 dots. HO 1/14/2010
    SGOptions.xyStableDots = 1; %1 if you want to remove moving dots
    SGOptions.PCA = 1; %1 if you want to select dots using PCA analysis
    SGOptions.MinThreshold = 0; %this is the conventional way that removes dots based on thresholds for each dot property defined in cP and dP HO 6/25/2010
    SGOptions.YesNo3D = 0; %1 if you want to generate 3D TestCrit and SGCrit for 3D yes no HO 1/14/2010
    if isfield(SG,'SGOptions')  % Check if SGOptions has been created
        if length(fieldnames(SGOptions)) == length(fieldnames(SG.SGOptions))
            SGOptions = SG.SGOptions; % Load SGOptions
        end
    end
    title = 'Which methods would you like to use? Enter 1 for yes, 0 for no.';
    SGOptions = LDSgetVars(SGOptions,title);  % Let user define first criteria
    SG.SGOptions = SGOptions;
    save([TPN 'find\SG.mat'],'SG')
    
    
    if (SGOptions.PCA == 1) || (SGOptions.MinThreshold == 1)
        % Collect dot properties %moved inside this if-end loop 2/9/2011 HO
        clear dP
        dP.iTMax = double(Dots.ITMax); % Minimum times a dot must pass iterative threshold
        dP.iTSum = Dots.ItSum; % Minimum sum of volume of each threshold pass
        dP.vol = Dots.Vol;  % Minimum volume to pass
        dP.meanBright = double(Dots.MeanBright);
        dP.compact = Dots.Round.Compact; % select for compact volumes
        dP.contrast=dP.iTMax./dP.meanBright;% Define Contrast as number of times passing threshold divided by mean brightness of puncta  or %max(1,(mB-(It*2)-gBG));
        dP.contrastVol=dP.iTMax./(dP.vol.^(1/3)); % Scale Contrast according to volume of puncta, or   %Contrast=It./max(mBG,1);
        for i = 1: Dots.Num
            %brights = Dots.Vox(i).Igm;
            brights = Dots.Vox(i).RawBright; %use RawBright instead of repeatedly-median-filtered Igm generated in anaRa HO 6/8/2010
            brights = sort(brights);
            LB = length(brights);
            numBrights = fix(LB/10)+1;
            topB = brights(LB-numBrights:LB); %top 10% + 1 voxels counted as max brightness
            bottomB = brights(1:numBrights); %bottom 10% counted as min brightness
            dP.internalContrast(i) = (mean(topB)-mean(bottomB));
        end
        % dP.dist2Dend = Dots.Dist2Dend; faces the wrong direction
        % dP.distToMask = Dots.DistToMask; faces the wrong direction
        if exist([TPN 'Dend.mat']) %in case you don't have dend channel, you can skip this part (DenCh will be 0). 7/9/2010 HO
            dP.ratio = Dots.Ratio;  %AB 20101103 added ratio for GRGM spots
            dP.deltaF = Dots.DF; % Delta F (fluorecence - predicted fluorescence)
            dP.deltaFOf = Dots.DFOf; % Delta F over F (predicted for red) for puncta
            dP.deltaFOfTop = Dots.DFOfTopHalf; % Minimum average delta F over f for the brightest %50 of voxels of a dot
            %The following calculation is not necessary, so commented out. HO
            %         gBG=Dots.Im.GreenBackGround;  % Get the Background value of the Green Channel
            %         mBG=dP.meanBright-gBG; % Get the Difference betteen dot brightness and background
            %         dist2CB = Dots.Dist2CB;
            %         if sum(isnan(dist2CB))
            %             dist2CB = Dots.Pos(:,3);
            %             'Using Z depth instead of CB',
            %         end
        end
        
        SG.dotProperties = dP;
        % Classify.dotProperties = dP; %remove Classify from the program. HO 2/9/2011
    end
    
    if SGOptions.PCA == 1;
        SGPCAOptions.iTMax = 1;
        SGPCAOptions.iTSum = 1;
        SGPCAOptions.vol = 1;
        SGPCAOptions.meanBright = 1;
        SGPCAOptions.compact = 1;
        SGPCAOptions.contrast = 1;
        SGPCAOptions.contrastVol = 1;
        SGPCAOptions.internalContrast = 1;
        SGPCAOptions.deltaF = 0; % Delta F (fluorecence - predicted fluorescence)
        SGPCAOptions.deltaFOf = 0; % Delta F over F (predicted for red) for puncta
        SGPCAOptions.deltaFOfTop = 0; % Minimum average delta F over f for the brightest %50 of voxels of a dot
        if isfield(Settings,'SGPCAOptions')  % Check if SGOptions has been created
            if length(fieldnames(SGPCAOptions)) == length(fieldnames(Settings.SGPCAOptions))
                SGPCAOptions = Settings.SGPCAOptions; % Load SGPCAOptions
            end
        end
        title = 'Which dot parameters would you like to use for PCA? Enter 1 for yes, 0 for no.';
        SGPCAOptions = LDSgetVars(SGPCAOptions,title);  % Let user define first criteria
        SG.SGPCAOptions = SGPCAOptions;
        save([TPN 'find\SG.mat'],'SG')
    end
    
    
    if SGOptions.MinThreshold == 1;
        SGMinThresholdOptions.iTMax = 0;
        SGMinThresholdOptions.iTSum = 0;
        SGMinThresholdOptions.vol = 0;
        SGMinThresholdOptions.meanBright = 0;
        SGMinThresholdOptions.compact = 0;
        SGMinThresholdOptions.contrast = 0;
        SGMinThresholdOptions.contrastVol = 0;
        SGMinThresholdOptions.internalContrast = 0;
        SGMinThresholdOptions.deltaF = 0; % Delta F (fluorecence - predicted fluorescence)
        SGMinThresholdOptions.deltaFOf = 0; % Delta F over F (predicted for red) for puncta
        SGMinThresholdOptions.deltaFOfTop = 0; % Minimum average delta F over f for the brightest %50 of voxels of a dot
        SGMinThresholdOptions.ratio = 1;
        if isfield(Settings,'SGMinThresholdOptions')  % Check if SGOptions has been created
            if length(fieldnames(SGMinThresholdOptions)) == length(fieldnames(Settings.SGMinThresholdOptions))
                SGMinThresholdOptions = Settings.SGMinThresholdOptions; % Load SGMinThresholdOptions
            end
        end
        title = 'Which dot parameters would you like to use for MinThreshold? Enter 1 for yes, 0 for no.';
        SGMinThresholdOptions = LDSgetVars(SGMinThresholdOptions,title);  % Let user define first criteria
        SG.SGMinThresholdOptions = SGMinThresholdOptions;
        save([TPN 'find\SG.mat'],'SG')
    end
    
    
    if SGOptions.MinThreshold == 1;
        %%Default dot criteria (first pass) Default values were inherited
        %%from Josh's original program, and iTSum was found to be most useful.
        %%Adam added ratio (ratio of dot channel over cell fill channel), which
        %%is also useful in some cases. HO 6/16/2011
        cP.iTMax = 0; % Minimum times a dot must pass iterative threshold. This used to be 3.
        cP.iTSum = 50; % Minimum sum of volume of each threshold pass. This used to be 50.
        cP.vol = 0;  % Minimum volume to pass. This used to be 3.
        cP.meanBright = 0; % Minimum mean brightness of puncta. This used to be 0.
        cP.compact = 0; % select for compact volumes. This used to be commented out.
        cP.contrast = 0; % measure of dot contrast (iterative threshold passes/ mean brightness). This used to be 0.05.
        cP.contrastVol = 0;  % Threshold for (ITMax)./(Dots.Vol.^(1/3)); This used to be 1.
        cP.internalContrast = 0; %This used to be 0.
        %cP.dist2Dend = 0; %distance of dot peak to dend filament
        %cP.distToMask = 1; %distance of dot peak to dend mask
        cP.deltaF = 0; % Delta F (fluorecence - predicted fluorescence). This used to be 2.
        cP.deltaFOf = 0; % Delta F over F (predicted for red) for puncta. This used to be 0.3.
        cP.deltaFOfTop = 0; % Minimum average delta F over f for the brightest %50 of voxels of a dot. This used to be 0.3.
        cP.ratio = 0.3; %AB 20101103 added ratio for GRGM spots, HO 2/16/2011 ~0.3 worked well for one of my PSD95CFP image.
        cP.reset = 0; % Resets Criteria to defaults
        
        defaultP.cP=cP; %Record as defaults
        
        if isfield(SG,'FirstThresholds')  % Check if FirstThresholds has been created
            if length(fieldnames(cP)) == length(fieldnames(SG.FirstThresholds))
                cP = SG.FirstThresholds; % Load FirstThresholds
            end
        end
        
        runFields = fieldnames(cP);
        for i = 1: length(runFields)
            if isfield(SGMinThresholdOptions,runFields{i}) % check if there is the same field in SGMinThresholdOptions for that property
                if SGMinThresholdOptions.(runFields{i}) ~= 1;
                    cP = rmfield(cP, runFields{i});
                end
            end
        end
        
        title = 'Enter minimum threshold for each parameter.';
        cP = LDSgetVars(cP,title);  % Let user define first criteria
        if cP.reset, cP = defaultP.cP; end % If reset was pressed, reset values to defalut
        
        SG.FirstThresholds = cP;
        Classify.UserThresholds = cP;
        save([TPN 'find\SG.mat'],'SG')
    end
    
    
    %% Applly Criteria to Puncta
    
    pass = ones(1,Dots.Num); %set up vector to record threshold passes
    
    %% Remove dots facing outside the mask or outside the image
    % -------------------------------------------------------------------------
    %  Created by Haruhisa Okawa on 2010-01-14
    % -------------------------------------------------------------------------
    
    if SGOptions.EdgeDotCut;
        if exist([TPN 'data\NonEdgeDots.mat']) %if-end loop added to remove the repetition 1/4/2011 HO
            EdgeDotCutRedoFlag = input('Dots facing outside already calculated. Redo? 1 = YES, 0 = Use pre-calculated.\n');
        else
            EdgeDotCutRedoFlag = 1;
        end
        
        if EdgeDotCutRedoFlag == 1; %if-end loop added to remove the repetition 1/4/2011 HO
            
            load([TPN 'D.mat']);
            D = uint8(D);
            D = bwperim(D, 6); %this operation will leave mask voxels facing 0 or outside the image as 1, and change the other mask voxels to 0.
            VoxIDMap = zeros(Dots.ImSize);
            for i=1:Dots.Num
                VoxIDMap(Dots.Vox(i).Ind)=i;
            end
            EdgeVoxIDMap = D.*VoxIDMap; %contour voxels located at the edge of the mask or image remains, and shows the dot ID#, other voxels are all 0.
            clear D;
            EdgeDotIDs = unique(EdgeVoxIDMap);
            if EdgeDotIDs(1) == 0;
                EdgeDotIDs(1)=[];
            end
            NonEdgeDots = ones(1,Dots.Num);
            NonEdgeDots(1, EdgeDotIDs)=0;
            save([TPN 'data\NonEdgeDots.mat'],'NonEdgeDots')
            pass = pass & NonEdgeDots; % Exclude edge dots
            
        else
            load([TPN 'data\NonEdgeDots.mat']);
            pass = pass & NonEdgeDots; % Exclude edge dots
        end
        disp(sprintf('Dots excluded because facing outside mask: %u', numel(NonEdgeDots) - numel(find(NonEdgeDots))));
    end
    
    %% Remove dots present in only one Z plane
    % -------------------------------------------------------------------------
    %  created by Haruhisa Okawa on 2010-01-14
    % -------------------------------------------------------------------------
    %
    % added the following line to exclude dots whose voxels spread only in one
    % z plane. This is not necessary for PSD95 dots but sometimes works well
    % with CtBP2 dots which include very dim noisy dots and speckling noise.
    %
    % Note: It is crucial to acquire image with Z-step smaller than dot size
    
    if SGOptions.SingleZDotCut;
        if exist([TPN 'data\NonSingleZDots.mat']) %if-end loop added to remove the repetition 1/4/2011 HO
            NonSingleZDotCutRedoFlag = input('Single Z plane dots already calculated. Redo? 1 = YES, 0 = Use pre-calculated.\n');
        else
            NonSingleZDotCutRedoFlag = 1;
        end
        
        if NonSingleZDotCutRedoFlag == 1; %if-end loop added to remove the repetition 1/4/2011 HO
            
            for i=1:Dots.Num
                zVoxNum(i) = length(unique(Dots.Vox(i).Pos(:,3)));
            end
            SingleZDotIDs = find(zVoxNum==1);
            NonSingleZDots = ones(1,Dots.Num);
            NonSingleZDots(1, SingleZDotIDs)=0;
            save([TPN 'data\NonSingleZDots.mat'],'NonSingleZDots')
            pass = pass & NonSingleZDots; % Exclude single Z dots
            
        else
            load([TPN 'data\NonSingleZDots.mat']);
            pass = pass & NonSingleZDots; % Exclude single Z dots
        end
        disp(sprintf('Dots excluded because present on a single Z plane: %u', numel(NonSingleZDots) - numel(find(NonSingleZDots))));
    end
    
    %% Remove Dots that are moving along Z stack
    % -------------------------------------------------------------------------
    %  Created by Luca Della Santina on 2011-10-06
    % -------------------------------------------------------------------------
    % Note: 3D median filter image before running dotfinder helps too
    
    if SGOptions.xyStableDots;
        if exist([TPN 'data\xyStableDots.mat']) %if-end loop added to remove the repetition 1/4/2011 HO
            NonSingleZDotCutRedoFlag = input('XY Stable Dots already computed. Redo? 1 = YES, 0 = Use pre-calculated.\n');
        else
            NonSingleZDotCutRedoFlag = 1;
        end
        
        if NonSingleZDotCutRedoFlag %   if-end loop added to remove the repetition 1/4/2011 HO
            xyStableDots = ones(1,Dots.Num);                                %Store Dots passing the test
            for i=1:Dots.Num
                zPlanes = sort(unique(Dots.Vox(i).Pos(:,3)));               % Find unique z planes and store value of their Z position
                xyPos = zeros(numel(zPlanes), 2);                           % Store here XY position of centroids for each Z plane
                
                for j=1:numel(zPlanes)                                      % For each z plane
                    zPosMask = Dots.Vox(i).Pos(:,3) == zPlanes(j);          % Find voxels belonging to the current Z plane
                    brightMask = Dots.Vox(i).RawBright == max(Dots.Vox(i).RawBright(zPosMask)); % Find brightest voxels
                    xyPos(j,:) = [mean(Dots.Vox(i).Pos(zPosMask & brightMask,1)),...
                        mean(Dots.Vox(i).Pos(zPosMask & brightMask,2))];    % Store position of centroid
                end
                
                if (std(xyPos(:,1))>2) | (std(xyPos(:,2))>2)                % Test centroids are moving for more than 2 st.dev on XY plane
                    xyStableDots(i)=0;                                      % False if not aligned
                end
            end
            save([TPN 'data\xyStableDots.mat'],'xyStableDots')
            pass = pass & xyStableDots;                                     % Trim moving dots away
        else
            load([TPN 'data\xyStableDots.mat']);
            pass = pass & xyStableDots;                                     % Trim moving dots away
        end
        disp(sprintf('Dots excluded because moving during aquisition: %u', numel(xyStableDots) - numel(find(xyStableDots))));
    end
    %% Plot PCA analysis and ask user to cluster good dots
    % -------------------------------------------------------------------------
    %  Created by Haruhisa Okawa on 2010-06-08
    % -------------------------------------------------------------------------
    
    if SGOptions.PCA;
        if exist([TPN 'data\PCAPassDots.mat']) %if-end loop added to remove the repetition 1/4/2011 HO
            PCARedoFlag = input('PCA was already run. Do you want to rerun PCA? Type 1 for rerunning, 0 for skipping.\n');
        else
            PCARedoFlag = 1;
        end
        
        if PCARedoFlag == 1; %if-end loop added to remove the repetition 1/4/2011 HO
            
            runFields = fieldnames(dP);
            dPMat=[];
            counter=0;
            for i = 1: length(runFields)
                if isfield(SGPCAOptions,runFields{i}) % check if there is the same field in SGPCAOptions for that property
                    if SGPCAOptions.(runFields{i}) == 1;
                        counter = counter+1;
                        dPMat(counter,:) = dP.(runFields{i})(pass)/std(dP.(runFields{i})(pass));
                        [runFields{i}  '  will be used for PCA.' ]
                    end
                end
            end
            dPMat = dPMat';
            
            [coef, zscore, pcvars] = princomp(dPMat);
            cumsum(pcvars./sum(pcvars) * 100) %check if two PCA components can account for >90% of variance, otherwise this analysis will not work good.
            coef %check each criteria is used well
            cfigure(20,20);
            scatter(zscore(:,1),zscore(:,2),'k.');
            hold on
            scatter(zscore(dP.iTSum(pass)<100,1),zscore(dP.iTSum(pass)<100,2),'r.');
            scatter(zscore(dP.iTSum(pass)<50,1),zscore(dP.iTSum(pass)<50,2),'b.');
            scatter(zscore(dP.iTSum(pass)<20,1),zscore(dP.iTSum(pass)<20,2),'m.');
            legend('100<=iTSum', '50<=iTSum<100', '20<=iTSum<50', 'iTSum<20');
            
            'press enter after circling PASS dots with a polygon'
            [X,Y]=getline(gcf,'closed');
            hpoly=plot(X,Y,'g');
            in=find(inpolygon(zscore(:,1),zscore(:,2),X,Y));
            scatter(zscore(in,1),zscore(in,2),'g');
            pause(3);
            clf
            
            TempPassDotNum = find(pass==1);
            PCAPassDotIDs = TempPassDotNum(in);
            PCAPassDots = zeros(1,Dots.Num);
            PCAPassDots(1, PCAPassDotIDs)=1;
            save([TPN 'data\PCAPassDots.mat'],'PCAPassDots')
            pass = pass & PCAPassDots; % Include PCA Pass dots
            
        else
            load([TPN 'data\PCAPassDots.mat']);
            pass = pass & PCAPassDots; % Include PCA Pass dots
        end
    end
    %% Save information about dots passing
  
    SG.pass1=pass';
    save([TPN 'find\SG.mat'],'SG')
    
    P=find(pass); %create list of passing puncta
    
    
else %if RedoFlag is not 1 (line46-51) 1/20/2010 HO
    load([TPN 'find\SG.mat']);
    pass=SG.pass1;
    P=find(pass);
end %if RedoFlag == 1 1/20/2010 HO

%% Draw image of the mask (maxRaw) and overlay dots on that

%%Create relevant image matricies
maxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
AllmaxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
maxPassed=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');
YXsize=Dots.ImSize(1)*Dots.ImSize(2); %?
DisAmAll=maxPassed; % Plot all Dot Voxels
DisAm=maxPassed;    % Plot voxels from Dots that pass criteria

%create pic of passed
for i = 1:Dots.Num
    %index is assigned from the top z plane to the bottom z plane in order
    %so for example the first assigned voxel in the second plane has an
    %index of YXsize+1. Thus, mod(Dots.Vox(i).Ind-1,YXsize)+1 will bring all
    %the voxels in the dot to the same z plane
    DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)=DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)+1;
    AllmaxID(mod(Dots.Vox(i).Ind-1,YXsize)+1)=i; %This will overwrite voxels with the same x y, but different z? 1/5/2010 HO
end

for i = 1: length(P)
    maxID(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=P(i);
    maxPassed(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200; % Assign number ( could be dot value)
    DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1;  %Stack labels
end

%%Assign different values for pixels corresponding to one vs more than one dot
OverLapped=DisAm;
OverLapped(DisAm==1)=150;
OverLapped(DisAm>1)=255;
%OverLapped at this point will be 150 in voxels where there was 1 passed dot voxel in
%that xy position throughout all z, 255 in the voxels where there were more
%than 1 passed dot voxels in that xy position throughout all z planes.
 

%load Raw
if exist([TPN 'images\maxRaw.mat'])
    load([TPN 'images\maxRaw.mat'])
elseif exist([TPN 'images\RawMax.tif']);
    maxRaw=imread([TPN 'images\RawMax.tif']);
    save([TPN 'images\maxRaw.mat'],'maxRaw');
else
    load([TPN 'Post.mat'])
    imgCutUp = max(Dots.Pos(SG.pass1,3))+3; %Dots.Pos(SG.pass1,3) will be the z positions of all the passed dots
    imgCutDown =  min(Dots.Pos(SG.pass1,3))-3;
    if imgCutDown < 1       %boarder gaurds
        imgCutDown = 1;           %bg
    end
    if imgCutUp > size(Post,3)  %bg
        imgCutUp = size(Post,3); %bg
    end
    
    %Mask the PostCut because imgCutUP doesn't help if the retina is
    %oblique 1/4/2011
    %PostCut = Post;
    load([TPN 'D.mat']);
    PostCut = Post.*D;
    
    PostCut(:,:,imgCutUp:end) = 0; % need to get rid of the noise in the cell body layers in the Post channel
    
    if exist([TPN 'Dend.mat']) %in case you don't have Dend channel 7/9/2010 HO
        load([TPN 'Dend.mat'])
        maxRaw=max(Dend,[],3);
    end
    maxRaw(:,:,2)=max(PostCut,[],3); %if Dend doesn't exist, this will pad 0 in the first z plane. 7/9/2010 HO
    save([TPN 'images\maxRaw.mat'],'maxRaw')
end

%%combine and save and image Comparison
TestID=uint16(maxRaw)*2^8;
TestCrit=maxRaw;
TestID(:,:,3)=maxID;
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
TestCrit(:,:,3)=OverLapped;
imwrite(TestCrit,[TPN 'find\TestCrit.tif'],'Compression','none')
%imwrite(TestID,[TPN 'images\TestID.tif'],'Compression','none')
subplot(1,2,1),image(TestCrit),pause(.01)


%% LUCA all of this seems shit, try to see what's good and what's not

%% Also save 3D tif file for the TestCrit 1/15/2010 HO
if SG.SGOptions.YesNo3D == 1;  
    VoxMap = zeros(Dots.ImSize);
    for i=1:length(P)
        VoxMap(Dots.Vox(P(i)).Ind)=1;
    end
    imwrite(VoxMap(:,:,1), [TPN 'find\TestCrit3D.tif'], 'tif', 'compression', 'none');
    for i=2:size(VoxMap,3)
        imwrite(VoxMap(:,:,i), [TPN 'find\TestCrit3D.tif'], 'compression', 'none', 'WriteMode', 'append');
    end

    %also save blank 3D tif file to provide a background for yes no when you do
    %it in Fiji.
    VoxMap = zeros(Dots.ImSize);
    imwrite(VoxMap(:,:,1), [TPN 'find\Blank3D.tif'], 'tif', 'compression', 'none');
    for i=2:size(VoxMap,3)
        imwrite(VoxMap(:,:,i), [TPN 'find\Blank3D.tif'], 'compression', 'none', 'WriteMode', 'append');
    end
    clear VoxMap;
end


%%  Check if the user information has been provided, terminate if not
if ~exist([TPN 'find\yes.tif']) || ~exist([TPN 'find\no.tif'])|| ~exist([TPN 'find\YesNo.mat'])
    'Please provide yes and no spots for further analysis'
end

%% Use user labeled images to Identify good and bad dots

Miss=[]; %Define vector to list artifacts
if exist([TPN 'find\YesNo.mat'])
    load ([TPN 'find\YesNo.mat'])
    Yunits = cell2mat(Settings.ImInfo(4));%convert cells to double
    NoIM = double(YesNo.No); 
    NoIM = fix(NoIM./Yunits); %convert from units (microns) to matrix values
    NoSize = size(AllmaxID);
    NoIMax = uint8(zeros(NoSize(1),NoSize(2)));
    NumNos = size(YesNo.No);
    for i = 1:NumNos(1);
        NoIMax(NoIM(i,1),NoIM(i,2)) = 1;
    end
    foundID = maxID(find(NoIMax>0));
    foundID = foundID(foundID>0);
    Miss = unique(foundID)';
else if exist([TPN 'find\no.tif'])
    NoIM=max(imread([TPN 'find\no.tif']),[],3);
    foundID = maxID(find(NoIM>0));
    foundID = foundID(foundID>0);
    Miss = unique(foundID)';
    end
end
        
Hit = []; %define vector to list identified puncta
if exist([TPN 'find\YesNo.mat'])
    %load ([TPN 'find\YesNo.mat'])
    %Yunits = cell2mat(Settings.ImInfo(4));%convert cells to double
    YesIM = double(YesNo.Yes); 
    YesIM = fix(YesIM./Yunits); %convert from units (microns) to matrix values
    YesSize = size(AllmaxID);
    YesIMax = uint8(zeros(YesSize(1),YesSize(2)));
    for i = 1:length(YesNo.Yes);
        YesIMax(YesIM(i,1),YesIM(i,2)) = 1;
    end
    FindYes = find(YesIMax>0);
    Hit = AllmaxID(FindYes);
    Hit = Hit(DisAmAll(FindYes)==1); 
    Hit = unique(Hit)';
else if exist([TPN 'find\yes.tif'])
    YesIM=max(imread([TPN 'find\yes.tif']),[],3);
    FindYes = find(YesIM>0);
    Hit = AllmaxID(FindYes);
    Hit = Hit(DisAmAll(FindYes)>0); %changed from FindYes == 0 to >0  there is a problem with overlapping dots in the z plane with NLG data
    Hit = unique(Hit)';
    Hit = Hit.*uint16(Hit>0); % exclude any zero values
    end
end

%If no hits are defined, define everything that is not a miss as a hit
if isempty(Hit)
    Hit = unique(AllmaxID);
    Hit = setdiff(Hit,Miss);
    Hit = Hit(Hit>0);
end

SG.manual.Miss=Miss;
SG.manual.Hit=Hit;

Classify.Miss = Miss;
Classify.Hit = Hit;



%run final pass adding manual Hits and Misses
passF=SG.pass1;
%SG.pass2=pass2'; %removed by luca
SG.passF=passF';


%% 3D version of forceCorrect HO 1/15/2010
Miss3D=[];
Hit3D=[];
if exist([TPN 'find\no3D.tif'])
    NoIM = zeros(Dots.ImSize);
    for i = 1:Dots.ImSize(3)
        NoIM(:,:,i)=imread([TPN 'find\no3D.tif'], i);
    end
    VoxIDMap = zeros(Dots.ImSize);
    for i=1:Dots.Num
        VoxIDMap(Dots.Vox(i).Ind)=i;
    end
    foundID = VoxIDMap(find(NoIM>0));
    foundID = foundID(foundID>0);
    Miss3D = unique(foundID);
    clear NoIM VoxIDMap;
end
if exist([TPN 'find\yes3D.tif'])
    YesIM = zeros(Dots.ImSize);
    for i = 1:Dots.ImSize(3)
        YesIM(:,:,i)=imread([TPN 'find\yes3D.tif'], i);
    end
    VoxIDMap = zeros(Dots.ImSize);
    for i=1:Dots.Num
        VoxIDMap(Dots.Vox(i).Ind)=i;
    end
    foundID = VoxIDMap(find(YesIM>0));
    foundID = foundID(foundID>0);
    Hit3D = unique(foundID);
    clear YesIM VoxIDMap;
end
if ~isempty(Miss3D) | ~isempty(Hit3D)
    if isfield(SG,'guideThresholds3D')  % Check if FirstThresholds has been created %switched to guideThresholds3D 1/15/2010 HO
        if length(fieldnames(gP)) == length(fieldnames(SG.guideThresholds3D)) %switched to guideThresholds3D 1/15/2010 HO
            gP = SG.guideThresholds3D; % Load FirstThresholds %switched to guideThresholds3D 1/15/2010 HO
        end
    end
    %%User defines which dot properties to use
    title = 'Which variables should change according to user labeling?';
    gP = getVars(gP,title);  % Let user define first criteria
    if gP.reset, gP = default.gP; end % If reset was pressed, reset values to defalut ****EL 041409
    if ~isempty(Miss3D)
        passF(Miss3D)=0;
    end
    if ~isempty(Hit3D)
        passF(Hit3D)=1;
    end
    %%Save Criteria
    SG.guideThresholds3D = gP; %switched to guideThresholds3D 1/15/2010 HO
    save([TPN 'find\SG.mat'],'SG')
    
    SG.passF=passF';
end

%% Save second order properties
save([TPN 'find\SG.mat'],'SG')

P=find(passF); %% list of passing puncta


%% Draw image
'Drawing images'
totalPassing = length(passF)
totalPassWithManual = length(P)        %042909 aab added readout for total passing puncta


maxSum=zeros(Dots.ImSize(1),Dots.ImSize(2));
maxID1=maxSum; maxID2 = maxSum; maxID3 = maxSum;


for i = 1: length(P)

    maxID1(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50; %rand outputs a random number between 0 and 1, so this line will produce a random number between 50 and 150 at the voxel positions of a passed dot
    maxID2(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;
    maxID3(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;
    maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1; %Isn't this the same as DisAm???

end

maxPassed2=uint8((maxSum>0)*200);

%OverLapped will show 60 in the voxels where there was 1 passed dot voxel in
%that xy position throughout all z, 100 in the voxels where there were more
%than 1 passed dot voxels in that xy position throughout all z planes, then
%add 150 more for all of these voxels. Why doing this way?
OverLapped=maxSum*0;
OverLapped(DisAm>0)=60;
OverLapped(DisAm>1)=100;
OverLapped(maxSum>0)=OverLapped(maxSum>0)+150;

MaxC=maxID1+(maxSum>1)*1000;
MaxC(:,:,2)=maxID2+(maxSum>1)*1000;
MaxC(:,:,3)=maxID3+(maxSum>1)*1000;
MaxC=uint8(MaxC);
subplot(1,2,2),image(MaxC),pause(.01)
% image(max(fullID,[],3)),P(i),pause

%combine and save and image Comparison
SGCrit=maxRaw;
SGCrit(:,:,1)=SGCrit(:,:,1);
SGCrit(:,:,2)=SGCrit(:,:,2);
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
SGCrit(:,:,3)=OverLapped;
imwrite(SGCrit,[TPN 'find\SGCrit.tif'],'Compression','none')
imwrite(MaxC,[TPN 'find\MaxCids.tif'],'Compression','none')

save([TPN 'Classify.mat'],'Classify')


%image(SGCrit), pause(.01)

%% Also save 3D tif file for the SGCrit 1/15/2010 HO
if SG.SGOptions.YesNo3D == 1;
    VoxMap = zeros(Dots.ImSize);
    for i=1:length(P)
        VoxMap(Dots.Vox(P(i)).Ind)=1;
    end
    imwrite(VoxMap(:,:,1), [TPN 'find\SGCrit3D.tif'], 'tif', 'compression', 'none');
    for i=2:size(VoxMap,3)
        imwrite(VoxMap(:,:,i), [TPN 'find\SGCrit3D.tif'], 'compression', 'none', 'WriteMode', 'append');
    end
end

