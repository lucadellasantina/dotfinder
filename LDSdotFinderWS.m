%% comments and notes
% [date] -line   'notes'.user
%--------------------------------------------------------------------------
% HOdotFinderWS was created by copying HOdotFinderInMaskWS and just
% commenting out the masking of Post by D (line 79-83) so that the program
% can run without having cell fill channel or mask. HO 3/22/2011
%
%
% 110309 -173 'added minimum ITMax size critereon'.adam 
% 110309 -194 'redefined cutOff to attempt to capture more accurate punta
%   volume'.adam
%
% 010510 HO Gmode is set to the peak of the noise distribution rather than
%   minimum voxel intensity.
% 010510 HO Instead of throwing away dots found in the OUTER HALF of the
%   buffer region, throw away dots found in the entire buffer region.
% 010510 HO Modified to apply mask before dot-finding so that only spaces 
%   around dendrites are searched for potential dots.
% 011810 HO Added the saving of Gmode in Settings.dotfinder to check if
%   the newly defined Gmode is not too high.
% 021210 HO Changed the way to calculate cutOff to estimate the contour of
%   dots. Now you can assign the range of cutOff (upper bound for smallest 
%   ITMax dots (actually ITMax=0 although not existing) and lower bound for 
%   largest possible ITMax).
% 051210 -260 'added watershed function to assign more appropriate regions
%   to multipeak iterative threshold dots'. AAB
% 060710 HO modified AAB's watershed, do watershed for the whole
%   thresholdMap instead of repeating it for all the multipeak puncta. This
%   will remove the whole processings for multipeak puncta in the previous
%   version because multipeak puncta are divided into single-peak puncta in
%   wsTMLabels before going into the contour loop.
% 060710 HO also previous modification for Gmode was stupid, now just take
%   the mode of non-zero voxels, and do iteration up to Gmode+1.
% 020811 HO changed line 230, instead of removing all the voxels that have
% iteration of < Settings.dotfinder.itMaxMin, set a different threshold for
% removing voxels. Since itMaxMin was confusing name for this purpose, I 
% set up Settings.dotfinder.itMin, and now voxels whose iteration was < 
% itMin, which I set to 2 for now, will be removed before taken to 
% watershed because those voxels that have only 1 iteration will be likely
% noise. Then, HO introduced Settings.dotfinder.minFinalDotITMax in line 
% 290 so that dots whose itMax is < minFinalDotITMax will be excluded from
% the registration intoDots. This way, if minFinalDotITMax is set to 4, you
% still keep voxels that had 2 iterations for dots whose itMax was 4, but 
% not for those whose itMax was 3.
% 020811 HO added NumPeak to adjust max dot size limitation for multi-peak
% object (around line 190). I added if-end loop before dot registration to
% limit the size of the individual dots because allowing larger size
% limitation for multi-peak object creates artificially large dots.
% 020911 HO Some voxels could be shared by multiple dots, which would not 
% happen if the stack was not divided into sub-matrix. This happens around
% the dividing lines of sub-matrices. Instead of removing less brighter
% overlapping dot later in Singles criterion in anaSG, re-assign those
% voxels to a single dot. This operation was added in the last part of the
% program.
% 021511 HO added around 305, which excludes dots whose final size is
% < Settings.dotfinder.minFinalDotSize from the registration into Dots.

%%
function[]=LDSdotFinderWS(TPN) %#ok<INUSL>
% 
if ~exist('TPN')
    TPN = GetMyDir;
%large files will be broken into smaller blocks and then recombined
end

%% Get file names
%tic

load([TPN 'Settings.mat'])
load([TPN 'Post.mat']);
[Rys Rxs Rzs] = size(Post);
subplot(2,1,1),colormap gray(255);image(max(Post,[],3)),pause(.01)

%colormap gray(255) %standard grey colormap
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps

%% Apply mask before dot finding 1/5/2010 HO
% load([TPN 'D.mat']);
% D = uint8(D);
% Post = Post.*D;
% subplot(2,1,2),colormap gray(255);image(max(Post,[],3)),pause(1)

%% Enter Variables

v = Settings.dotfinder;
% 
% %%Image Variables
% % xyum=.103;
% % zum=.3;
% % % aspect=zum/xyum;% ratio of z to xy dimentions
% % channels=3;
% v.blockSize=200; %aproximate size of individual processing blocks
% v.blockBuffer=30; %amount of block to be cut off as edge buffer
% 
% %%Dot criteria
% step=2; %Sensitivity=grey value step of iterative threshold (2 was standard)
% MaxDot=7^3;  %Maximum Dot Volume (in pixels)= maximum dot size for iterative threshold (6^3 was standard)
% MinDot=5;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
% %MinDot=3;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
% % PunctaThreshold=1; %Minimum Number of Steps Passed = centroids less than puncta threshold are zeroed. 
% % EdgeOfPeak=.5; %Determine Dot Edge = ratio of edge brightness to peak brightness (.5 was standard)
% % minFilledVolume=3; %minimum number of pixels in final object contour (20 was standard)
% % RoundThreshold=50; %minimum roundness threshold (60 was standard)


%% Find Image Stats
%sort command during mode operation causes 'out of memory' error for my PC,
%so go back to the way not using mode command HO 3/25/2011
edges = -0.5:1:255.5;
N = histc(Post(:), edges);
Gmode = find(N == max(N)) - 1;

%Gmode = mode(double(Post(Post>0))); %060410 HO I made the line simpler
Gmax = double(max(Post(Post>0))); % 050510 AAB define the maximum values of post synaptic channel to scale ITMax by later (see line 221 etc...)


%% SubSample
if Rxs>v.blockSize
    NumBx=round(Rxs/v.blockSize);
    Bxc=fix(Rxs/NumBx); 
else
    Bxc=Rxs;
    NumBx=1;
end

if Rys>v.blockSize
    NumBy=round(Rys/v.blockSize);
    Byc=fix(Rys/NumBy); 
else
    Byc=Rys; 
    NumBy=1; 
end

if Rzs>v.blockSize
    NumBz=round(Rzs/v.blockSize);
    Bzc=fix(Rzs/NumBz); 
else
    Bzc=Rzs; 
    NumBz=1; 
end

%% Run Blocks
tic
%counters
countMultPeaks = 0;
maxNPeaks = 0;
%
for Bz=1:NumBz, for By=1:NumBy, for Bx=1:NumBx      %#ok<ALIGN>

PercentBlocksDone = ((Bz-1)*NumBy*NumBx+(By-1)*NumBx+Bx)/(NumBz*NumBx*NumBy)  %#ok<NASGU,NOPRT>

%Find real territory
Tystart=(By-1)*Byc+1;
Txstart=(Bx-1)*Bxc+1;
Tzstart=(Bz-1)*Bzc+1;
if By<Byc, Tyend=By*Byc; else Tyend=Rys; end
if Bx<Bxc, Txend=Bx*Bxc; else Txend=Rxs; end
if Bz<Bzc, Tzend=Bz*Bzc; else Tzend=Rzs; end

%Find buffered Borders (extend to image boarders for last blocks in row and column)
yStart=Tystart-v.blockBuffer;
yStart(yStart<1)=1;
yEnd=Tyend+v.blockBuffer;
yEnd(yEnd>Rys)=Rys;
xStart=Txstart-v.blockBuffer;
xStart(xStart<1)=1;
xEnd=Txend+v.blockBuffer;
xEnd(xEnd>Rxs)=Rxs;
zStart=Tzstart-v.blockBuffer;
zStart(zStart<1)=1;
zEnd=Tzend+v.blockBuffer;
zEnd(zEnd>Rzs)=Rzs;

%load([TPN 'Post.mat']); 
Igm = single(Post(yStart:yEnd,xStart:xEnd,zStart:zEnd)); 
%clear Post
[ys,xs,zs]=size(Igm);

%% FIND DOTS Green Channel%%
peakMap = zeros(ys,xs,zs,'uint8');  %set up matrix to map peaks
thresholdMap = zeros(ys,xs,zs,'uint8');   %set up matrix to sum passed thresholds
% indVect = 1:ys*xs*zs;
maxIntensity = uint16(max(Igm(:)));

%make threshold steps same between different blocks of image 1/5/2010 HO, a little modification on 6/4/2010 HO
if mod(maxIntensity, v.thresholdStep) ~= mod(Gmode+1, v.thresholdStep); %do iteration up to Gmode+1 for all the blocks, Gmode+1 is fine and the program doesn't need to iterate up to obviously noise level.
    maxIntensity = maxIntensity+1;
end

%for i = maxIntensity:-v.thresholdStep:Gmode+1 %do iteration up to Gmode+1 for all the blocks, Gmode+1 is fine and the program doesn't need to iterate up to obviously noise level.
%The above way to create loop didn't work with my PC. It just skipped the
%loop, so use the following way instead. 3/25/2011 HO
for iter = 1:(maxIntensity-Gmode-1)/v.thresholdStep;
    i = maxIntensity-(iter-1)*v.thresholdStep;
    %run thresholds through all relevant intensities
    clear Igl labels
    [Igl,labels] = bwlabeln(Igm>i,6);%label each area to check

    %reduce bitdepth if possible
    if labels<65536
        Igl=uint16(Igl); 
    end
    if labels <= 1
        labels =2;
    end
    nPixel = hist(Igl(Igl>0), 1:labels);
    %run all lables
    for p=1:labels
        pixelIndex = find(Igl==p);
        
        %HO added NumPeak to adjust max dot size limitation for multi-peak object 2/8/2011
        NumPeak = sum(peakMap(pixelIndex));
        if NumPeak == 0;
            NumPeak = 1;
        end

        %HO introduced v.MultiPeakDotSizeCorrectionFactor to adjust max dot size limitation for multi-peak object 2/8/2011
        if (nPixel(p) < v.maxDotSize+v.maxDotSize*(NumPeak-1)*v.MultiPeakDotSizeCorrectionFactor) && (nPixel(p) > v.minDotSize) % Morphology Filter, Puncta size criteria
            %identify peak in labeled object
            if sum(peakMap(pixelIndex))== 0 % sets up critireon to limit one peak per labeled field "Igl"
                peakValue = max(Igm(pixelIndex));
                peakIndex = find(Igl==p & Igm==peakValue);
                if numel(peakIndex) > 1
                    peakIndex = peakIndex(round(numel(peakIndex)/2));
                end
                [y,x,z] = ind2sub([ys xs zs], peakIndex);
                %Register in peak map
                peakMap(y,x,z) = 1;
            end
        else
            Igl(pixelIndex)=0;
        end
    end
    %%Add all passing labeled objects to thresholdMap
    thresholdMap(Igl>0)=thresholdMap(Igl>0)+1;
end 

clear Igl peakIndex


%% FIND DOT CONTOUR AND DIVIDE IF MULTIPLE PEAKS WITHIN

%thresholdMap(thresholdMap<v.itMaxMin) = 0; % 110309 Adam added minimum ITMax size critereon; only puncta with more than 1 IT pass are analyzed %I changed it to 2 because later in SG the min ITMax criterion is 3. HO 1/5/2010
thresholdMap(thresholdMap<v.itMin) = 0; %Since itMaxMin is confusing naming to me, I changed it to v.itMin, which I set to 2 for now because you want to save voxels of IT = 2 for the dot of ITmax = 4, for example. HO 2/8/2011
% [thresholdLabel nLabels] = bwlabeln(thresholdMap, 6); % find connectivity of objects in thresholdMap
% 
% if nLabels<65536 %reduce bit depth if possible
%     thresholdLabel = uint16(thresholdLabel);
% end


% watershed to find local minima and corresponding area AAB 2010_05_12, use uint8 for all the operations HO 6/4/2010
wsThresholdMapBin = uint8(thresholdMap>0);                                               % binary map of threshold
preKernel = ones(3,3,3);                                                            % create 26 connectivity kernel
wsThresholdMapBinOpen = imdilate(wsThresholdMapBin, preKernel);                     % dilate Bin map (dilated perimeter acts like ridges between background and ROIs
wsThresholdMapComp = imcomplement(thresholdMap);                                  % complement (invert) image. watershed fills holes not mountains. imcomplement creates complement using the entire range of the class, so for uint8, 0 becomes 255 and 255 becomes 0, but for double 0 becomes 1 and 255 becomes -254. 
%minTMPeaks = false(size(wsThresholdMap));                                           % create matrix of peaks 
%minTMPeaks(sub2ind(size(wsThresholdMap),yPeak(:), xPeak(:), zPeak(:))) = true;      % create matrix of peaks
%wsTMMod = imimposemin(wsThresholdMapComp, minTMPeaks | ~wsThresholdMapBinOpen);     % force local minima to be peaks and background outside of dilated region. imimposemin forces minimum to the masked pixels (0 for uint8, -inf for double), and remove local minima in the nonmasked pixels.
wsTMMod = wsThresholdMapComp.*wsThresholdMapBinOpen;     % this will force background outside of dilated region to 0, and leaves walls of 255 between puncta and background.
wsTMLabels = watershed(wsTMMod, 6);                                                 % 6 voxel connectivity watershed, this will fill background with 1, ridges with 0 and puncta with 2,3,4,... in double format
BackgroundLabel = mode(wsTMLabels(:));
wsTMLabels(wsTMLabels == BackgroundLabel) = 0; %seems that sometimes Background can get into puncta... so the next line was not good enough to remove all the background labels.
wsTMLabels = wsTMLabels.*double(wsThresholdMapBin); %masking out non-puncta voxels, this makes background and dilated voxels to 0. This also prevents trough voxels from being added back somehow with background. HO 6/4/2010
% add back the zero ridges in thresholdMap to their most similar neighbors
wsTMLZeros = find(wsTMLabels == 0 & thresholdMap > 0);                            % find zeros of watersheds inside of thresholdmap
if ~isempty(wsTMLZeros) % if there exist zeros in the map
    [wsTMLZerosY,wsTMLZerosX,wsTMLZerosZ] = ind2sub(size(thresholdMap),wsTMLZeros); %6/4/2010 HO
    clear nZeroID;
    for j = 1:length(wsTMLZeros)     
        %tempZeroMat = false(size(wsTMLabels));                                      % create temporary matrix for first zero position
        %tempZeroMat(wsTMLZeros(j))= true;                                           % define zero postision as true
        %tempZMID =  wsTMLabels(imdilate(tempZeroMat == 1, preKernel));              % create dilated matrix to examine neighbor connectivity
        %imdilate for whole block of image for all the zero ridges voxels
        %is too slow, try this instead HO 6/4/2010
        tempZMID =  wsTMLabels(max(1,wsTMLZerosY(j)-1):min(ys,wsTMLZerosY(j)+1), max(1,wsTMLZerosX(j)-1):min(xs,wsTMLZerosX(j)+1), max(1,wsTMLZerosZ(j)-1):min(zs,wsTMLZerosZ(j)+1)); %HO 6/4/2010
        nZeroID = mode(tempZMID(tempZMID~=0));                                      % find most common neighbor value (watershed) not including zero
        wsTMLabels(wsTMLZeros(j)) = nZeroID;                                        % re-define zero with new watershed ID (this process will act similar to watershed by making new neighboring voxels feed into the decision of subsequent zero voxels)
    end    
end
wsTMLabels = uint16(wsTMLabels);                                                    

wsLabelList = unique(wsTMLabels);
wsLabelList(1) = []; %remove background (now labeled as 0) from the list
nLabels = length(wsLabelList); 


for i = 1:nLabels
    peakIndex = find(wsTMLabels==wsLabelList(i) & peakMap>0); %this line adjustd for watershed HO 6/7/2010
    thresholdPeak = thresholdMap(peakIndex);
    nPeaks = numel(peakIndex);
    %following three lines commented out HO 6/7/2010
%     if nPeaks>maxNPeaks   %what is the maximum number of watersheds in a connected area
%         maxNPeaks = nPeaks;
%     end
    [yPeak xPeak zPeak] = ind2sub([ys xs zs], peakIndex);
    if nPeaks == 1
        if  (yPeak <= v.blockBuffer && yStart > 1) ||... %remove all dots found in buffer region instead of those found in the OUTER HALF of buffer region 1/5/2010 HO
                (xPeak <= v.blockBuffer && xStart > 1) ||...
                (zPeak <= v.blockBuffer && zStart > 1) ||...
                (ys-yPeak < v.blockBuffer && yEnd < Rys) ||...
                (xs-xPeak < v.blockBuffer && xEnd < Rxs) ||...
                (zs-zPeak < v.blockBuffer && zEnd < Rzs) ||...
                (thresholdPeak < v.minFinalDotITMax) %HO 2/8/2011 added excluding dots that did not reach minFinalDotITMax criterion
        else
            %cutOff = 0.5 * thresholdPeak; %original criterion
            %cutOff = (v.peakCutoff-double(thresholdPeak)/(255-Gmode))*thresholdPeak; %nonlinear criterion, changed from 255 to (255-Gmode) because iteration stops at Gmode HO 1/14/2010
            %changed to define the cutOff range by yourself (upper and lower
            %bound rather than just defining the upper bound previously) HO 2/12/2010
            PossibleMaxITMax = (Gmax - Gmode)/v.thresholdStep; %HO 2/12/2010 AAB 05/
            thresholdPeak = double(thresholdPeak); %HO 3/1/2010 if thresholdPeak is uint8, the next line always spit out zero, so thresholdPeak must be double.
            cutOff = (v.peakCutoffUpperBound - (v.peakCutoffUpperBound-v.peakCutoffLowerBound)*thresholdPeak/PossibleMaxITMax)*thresholdPeak; %HO 2/12/2010
            thresholdPeak = uint8(thresholdPeak); %HO 3/1/2010 convert thresholdPeak back to uint8.
            contourIndex = find(wsTMLabels==wsLabelList(i) & thresholdMap>=cutOff); %this line adjustd for watershed HO 6/7/2010
            
            if numel(contourIndex) < v.minFinalDotSize %HO 2/15/2011 added excludeing dots that did not reach minFinalDotSize criterion 
            
            else
                %2/8/2011 HO added this if-end loop to impose maxDotSize limit
                %to watershed-separated individual dots because I changed max
                %dot size limitation for multi-peak object during iterative
                %thresholding to maxDotSize*NumPeak (see line 193). This way,
                %you prevent individual dot from becoming too large.
                if numel(contourIndex) > v.maxDotSize
                    ITList = thresholdMap(contourIndex);
                    NewcutOff = thresholdPeak+1;
                    VoxSum = 0;
                    while VoxSum < v.maxDotSize
                        NewcutOff = NewcutOff-1;
                        VoxSum = length(find(ITList>=NewcutOff));
                    end
                    cutOff = NewcutOff + 1;
                    contourIndex = find(wsTMLabels==wsLabelList(i) & thresholdMap>=cutOff);
                end

                [yContour xContour zContour] = ind2sub([ys xs zs], contourIndex);

                if ~exist('Dots','var')
                    lastEntry = 0;
                else
                    [lastEntry dummy] = size(Dots.Pos);  %#ok<NASGU>
                end
                next = lastEntry+1;
                Dots.Pos(next,:) = [yPeak+yStart-1, xPeak+xStart-1,...
                    zPeak+zStart-1];
                Dots.Vox(next).Pos = [yContour+yStart-1,...
                    xContour+xStart-1, zContour+zStart-1];
                Dots.Vox(next).Ind = sub2ind([Rys Rxs Rzs],...
                    Dots.Vox(next).Pos(:,1), Dots.Vox(next).Pos(:,2),...
                    Dots.Vox(next).Pos(:,3));
                Dots.Vol(next) = numel(contourIndex);
                Dots.ITMax(next) = thresholdPeak;
                Dots.ItSum(next) = sum(thresholdMap(contourIndex));
                Dots.Vox(next).RawBright = Igm(contourIndex);
                Dots.Vox(next).IT = thresholdMap(contourIndex); %save also the iteration of each voxel 2/9/2011 HO
                Dots.MeanBright(next) = mean(Igm(contourIndex));
            end      
        end
                             
    else
        'more than 1 peaks were found!'
%         %% watershed to find local minima and corresponding area AAB 2010_05_12 
%         wsThresholdMap = zeros(size(thresholdMap)); %wsThresholdMap is double and thresholdMap is uint8
%         wsThresholdMap(find(thresholdLabel==i)) = thresholdMap(find(thresholdLabel==i));    % make threshold map for double peak dot in question
%         wsThresholdMapBin = wsThresholdMap>0;                                               % binary map of threshold
%         preKernel = ones(3,3,3);                                                            % create 26 connectivity kernel
%         wsThresholdMapBinOpen = imdilate(wsThresholdMapBin, preKernel);                     % dilate Bin map (dilated perimeter acts like ridges between background and ROIs
%         wsThresholdMapComp = imcomplement(wsThresholdMap);                                  % complement (invert) image. watershed fills holes not mountains. imcomplement creates complement using the entire range of the class, so for uint8, 0 becomes 255 and 255 becomes 0, but for double 0 becomes 1 and 255 becomes -254. 
%         minTMPeaks = false(size(wsThresholdMap));                                           % create matrix of peaks 
%         minTMPeaks(sub2ind(size(wsThresholdMap),yPeak(:), xPeak(:), zPeak(:))) = true;      % create matrix of peaks
%         wsTMMod = imimposemin(wsThresholdMapComp, minTMPeaks | ~wsThresholdMapBinOpen);     % force local minima to be peaks and background outside of dilated region. imimposemin forces minimum to the masked pixels (0 for uint8, -inf for double), and remove local minima in the nonmasked pixels.
%         wsTMLabels = watershed(wsTMMod, 6);                                                 % 6 voxel connectivity watershed
%         % add back the zero ridges in thresholdMap to their most similar neighbors
%         wsTMLZeros = find(wsTMLabels == 0 & wsThresholdMap > 0);                            % find zeros of watersheds inside of thresholdmap
%         if ~isempty(wsTMLZeros)                                                             % if there exist zeros in the map
%             for j = 1:length(wsTMLZeros)     
%                 tempZeroMat = false(size(wsTMLabels));                                      % create temporary matrix for first zero position
%                 tempZeroMat(wsTMLZeros(j))= true;                                           % define zero postision as true
%                 tempZMID =  wsTMLabels(imdilate(tempZeroMat == 1, preKernel));              % create dilated matrix to examine neighbor connectivity
%                 nZeroID = mode(tempZMID(tempZMID~=0));                                      % find most common neighbor value (watershed) not including zero
%                 wsTMLabels(wsTMLZeros(j)) = nZeroID;                                        % re-define zero with new watershed ID (this process will act similar to watershed by making new neighboring voxels feed into the decision of subsequent zero voxels)
%             end
%         end 
%         %
%         wsTMLabels = uint16(wsTMLabels);                                                    
%         wsPeakLabelsID = wsTMLabels(sub2ind(size(wsThresholdMap),yPeak(:), xPeak(:), zPeak(:))); % determine the watershed ID for each peak in question
%         if min(wsPeakLabelsID) < 1                                                               % error readout in case any peak ends up in a zero watershed (shouldnt happen)
%             testVoxelMaps = wsThresholdMap;
%             testDots.wsLabels(i).Vox = testVoxelMaps;
%             testDotID = sub2ind(size(wsThresholdMap),yPeak(:), xPeak(:), zPeak(:));
%             testDots(i).wsPos = testDotID;
%             save ([TPN 'testDots.mat'],'testDots')
%             'error - peaks in troughs'
%         end
%         for k = 1:nPeaks
%             PossibleMaxITMax = (Gmax - Gmode)/v.thresholdStep;                              %HO 2/12/2010
%             thresholdPeak = double(thresholdPeak);                                          %HO 3/1/2010 if thresholdPeak is uint8, the next line always spit out zero, so thresholdPeak must be double.
%             cutOff = (v.peakCutoffUpperBound - (v.peakCutoffUpperBound-v.peakCutoffLowerBound)*thresholdPeak(k)/PossibleMaxITMax)*thresholdPeak(k); %HO 2/12/2010
%             thresholdPeak = uint8(thresholdPeak);                                           %HO 3/1/2010 convert thresholdPeak back to uint8.
%             contourIndex = find(wsTMLabels==(wsPeakLabelsID(k)) & wsThresholdMap>=cutOff);
%             [yContour xContour zContour]= ind2sub([ys xs zs], contourIndex); 
%             peak(k).contour = [yContour xContour zContour];
%             peak(k).contourIndex = contourIndex;
%             
%             if  (yPeak(k) <= v.blockBuffer && yStart > 1) ||... %remove all dots found in buffer region instead of those found in the OUTER HALF of buffer region 1/5/2010 HO
%                     (xPeak(k) <= v.blockBuffer && xStart > 1) ||...
%                     (zPeak(k) <= v.blockBuffer && zStart > 1) ||...
%                     (ys-yPeak(k) < v.blockBuffer && yEnd < Rys) ||...
%                     (xs-xPeak(k) < v.blockBuffer && xEnd < Rxs) ||...
%                     (zs-zPeak(k) < v.blockBuffer && zEnd < Rzs)
%             else
%             
%                 if ~exist('Dots','var')
%                     lastEntry = 0;
%                 else
%                     [lastEntry dummy] = size(Dots.Pos); %#ok<NASGU>
%                 end
%                 next = lastEntry+1;
%                 Dots.Pos(next,:) = [yPeak(k)+yStart-1, xPeak(k)+xStart-1,...
%                     zPeak(k)+zStart-1];
%                 Dots.Vox(next).Pos = [peak(k).contour(:,1)+yStart-1,...
%                     peak(k).contour(:,2)+xStart-1,...
%                     peak(k).contour(:,3)+zStart-1];
%                 Dots.Vox(next).Ind = sub2ind([Rys Rxs Rzs],...
%                     Dots.Vox(next).Pos(:,1), Dots.Vox(next).Pos(:,2),...
%                     Dots.Vox(next).Pos(:,3));
%                 Dots.Vol(next) = numel(peak(k).contourIndex);
%                 Dots.ITMax(next) = thresholdPeak(k);
%                 Dots.ItSum(next) = sum(thresholdMap(peak(k).contourIndex));
%                 Dots.Vox(next).RawBright = Igm(peak(k).contourIndex);
%                 Dots.MeanBright(next) = mean(Igm(peak(k).contourIndex));
%                 countMultPeaks = countMultPeaks+1;                          % counter to see how many times an iterative region contains multiple peaks
%             end
%         end
    end  % if peaks
end  %All labels
        end
    end
end
% Dots.maxNPeaks = maxNPeaks; %commented out 6/7/2010 HO
% Dots.multPeaks = countMultPeaks; %commented out 6/7/2010 HO
Dots.ImSize = [Rys Rxs Rzs];
Dots.Num = size(Dots.Pos,1); %#ok<NASGU> 
save([TPN 'Dots.mat'],'Dots')
disp('iterative threshold done')

%Some voxels could be shared by multiple dots, which would not happen if
%the stack was not divided into sub-matrix. This happens around the
%dividing lines of sub-matrices. Instead of removing less brighter
%overlapping dot later in Singles criterion in anaSG, re-assign those voxels
%to a single dot. 2/9/2011 HO
disp('Re-assigning redundant voxels...')
VoxMap = uint8(zeros(Dots.ImSize));
VoxIDMap = zeros(Dots.ImSize);
[ys xs zs] = size(VoxMap);
TotalNumOverlapDots = 0;
TotalNumOverlapVoxs = 0;
for i = 1:Dots.Num
    OverlapVoxInd = find((VoxMap(Dots.Vox(i).Ind) > 0));
    if ~isempty(OverlapVoxInd)
        i
        TotalNumOverlapDots = TotalNumOverlapDots+1;
        TotalNumOverlapVoxs = TotalNumOverlapVoxs+length(OverlapVoxInd);
        OverlapVoxInds = Dots.Vox(i).Ind(OverlapVoxInd);
        OverlapVoxIDs = VoxIDMap(OverlapVoxInds);
        VoxMap(Dots.Vox(i).Ind) = 1;
        VoxIDMap(Dots.Vox(i).Ind) = i;
        VoxIDMap(Dots.Vox(i).Ind(OverlapVoxInd)) = 0; %bring the IDs of the shared voxels back to zero
        %Re-assign each overlapping voxel with either ID fighting over the voxel
        %using the same way as filling dot edges of watershed, HO 2/9/2011
        [OverlapVoxY, OverlapVoxX, OverlapVoxZ] = ind2sub(size(VoxMap), OverlapVoxInds);
        for k = 1:length(OverlapVoxInds)
            SurroudingIDs =  VoxIDMap(max(1,OverlapVoxY(k)-1):min(ys,OverlapVoxY(k)+1), max(1,OverlapVoxX(k)-1):min(xs,OverlapVoxX(k)+1), max(1,OverlapVoxZ(k)-1):min(zs,OverlapVoxZ(k)+1));
            WinningID = mode(SurroudingIDs((SurroudingIDs==i) | (SurroudingIDs==OverlapVoxIDs(k)))); % find which of two dots fighting over the voxel is more common neighbor.
            if WinningID == i
                LosingID = OverlapVoxIDs(k);
            else
                LosingID = i;
            end  
            VoxIDMap(OverlapVoxInds(k)) = WinningID;
            LosingVox = find(Dots.Vox(LosingID).Ind == OverlapVoxInds(k));
            Dots.Vox(LosingID).Pos(LosingVox,:) = [];
            Dots.Vox(LosingID).Ind(LosingVox) = [];
            Dots.Vox(LosingID).RawBright(LosingVox) = [];
            Dots.Vox(LosingID).IT(LosingVox) = [];
            Dots.Vol(LosingID) = Dots.Vol(LosingID)-1;
            Dots.ITMax(LosingID) = max(Dots.Vox(LosingID).IT);
            Dots.ItSum(LosingID) = sum(Dots.Vox(LosingID).IT);
            Dots.MeanBright(LosingID) = mean(Dots.Vox(LosingID).RawBright); 
        end
    else                
        VoxMap(Dots.Vox(i).Ind) = 1;
        VoxIDMap(Dots.Vox(i).Ind) = i;
    end
end

TotalNumOverlapDots
TotalNumOverlapVoxs
    
Dots.TotalNumOverlapDots = TotalNumOverlapDots;
Dots.TotalNumOverlapVoxs = TotalNumOverlapVoxs;

save([TPN 'Dots.mat'], 'Dots');
    
toc       

Settings.dotfinder.Gmode = Gmode;                   % added the saving of Gmode to check if Gmode is not too high. 1/18/2010
Settings.dotfinder.Gmax = Gmax;                     % Gmax for scaling the image if you didnt use the full bitdepth when creating tif stacks for 'I'
save([TPN 'Settings.mat'],'Settings')

%% nothing


