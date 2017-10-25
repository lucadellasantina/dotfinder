%% Find objects using an iterative thresholding approach in a masked volume
% -------------------------------------------------------------------------
% The object finder searches for individual objects within the volume by
% first dividing the volume in multiple sub-blocks and then applying the
% following steps to each search block:
%
% Step 1: Iteratively scan the volume to create a scoring map
% Step 2: Segment objects with multiple scoring peaks using watershed
% Step 3: Calculate volume, position, scoring and intensity of each object
% Step 4: Resolve objects laying across multiple search blocks 
% Tested object are synaptic puncta (ellipsoids)
% 
% Version 2.0                               2017-08-02  Luca Della Santina
%
% + Split computation in 4 fundamental steps
% + Improved ~15% speed by switching from bwlabeln to conncomp+labelmatrix
% % Relevant settings are passed to the function as parameters
%
% Version 1.0 formerly known as DotFinder   2010-xx-xx  Haruhisa Okawa
%                                           2010-xx-xx  Adam Bleckert
%                                           2008-xx-xx  Josh Morgan
%
% input: TPN: path to working directory
%        Post: 3D image stack of channel to search for objects
%        D: binary mask to restrict the search within
%        Settings: settings for the dofinding process
% output: Dots (and relative informations saved to TPN path)
% dependencies: txtBar.m
%               image processing toolbox
%
% -------------------------------------------------------------------------

function[Dots]=findObjects(Post, D, Settings)

disp('** Finding objects in the masked volume **');
debug=0;
[Rys, Rxs, Rzs] = size(Post);
if debug subplot(2,1,1),colormap gray(255);image(max(Post,[],3)),pause(.01); end

% Find intensity of noise in the image (Gmode),
% as most common pixel intensity other than zero.
% Dotfinder will search only above the noise intensity level
% up to the maximum intensity in the signal (Gmax)
Gmode = double(mode(Post(Post>0)));
Gmax = double(max(Post(Post>0)));

% Apply mask before dot finding 1/5/2010 HO
D = uint8(D);
Post = Post.*D;
if debug subplot(2,1,2),colormap gray(255),image(max(Post,[],3)),pause(1); end

% Enter Variables
v = Settings.dotfinder;

% Calculate the block size to Subsample the volume
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

% -- STEP 1: scan the volume and find areas crossing local contrast threshold with a progressively coarser intensity filter --
clear tmpBlocks;
tmpBlocks(NumBx, NumBy, NumBz) = struct;

txtBar('   STEP 1/4: Search candidate dots using iterarive threshold ... ');
tic;
for Bz=1:NumBz
    for By=1:NumBy
        for Bx=1:NumBx
            txtBar(100*((Bz-1)*NumBy*NumBx+(By-1)*NumBx+Bx)/(NumBz*NumBx*NumBy));
            
            %Find real territory
            Tystart=(By-1)*Byc+1;
            Txstart=(Bx-1)*Bxc+1;
            Tzstart=(Bz-1)*Bzc+1;
            if By<Byc, Tyend=By*Byc; else Tyend=Rys; end
            if Bx<Bxc, Txend=Bx*Bxc; else Txend=Rxs; end
            if Bz<Bzc, Tzend=Bz*Bzc; else Tzend=Rzs; end
            
            %Find buffered Borders (if last block, extend to image borders)
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
            
            % Slice the raw image into a block (Igm) to analyze and initialize matrixes based on its size
            Igm = single(Post(yStart:yEnd,xStart:xEnd,zStart:zEnd));
            [ys,xs,zs]=size(Igm);
            peakMap = zeros(ys,xs,zs,'uint8');      % Initialize matrix to map peaks found
            thresholdMap = zeros(ys,xs,zs,'uint8'); % Initialize matrix to sum passed thresholds
            maxIntensity = uint16(max(Igm(:)));     % Maximum intensity value found in the block
            
            % Make threshold steps same between different blocks of image 1/5/2010 HO
            if mod(maxIntensity, v.thresholdStep) ~= mod(Gmode+1, v.thresholdStep); % Iterate up to Gmode+1 for all blocks, Gmode+1 is fine and the program doesn't need to iterate up to obviously noise level.
                maxIntensity = maxIntensity+1;
            end
            
            % -- STEP 1: scan the volume and find areas crossing local contrast threshold with a progressively coarser intensity filter --
            for i = maxIntensity:-v.thresholdStep:Gmode+1 % Iterate from maxIntensity down to noise level (Gmode+1) within the block
                % Label all areas in the block (Igl) that crosses the intensity threshold "i"
                %[Igl,labels] = bwlabeln(Igm>i,6); OLD version, slower
                CC = bwconncomp(Igm>i,6); % using bwconncomp+labelmatric instead of bwlabeln increases speed about 10% LDS 2017-08-02
                labels = CC.NumObjects;
                Igl = labelmatrix(CC);
                if labels == 0
                    continue;
                elseif labels <= 1
                    labels = 2;
                end
                if labels<65536 Igl=uint16(Igl);end % Reduce bitdepth if possible
                
                nPixel = hist(Igl(Igl>0), 1:labels);
                
                for p=1:labels
                    pixelIndex = find(Igl==p);
                    
                    % Adjust max dot size if multi-peak object is found HO 2/8/2011
                    NumPeaks = sum(peakMap(pixelIndex));
                    if NumPeaks == 0
                        maxCompoundDotSize = v.maxDotSize;
                    else
                        maxCompoundDotSize = v.maxDotSize + v.maxDotSize*(NumPeaks-1)*v.MultiPeakDotSizeCorrectionFactor;
                    end
                    
                    % Identify the peak location in each labeled object whose size is within min and max DotSize
                    if (nPixel(p) < maxCompoundDotSize) && (nPixel(p) > v.minDotSize)
                        if sum(peakMap(pixelIndex))== 0 % sets up critireon to limit one peak per labeled field "Igl"
                            peakValue = max(Igm(pixelIndex));
                            peakIndex = find(Igl==p & Igm==peakValue);
                            if numel(peakIndex) > 1
                                peakIndex = peakIndex(round(numel(peakIndex)/2));
                            end
                            [y,x,z] = ind2sub([ys xs zs], peakIndex);
                            peakMap(y,x,z) = 1; % Register the peak position in peak map
                        end
                    else
                        Igl(pixelIndex)=0;
                    end
                end
                thresholdMap(Igl>0) = thresholdMap(Igl>0)+1; % Add 1 to the threshold score of all peaks that passed this iteration
                
            end % for all intensities
            thresholdMap(thresholdMap<v.itMin) = 0;             % only puncta with more than 2 IT pass are analyzed, I set to 2 for now because you want to save voxels of IT = 2 for the dot of ITmax = 4, for example. HO 2/8/2011
            tmpBlocks(Bx,By,Bz).sizeIgm = [ys, xs, zs];                 % Store information for later
            tmpBlocks(Bx,By,Bz).startPos = [yStart, xStart, zStart];    % Store information for later
            tmpBlocks(Bx,By,Bz).endPos = [yEnd, xEnd, zEnd];            % Store information for later
            tmpBlocks(Bx,By,Bz).thresholdMap = thresholdMap;            % Store information for later
            tmpBlocks(Bx,By,Bz).peakMap = peakMap;                      % Store information for later
            
        end % for all x blocks
    end % for all y blocks
end % for all z blocks
txtBar(['DONE in ' num2str(toc) ' seconds']);

%%
txtBar('   STEP 2/4: Watershed segmentation of multi-peak dots ... ');
tic;
% -- STEP 2: Find the countour of each dot and split it using watershed if multiple peaks are found within the same dot --
for Bz=1:NumBz
    for By=1:NumBy
        for Bx=1:NumBx
            txtBar(100*((Bz-1)*NumBy*NumBx+(By-1)*NumBx+Bx)/(NumBz*NumBx*NumBy));
            ys = tmpBlocks(Bx,By,Bz).sizeIgm(1);                              % retrieve stored values
            xs = tmpBlocks(Bx,By,Bz).sizeIgm(2);                              % retrieve stored values
            zs = tmpBlocks(Bx,By,Bz).sizeIgm(3);                              % retrieve stored values
            thresholdMap = tmpBlocks(Bx, By, Bz).thresholdMap;                % retrieve stored values
            
            wsThresholdMapBin = uint8(thresholdMap>0);                        % binary map of threshold
            wsThresholdMapBinOpen = imdilate(wsThresholdMapBin, ones(3,3,3)); % dilate Bin map with a 3x3x3 kernel (dilated perimeter acts like ridges between background and ROIs
            wsThresholdMapComp = imcomplement(thresholdMap);                  % complement (invert) image. Required because watershed() separate holes, not mountains. imcomplement creates complement using the entire range of the class, so for uint8, 0 becomes 255 and 255 becomes 0, but for double 0 becomes 1 and 255 becomes -254.
            wsTMMod = wsThresholdMapComp.*wsThresholdMapBinOpen;              % Force background outside of dilated region to 0, and leaves walls of 255 between puncta and background.
            wsTMLabels = watershed(wsTMMod, 6);                               % 6 voxel connectivity watershed, this will fill background with 1, ridges with 0 and puncta with 2,3,4,... in double format
            wsBackgroundLabel = mode(double(wsTMLabels(:)));                  % calculate background level
            wsTMLabels(wsTMLabels == wsBackgroundLabel) = 0;                  % seems that sometimes Background can get into puncta... so the next line was not good enough to remove all the background labels.
            wsTMLabels = double(wsTMLabels).*double(wsThresholdMapBin);       % masking out non-puncta voxels, this makes background and dilated voxels to 0. This also prevents trough voxels from being added back somehow with background. HO 6/4/2010
            wsTMLZeros = find(wsTMLabels == 0 & thresholdMap > 0);            % find zeros of watersheds inside of thresholdmap (add back the zero ridges in thresholdMap to their most similar neighbors)
            
            %LDS figure out why X and Y are stored inverted, Y first then X
            if ~isempty(wsTMLZeros) % if there exist zeros in the map
                [wsTMLZerosY,wsTMLZerosX,wsTMLZerosZ] = ind2sub(size(thresholdMap),wsTMLZeros); %6/4/2010 HO
                clear nZeroID;
                for j = 1:length(wsTMLZeros) % create a dilated matrix to examine neighbor connectivity around the zero position
                    tempZMID =  wsTMLabels(max(1,wsTMLZerosY(j)-1):min(ys,wsTMLZerosY(j)+1), max(1,wsTMLZerosX(j)-1):min(xs,wsTMLZerosX(j)+1), max(1,wsTMLZerosZ(j)-1):min(zs,wsTMLZerosZ(j)+1)); %HO 6/4/2010
                    nZeroID = mode(tempZMID(tempZMID~=0));                                      % find most common neighbor value (watershed) not including zero
                    wsTMLabels(wsTMLZeros(j)) = nZeroID;                                        % re-define zero with new watershed ID (this process will act similar to watershed by making new neighboring voxels feed into the decision of subsequent zero voxels)
                end
            end
            wsTMLabels = uint16(wsTMLabels);
            wsLabelList = unique(wsTMLabels);
            wsLabelList(1) = []; %remove background (now labeled as 0) from the list
            nLabels = length(wsLabelList);
            
            tmpBlocks(Bx, By, Bz).nLabels = nLabels; % Store information for the next step
            tmpBlocks(Bx, By, Bz).wsTMLabels = wsTMLabels; % Store information for the nect step
            tmpBlocks(Bx, By, Bz).wsLabelList = wsLabelList; % Store information for the nect step
        end
    end
end
txtBar(['DONE in ' num2str(toc) ' seconds']);

%%
txtBar('   STEP 3/4: Accumulating information of valid dots ... ');
tic;
% -- STEP 3: calculate properties of each dot and store them into Dots struct array --
clear tmpDots;
tmpDots.Pos=[0,0,0];
tmpDots.Vox.Pos=[0,0,0];
tmpDots.Vox.Ind=[0,0,0];
tmpDots.Vol = 0;
tmpDots.ITMax = 0;
tmpDots.ItSum = 0;
tmpDots.Vox.RawBright = 0;
tmpDots.Vox.IT = 0;
tmpDots.MeanBright = 0;

for Bz=1:NumBz
    for By=1:NumBy
        for Bx=1:NumBx
            txtBar(100*((Bz-1)*NumBy*NumBx+(By-1)*NumBx+Bx)/(NumBz*NumBx*NumBy));
            
            ys = tmpBlocks(Bx,By,Bz).sizeIgm(1);             % get stored values
            xs = tmpBlocks(Bx,By,Bz).sizeIgm(2);             % get stored values
            zs = tmpBlocks(Bx,By,Bz).sizeIgm(3);             % get stored values
            yStart = tmpBlocks(Bx,By,Bz).startPos(1);        % get stored values
            xStart = tmpBlocks(Bx,By,Bz).startPos(2);        % get stored values
            zStart = tmpBlocks(Bx,By,Bz).startPos(3);        % get stored values
            yEnd = tmpBlocks(Bx,By,Bz).endPos(1);            % get stored values
            xEnd = tmpBlocks(Bx,By,Bz).endPos(2);            % get stored values
            zEnd = tmpBlocks(Bx,By,Bz).endPos(3);            % get stored values
            peakMap = tmpBlocks(Bx,By,Bz).peakMap;           % get stored values
            wsTMLabels = tmpBlocks(Bx,By,Bz).wsTMLabels;     % get stored values
            wsLabelList = tmpBlocks(Bx,By,Bz).wsLabelList;   % get stored values
            thresholdMap = tmpBlocks(Bx,By,Bz).thresholdMap; % get stored values
            nLabels = tmpBlocks(Bx,By,Bz).nLabels;           % get stored values
            Igm = single(Post(yStart:yEnd,xStart:xEnd,zStart:zEnd)); % get stored values
            
            for i = 1:nLabels
                peakIndex = find(wsTMLabels==wsLabelList(i) & peakMap>0); % this line adjusted for watershed HO 6/7/2010
                thresholdPeak = thresholdMap(peakIndex);
                nPeaks = numel(peakIndex);
                
                [yPeak xPeak zPeak] = ind2sub([ys xs zs], peakIndex);
                
                if nPeaks == 1
                    %remove all dots found in buffer region instead of those found in the OUTER HALF of buffer region 1/5/2010 HO
                    if  (yPeak <= v.blockBuffer && yStart > 1) ||...
                            (xPeak <= v.blockBuffer && xStart > 1) ||...
                            (zPeak <= v.blockBuffer && zStart > 1) ||...
                            (ys-yPeak < v.blockBuffer && yEnd < Rys) ||...
                            (xs-xPeak < v.blockBuffer && xEnd < Rxs) ||...
                            (zs-zPeak < v.blockBuffer && zEnd < Rzs) ||...
                            (thresholdPeak < v.minFinalDotITMax) %HO 2/8/2011 added excluding dots that did not reach minFinalDotITMax criterion
                        %disp('dot in buffer region');
                    else
                        % Changed to define the cutOff range by yourself (upper and lower-bound rather than just defining the upper bound previously) HO 2/12/2010
                        PossibleMaxITMax = (Gmax - Gmode)/v.thresholdStep; %HO 2/12/2010 AAB 05/
                        thresholdPeak = double(thresholdPeak); %HO 3/1/2010 if thresholdPeak is uint8, the next line always spit out zero, so thresholdPeak must be double.
                        cutOff = (v.peakCutoffUpperBound - (v.peakCutoffUpperBound-v.peakCutoffLowerBound)*thresholdPeak/PossibleMaxITMax)*thresholdPeak; %HO 2/12/2010
                        thresholdPeak = uint8(thresholdPeak); %HO 3/1/2010 convert thresholdPeak back to uint8.
                        contourIndex = find(wsTMLabels==wsLabelList(i) & thresholdMap>=cutOff); %this line adjusted for watershed HO 6/7/2010
                        
                        if numel(contourIndex) >= v.minFinalDotSize %HO 2/15/2011 added excluding dots that did not reach minFinalDotSize criterion
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
                            
                            %LDS figure out why X and Y are stored inverted, Y first then X
                            tmpDot = struct;
                            tmpDot.Pos = [yPeak+yStart-1, xPeak+xStart-1, zPeak+zStart-1];
                            tmpDot.Vox.Pos = [yContour+yStart-1, xContour+xStart-1, zContour+zStart-1];
                            tmpDot.Vox.Ind = sub2ind([Rys Rxs Rzs], tmpDot.Vox.Pos(:,1), tmpDot.Vox.Pos(:,2), tmpDot.Vox.Pos(:,3));
                            tmpDot.Vol = numel(contourIndex);
                            tmpDot.ITMax = thresholdPeak;
                            tmpDot.ItSum = sum(thresholdMap(contourIndex));
                            tmpDot.Vox.RawBright = Igm(contourIndex);
                            tmpDot.Vox.IT = thresholdMap(contourIndex); %save also the iteration of each voxel 2/9/2011 HO
                            tmpDot.MeanBright = mean(Igm(contourIndex));
                            tmpDots(end+1) = tmpDot; % LDS figure out why this does not work with parfor
                            
                        end % if dot size bigger than min dot size
                    end % if dot does not reside in the border region between blocks
                end  % if only 1 peak found
            end  % for all labels
        end
    end
end
tmpDots(1)=[]; % Remove first tmpDot used for initialization of the struct array

% Convert tmpDots into the deprecated "Dots" structure
for i=1:numel(tmpDots)
    Dots.Pos(i,:) = tmpDots(i).Pos;
    Dots.Vox(i).Pos = tmpDots(i).Vox.Pos;
    Dots.Vox(i).Ind = tmpDots(i).Vox.Ind;
    Dots.Vol(i) = tmpDots(i).Vol;
    Dots.ITMax(i) = tmpDots(i).ITMax;
    Dots.ItSum(i) = tmpDots(i).ItSum;
    Dots.Vox(i).RawBright = tmpDots(i).Vox.RawBright;
    Dots.Vox(i).IT = tmpDots(i).Vox.IT;
    Dots.MeanBright(i) = tmpDots(i).MeanBright;
end
Dots.ImSize = [Rys Rxs Rzs];
Dots.Num = size(Dots.Pos,1);
txtBar(['DONE in ' num2str(toc) ' seconds']);

%%
txtBar('   STEP 4/4: Resolving duplicate dots on the border of search blocks ... ');
tic;
% -- STEP 4: resolve dots spanning across border region of analyzed blocks --
% Some voxels could be shared by multiple dots, which would not happen if the stack was not divided into sub-matrix.
% This happens around the dividing lines of sub-matrices.
% Instead of removing less brighter overlapping dot later in Singles criterion in anaSG, re-assign those voxels to a single dot. 2/9/2011 HO
VoxMap = uint8(zeros(Dots.ImSize));
VoxIDMap = zeros(Dots.ImSize);
[ys xs zs] = size(VoxMap);
TotalNumOverlapDots = 0;
TotalNumOverlapVoxs = 0;
for i = 1:Dots.Num
    txtBar(100*(i/Dots.Num));
    OverlapVoxInd = find((VoxMap(Dots.Vox(i).Ind) > 0));
    if ~isempty(OverlapVoxInd)
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
            SurroudingIDs(isnan(SurroudingIDs)) = 0 ; % Convert NaN to 0 if present in the matrix LDS fix 7-25-2017
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

Dots.TotalNumOverlapDots = TotalNumOverlapDots;
Dots.TotalNumOverlapVoxs = TotalNumOverlapVoxs;
txtBar(['DONE in ' num2str(toc) ' seconds']);

Settings.dotfinder.Gmode = Gmode;      % added the saving of Gmode to check if Gmode is not too high. 1/18/2010
Settings.dotfinder.Gmax = Gmax;        % Gmax for scaling the image if you didnt use the full bitdepth when creating tif stacks for 'I'
save([Settings.TPN 'Settings.mat'],'Settings')

clear Bx* By* Bz* CC contour* cutOff debug Gm* i j k Ig* labels Losing*
clear max* n* Num* Overlap* p peak* Possible* Rxs Rys Rzs Surrouding*
clear temp* tmp* threshold* Total* Tx* Ty* Tz* v Vox* Winning* ws* x* y* z*
end