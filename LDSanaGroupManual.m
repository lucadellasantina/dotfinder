function[] = LDSanaGroupManual(TPN)

%1/15/2010 HO just added comments to make it easy to understand each line
%6/28/2010 HO changed nearby from 10 voxels to 1/Settings.ImInfo.xyum so
%that I can process images with different xy resolutions.
%6/29/2010 HO converted the structural elements in Grouped the same as Dots
%and added the same structural elements of Dots in the Grouped except those
%obtained by anaRa, anaMa, anaCB, anaRd. Now Grouped is similar to Dots but
%just grouping some dots together so that it's easy to understand. Plus,
%you can feed Grouped into anaMa, anaCB and anaRd to get those parameters
%for grouped dots.
%6/29/2010 HO removed repeating the drawing of groupI 10 times in the last
%part of the program.
%The original HOanaGroup was saved as HOanaGroupOriginal.m.
%1/5/2011 to 1/7/2011 HO added another criterion to group, which is based
%on how much raw brightness difference exists between MaxRawBright of
%dimmer dot and the brightest voxel among those facing each other. HO also
%added manual decision of grouping for dots that are in gray zone (can be
%grouped or not grouped). If the number of dots in gray zone is not many,
%like PSD95 dots, you can do it all manually. For dots like CtBP2, there
%are too many dots, so I made the sampling of some dots in the gray zone, 
%and you can decide where to draw lines to separate grouped vs non-grouped
%dots based on the distribution of sampled dots in the graph. HO also changed 
%the definition of contactRatio to compare face vs the volume of smaller 
%dot, not the average of both dots. 
%This new version was saved as HOanaGroupManual. The previous version was 
%left as HOanaGroup.
%6/16/2011 HO introduced SG.passI (using Imaris for final yes no), so use
%passI if it exists, otherwise use passF.

%% 
%colormap colorcube

%load data
load([TPN 'Dots.mat'])
load([TPN 'Settings.mat']) %HO 6/28/2010
load([TPN 'find\SG.mat'])

load([TPN 'Post.mat']) %HO 1/5/2011

%collect information from dots that passed SG thresholding
if isfield(SG,'passI')
    Pass = SG.passI;
else
    Pass = SG.passF;
end
IDs=find(Pass);
if size(IDs,1) == 1; %IDs need to be n*2 matrix for the rest of the code to work.
    IDs = IDs';
end

ImSize=Dots.ImSize;
Vox = Dots.Vox;

nearby = round(1/Settings.ImInfo.xyum);%10;  %Distance at which faces are looked at %look up to 1um nearby, changed from 10 to 1/Settings.ImInfo.xyum to deal with images with different resolution HO 6/28/2010
%maxContact = 0.3; % maximum faces/volume to be independent dot

%% Determine distance from dot to dot
Dists = zeros(length(IDs),length(IDs));
for i = 1: length(IDs)
   Dists(i,:) = dist(Dots.Pos(IDs,:),Dots.Pos(IDs(i),:)); %[i,j] element in Dists will be the distance between i-th dot and j-th dot in a unit of voxel number.
  
end
[x y] = find((Dists > 0) & (Dists < nearby)); %pick pair of dots whose distance between is >0 and <nearby.

%The following 3 lines remove redundant mirror pairs 1/5/2011 HO
x2 = x(find(x<y)); 
y2 = y(find(x<y));
x=x2;y=y2;clear x2 y2;

idX=IDs(x); %idX will be IDs of one dot in each pair of dots
idY=IDs(y); %idY will be IDs of the other dot in each pair of dots

%% Count shared faces of nearby dots
faces = [];
facesMaxRawBright = []; %1/6/2011 HO
for i = 1: length(idX) %go through individual idX-idY pair
   cDotX = Dots.Vox(idX(i)).Pos;
   cDotY = Dots.Vox(idY(i)).Pos;
   face = zeros(size(cDotX,1),1);
   faceRawBright = zeros(size(cDotX,1),1); %1/6/2011 HO
   for v = 1: size(cDotX,1) %go through individual voxels in the dot specified by idX(i) and calculate faces with all the voxels in the dot specified by idY(i)
        diff = [cDotY(:,1)-cDotX(v,1)  cDotY(:,2)-cDotX(v,2) cDotY(:,3)-cDotX(v,3)];
        face(v) = sum(sum(abs(diff),2)<2); %this condition will be met only when one voxel sits right next to the other in 6 primany neighboring directions

        %HO added the following to extract the max raw bright among facing
        %voxels. 1/6/2011
        facingVoxNumDotY = find(sum(abs(diff),2)<2);
        if ~isempty(facingVoxNumDotY)
            facingVoxRawBrightDotY = Dots.Vox(idY(i)).RawBright(facingVoxNumDotY);
            facingVoxRawBrightDotX = Dots.Vox(idX(i)).RawBright(v);
        else
            facingVoxRawBrightDotY = 0;
            facingVoxRawBrightDotX = 0;
        end
        faceRawBright(v) = max([facingVoxRawBrightDotY; facingVoxRawBrightDotX]);
    
   end  %so face(v) counts number of voxels in the idY(i) dot that faces a voxel 'v' in the idX(i) dot
    faces(i) = sum(face); %faces(i) will be the sum of all of these countings for idX(i) dot
    
    facesMaxRawBright(i) = max(faceRawBright); %1/6/2011 HO
end

for i = 1:length(faces);
   %Vol = Dots.Vol(idX(i)); %Josh original Vol definition
   Vol = min(Dots.Vol(idX(i)), Dots.Vol(idY(i))); %HO 1/5/2011 new definition
   %Vol = (Dots.Vol(idX(i))+Dots.Vol(idY(i)))/2; %HO 6/8/2010 definition
   contactRat(i) = faces(i)/Vol; %calculate the ratio of faces over the volume of the puncta
   
   %HO added another criterion 1/6/2011
   DiffRawBright(i) = min([max(Dots.Vox(idX(i)).RawBright), max(Dots.Vox(idY(i)).RawBright)]) - facesMaxRawBright(i);
   
end

plot(contactRat, DiffRawBright, 'k.')
xlabel('contactRat');ylabel('DiffRawBright');

%% form new groups if too much contact
%HO added 1/8/2011
GroupParam.ContactUpperBound = 0.35;
GroupParam.ContactLowerBound = 0.03;
GroupParam.DiffRawBrightUpperBound = 8.5;
GroupParam.DiffRawBrightLowerBound = 0.5;
if exist([TPN 'Grouped']);
    load([TPN 'Grouped']);
    if isfield(Grouped,'GroupParam')  % Check if GroupParam has been created
        if length(fieldnames(GroupParam)) == length(fieldnames(Grouped.GroupParam))
            GroupParam = Grouped.GroupParam;
        end
    end
    clear Grouped;
end
title = 'You can make Contact upper and lower bounds the same to avoid manual grouping.';
GroupParam = getVars(GroupParam,title);  % Let user define first criteria

ContactUpperBound = GroupParam.ContactUpperBound;
ContactLowerBound = GroupParam.ContactLowerBound;
DiffRawBrightUpperBound = GroupParam.DiffRawBrightUpperBound;
DiffRawBrightLowerBound = GroupParam.DiffRawBrightLowerBound;

GroupedFlag = zeros(1,length(faces));
GrayZoneFlag = zeros(1,length(faces));
ManualGroupingFlag = zeros(1,length(faces));
ManualGroupedFlag = zeros(1,length(faces));
ManualNonGroupedFlag = zeros(1,length(faces));

%newGroup = [];
for i = 1:length(faces);
   if ((contactRat(i) > ContactUpperBound) && (DiffRawBright(i) <= DiffRawBrightUpperBound)) || (DiffRawBright(i) <= DiffRawBrightLowerBound)
       %newGroup(size(newGroup,1)+1,:) = [idX(i) idY(i)];
       GroupedFlag(i) = 1;
   elseif (contactRat(i) > ContactLowerBound) && (contactRat(i) <= ContactUpperBound) && (DiffRawBright(i) > DiffRawBrightLowerBound) && (DiffRawBright(i) <= DiffRawBrightUpperBound)
       GrayZoneFlag(i) = 1;
   end 
end

plot(contactRat, DiffRawBright, 'k.')
hold on
plot(contactRat(contactRat>ContactUpperBound), DiffRawBright(contactRat>ContactUpperBound), 'r.')
plot(contactRat(contactRat<=ContactLowerBound), DiffRawBright(contactRat<=ContactLowerBound), 'g.')
plot(contactRat(DiffRawBright>DiffRawBrightUpperBound), DiffRawBright(DiffRawBright>DiffRawBrightUpperBound), 'g.')
plot(contactRat(DiffRawBright<=DiffRawBrightLowerBound), DiffRawBright(DiffRawBright<=DiffRawBrightLowerBound), 'r.')
xlabel('contactRat');ylabel('DiffRawBright');
legend('dots in gray zone', 'grouped dots', 'NOT grouped dots')
'Number of pairs of dots in the gray zone is'
length(find(GrayZoneFlag == 1))
%NumGrayZoneDots = length(find((contactRat >= ContactLowerBound) & (contactRat <= ContactUpperBound) & (DiffRawBright >= DiffRawBrightLowerBound) & (DiffRawBright <= DiffRawBrightUpperBound)))


%% Manual grouping added 1/8/2011 HO
%Make the list of dots that needs to be tested manually
if isempty(find(GrayZoneFlag == 1));
    ManualAllDotsFlag = 0;
else
    ManualAllDotsFlag = input('Manual grouping of dots in the gray zone. Type 1 for going through all the dots, 0 for sampling some dots and optimizing. \n');
    if ManualAllDotsFlag == 1
        ManualGroupingFlag = GrayZoneFlag;  
    else
        ManualAllDotsFlag = 2;
        ContactBoundaries = ContactLowerBound:0.03:ContactUpperBound+0.03;
        DiffRawBrightBoundaries = DiffRawBrightLowerBound:2:DiffRawBrightUpperBound+2;
        for i=1:length(ContactBoundaries)-1;
            for j=1:length(DiffRawBrightBoundaries)-1;
                SampleCandidates = find((GrayZoneFlag == 1) & (contactRat > ContactBoundaries(i)) & (contactRat <= ContactBoundaries(i+1)) & (DiffRawBright > DiffRawBrightBoundaries(j)) & (DiffRawBright <= DiffRawBrightBoundaries(j+1)));
                if ~isempty(SampleCandidates)
                    ManualGroupingFlag(SampleCandidates(1)) = 1; %just sample the first one
                end
            end
        end
    end
end

TotalNumDotsManualGrouping = length(find(ManualGroupingFlag == 1));
NumDotsDone = 0;

% cfigure(32,23);
% imageAX = axes('position',[.02 .02 .95 .95]);
% colormap gray(255);
close all;
i=1;
while i <= length(faces);
    if ManualGroupingFlag(i) == 1;
    %if (contactRat(i) >= ContactLowerBound) && (contactRat(i) <= ContactUpperBound) && (DiffRawBright(i) >= DiffRawBrightLowerBound) && (DiffRawBright(i) <= DiffRawBrightUpperBound)
       pDotX = Dots.Pos(idX(i),:);
       pDotY = Dots.Pos(idY(i),:);
       pDiff = pDotX - pDotY;
       pDiff = abs(pDiff);
       PostCuty1 = max(min(pDotX(1),pDotY(1))-floor((30-pDiff(1))/2), 1);
       PostCuty2 = min(max(pDotX(1),pDotY(1))+floor((30-pDiff(1))/2), ImSize(1));
       PostCutx1 = max(min(pDotX(2),pDotY(2))-floor((30-pDiff(2))/2), 1);
       PostCutx2 = min(max(pDotX(2),pDotY(2))+floor((30-pDiff(2))/2), ImSize(2));
       PostCutz1 = max(min(pDotX(3),pDotY(3))-floor((20-pDiff(3))/2), 1);
       PostCutz2 = min(max(pDotX(3),pDotY(3))+floor((20-pDiff(3))/2), ImSize(3));
       PostCut = Post(PostCuty1:PostCuty2, PostCutx1:PostCutx2, PostCutz1:PostCutz2);
       MaxRawBright = max([Dots.Vox(idX(i)).RawBright; Dots.Vox(idY(i)).RawBright]);
       ScalingFactor = double(round((255-mod(255, MaxRawBright))/MaxRawBright));
       PostCut = PostCut.*ScalingFactor;
       IDMapCut = PostCut.*0;
       IDMapCut1 = IDMapCut;
       IDMapCut2 = IDMapCut;
       VoxPosDotX = Dots.Vox(idX(i)).Pos;
       VoxPosDotY = Dots.Vox(idY(i)).Pos;
       VoxPosDotXShifted = [VoxPosDotX(:,1)-PostCuty1+1, VoxPosDotX(:,2)-PostCutx1+1, VoxPosDotX(:,3)-PostCutz1+1];
       VoxPosDotYShifted = [VoxPosDotY(:,1)-PostCuty1+1, VoxPosDotY(:,2)-PostCutx1+1, VoxPosDotY(:,3)-PostCutz1+1];
       VoxPosDotXShiftedInd = sub2ind(size(IDMapCut), VoxPosDotXShifted(:,1), VoxPosDotXShifted(:,2), VoxPosDotXShifted(:,3));
       VoxPosDotYShiftedInd = sub2ind(size(IDMapCut), VoxPosDotYShifted(:,1), VoxPosDotYShifted(:,2), VoxPosDotYShifted(:,3)); 
       IDMapCut1(VoxPosDotXShiftedInd) = 80;
       IDMapCut2(VoxPosDotYShiftedInd) = 175;
       IDMapCut = IDMapCut1 + IDMapCut2;
%        for slice = 1:10;
%            z = slice + 5; %start from 6th z plane
%            ImagePosy = mod(slice-1,5);
%            ImagePosx = (slice-1-ImagePosy)/5;
%            subplot('Position',[ImagePosx*2*0.16+0.02, 0.95-(ImagePosy+1)*0.18, 0.16, 0.18]);image(PostCut(:,:,z));axis off;colormap gray(255);
%            subplot('Position',[(ImagePosx*2+1)*0.16+0.02, 0.95-(ImagePosy+1)*0.18, 0.16, 0.18]);image(IDMapCut(:,:,z));axis off;colormap gray(255);
%        end

       ImStk = cat(2, PostCut, IDMapCut);
       colmap = 'gray(256)';
       FrmPerSec = 5;
       VideoWindowSize = [0 0.03 0.65 0.90];
       HOvideofig(size(ImStk,3), @(frm) HOredraw(frm, ImStk, colmap), FrmPerSec, [], [], VideoWindowSize); % if user bins is 10 this results in 100x real speed
       HOredraw(1, ImStk, colmap);

       set(gcf,'units','centimeters');
       % Place the figure
       set(gcf,'position',[1 3 25 17]);
       % Set figure units back to pixels
       set(gcf,'units','pixel');

       PostCutMaxZ = max(PostCut(:,:,min([VoxPosDotXShifted(:,3);VoxPosDotYShifted(:,3)]):max([VoxPosDotXShifted(:,3);VoxPosDotYShifted(:,3)])),[],3);
       PostCutMaxY = max(PostCut(min([VoxPosDotXShifted(:,1);VoxPosDotYShifted(:,1)]):max([VoxPosDotXShifted(:,1);VoxPosDotYShifted(:,1)]),:,:),[],1);
       PostCutMaxY = permute(PostCutMaxY, [3 2 1]);
       PostCutMaxX = max(PostCut(:,min([VoxPosDotXShifted(:,2);VoxPosDotYShifted(:,2)]):max([VoxPosDotXShifted(:,2);VoxPosDotYShifted(:,2)]),:),[],2);
       PostCutMaxX = permute(PostCutMaxX, [3 1 2]);
       IDMapCutMaxZ = max(IDMapCut1,[],3)+max(IDMapCut2,[],3);
       IDMapCutMaxY = max(IDMapCut1,[],1)+max(IDMapCut2,[],1);IDMapCutMaxY = permute(IDMapCutMaxY, [3 2 1]);
       IDMapCutMaxX = max(IDMapCut1,[],2)+max(IDMapCut2,[],2);IDMapCutMaxX = permute(IDMapCutMaxX, [3 1 2]);
       subplot('Position',[0.66, 0.77, 0.16, 0.18]);image(PostCutMaxZ);axis off;colormap gray(255);
       subplot('Position',[0.82, 0.77, 0.16, 0.18]);image(IDMapCutMaxZ);axis off;colormap gray(255);
       subplot('Position',[0.66, 0.59, 0.16, 0.18]);image(PostCutMaxY);axis off;colormap gray(255);
       subplot('Position',[0.82, 0.59, 0.16, 0.18]);image(IDMapCutMaxY);axis off;colormap gray(255);
       subplot('Position',[0.66, 0.41, 0.16, 0.18]);image(PostCutMaxX);axis off;colormap gray(255);
       subplot('Position',[0.82, 0.41, 0.16, 0.18]);image(IDMapCutMaxX);axis off;colormap gray(255);
       
       uicontrol('Style','text','Units','normalized','position',[.1,.96,.1,.02],'String','TotalNumDots');
       uicontrol('Style','text','Units','normalized','position',[.21,.96,.05,.02],'String',TotalNumDotsManualGrouping);
       uicontrol('Style','text','Units','normalized','position',[.29,.96,.1,.02],'String','NumDotsDone');
       uicontrol('Style','text','Units','normalized','position',[.4,.96,.05,.02],'String',NumDotsDone);
       uicontrol('Style','text','Units','normalized','position',[.5,.96,.1,.02],'String','contactRat');
       uicontrol('Style','text','Units','normalized','position',[.61,.96,.05,.02],'String',contactRat(i));
       uicontrol('Style','text','Units','normalized','position',[.69,.96,.1,.02],'String','DiffRawBright');
       uicontrol('Style','text','Units','normalized','position',[.8,.96,.05,.02],'String',DiffRawBright(i));
       uicontrol('Style','Pushbutton','Units','normalized','position',[.70,.2,.28,.08],...
            'String','Group','CallBack','ManualGroupedFlag(i) = 1; ManualNonGroupedFlag(i) = 0; uiresume');
       uicontrol('Style','Pushbutton','Units','normalized','position',[.70,.1,.28,.08],...
            'String','Not group','CallBack','ManualGroupedFlag(i) = 0; ManualNonGroupedFlag(i) = 1; uiresume');
        
       uiwait;
       
       close all;
       
       NumDotsDone = NumDotsDone + 1;
       
%        GroupingFlag = input('Type 1 for grouping, 0 for not grouping. \n');
%        if GroupingFlag == 1;
%            ManualGroupedFlag(i) = 1;
%            ManualNonGroupedFlag(i) = 0; %just to ensure when you repeat the program.
%            NumDotsDone = NumDotsDone + 1;
%        elseif GroupingFlag == 0;
%            ManualGroupedFlag(i) = 0; %just to ensure when you repeat the program.
%            ManualNonGroupedFlag(i) = 1;
%            NumDotsDone = NumDotsDone + 1;
%        else
%            i = i-1; %mistyping. go back to the same dot.
%        end
    end
    i=i+1; %go to the next dot
end

close all;
plot(contactRat, DiffRawBright, 'k.')
hold on
plot(contactRat(ManualGroupedFlag == 1), DiffRawBright(ManualGroupedFlag == 1), 'm.')
plot(contactRat(ManualNonGroupedFlag == 1), DiffRawBright(ManualNonGroupedFlag == 1), 'c.')
plot(contactRat(contactRat>ContactUpperBound), DiffRawBright(contactRat>ContactUpperBound), 'r.')
plot(contactRat(contactRat<=ContactLowerBound), DiffRawBright(contactRat<=ContactLowerBound), 'b.')
plot(contactRat(DiffRawBright>DiffRawBrightUpperBound), DiffRawBright(DiffRawBright>DiffRawBrightUpperBound), 'b.')
plot(contactRat(DiffRawBright<=DiffRawBrightLowerBound), DiffRawBright(DiffRawBright<=DiffRawBrightLowerBound), 'r.')
xlabel('contactRat');ylabel('DiffRawBright');
legend('dots in gray zone', 'dots in gray zone manually grouped', 'dots in gray zone manually NOT grouped', 'grouped dots', 'NOT grouped dots')
'Number of GROUPED pairs of dots in the gray zone is'
length(find(ManualGroupedFlag == 1))

if ManualAllDotsFlag == 2; %set arbitorary lines that divide grouped and not grouped dots.
    GroupGrayZoneParam.contactRatBorder = input('Where do you want to set a border line along contactRat axis for gray zone dots?\n');
    GroupGrayZoneParam.DiffRawBrightBorder = input('Where do you want to set a border line along DiffRawBright axis for gray zone dots?\n');
    
    ManualGroupedFlag((GrayZoneFlag == 1) & (contactRat > GroupGrayZoneParam.contactRatBorder) & (ManualNonGroupedFlag ~= 1)) = 1;
    ManualGroupedFlag((GrayZoneFlag == 1) & (DiffRawBright <= GroupGrayZoneParam.DiffRawBrightBorder) & (ManualNonGroupedFlag ~= 1)) = 1;
    plot(contactRat(ManualGroupedFlag == 1), DiffRawBright(ManualGroupedFlag == 1), 'ro');
    set(gca,'Title',text('String','Dots grouped in the gray zone were circled in red.'));
end

pause();
    
%% Organize groups
GroupedFlag = GroupedFlag + ManualGroupedFlag; %1/8/2011 HO
newGroup = [idX(GroupedFlag == 1) idY(GroupedFlag == 1)]; %1/8/2011 HO

depleteIDs= newGroup;
groups = {};
gNum = 0;
while sum(depleteIDs(:))
   gNum = gNum+1;
   newest = depleteIDs(find(depleteIDs>0,1)); %do the first one in the remaining list
   foundNew = 1;
   groups{gNum} = newest;
   while foundNew
       foundNew = 0;
       grby = [];
       for n = 1: length(newest)
            [grby grbx] = find(depleteIDs == newest(n)); %more than one newest(n) could be in deplateIDs, and grby and grbx will be the row and the column indices
            graby = [grby;grby]; %this will repeat the grby, make the column vector twice longer
       end
       if ~isempty(graby)
           grabbed = newGroup(graby,:); 
           grabbed = unique(grabbed(:)); %unique!? what is the grby repeat for?
           grabbed = setdiff(grabbed, newest); %find the ID of all the dots that have significant amt of facing with the 'newst' dot
           groups{gNum} = [groups{gNum} grabbed'];   
           depleteIDs(graby,:)=0; %convert already-checked dot pairs to 0 so that it avoids find(depleteIDs>0,1) line above
           newest = grabbed;
           foundNew = 1;
       end
   end
end

nonGroupedIDs = setdiff(IDs,newGroup(:));
for i = 1:length(nonGroupedIDs)
   groups{length(groups)+1} = nonGroupedIDs(i); 
end



%% Create merged dots
%6/29/2010 HO modified the structure of Grouped completely.
for i = 1:length(groups)
    Grouped.ids{i} = [];
    Grouped.Vox(i).Pos = [];
    Grouped.Vox(i).Ind = [];
    Grouped.Vox(i).RawBright = [];
end
Grouped.Pos = zeros(length(groups),3);
Grouped.Vol = zeros(1,length(groups));
Grouped.ITMax = zeros(1,length(groups));
Grouped.ItSum = zeros(1,length(groups));
Grouped.MeanBright = zeros(1,length(groups));
    
for i = 1:length(groups)
    Grouped.ids{i} = groups{i};
    for g= 1:length(groups{i})
        Grouped.Vox(i).Pos=[Grouped.Vox(i).Pos ; Dots.Vox(groups{i}(g)).Pos];
        Grouped.Vox(i).Ind=[Grouped.Vox(i).Ind ; Dots.Vox(groups{i}(g)).Ind];
        Grouped.Vox(i).RawBright=[Grouped.Vox(i).RawBright ; Dots.Vox(groups{i}(g)).RawBright];
    end
    Grouped.Pos(i,:)=mean(Grouped.Vox(i).Pos,1); %mean of all the voxels is Pos (so, it's center of mass.)
    Grouped.Vol(i) = length(Grouped.Vox(i).Ind); %6/29/2010 HO
    Grouped.ITMax(i) = max(Dots.ITMax(groups{i})); %6/29/2010 HO
    Grouped.ItSum(i) = sum(Dots.ItSum(groups{i})); %6/29/2010 HO
    Grouped.MeanBright(i) = mean(Grouped.Vox(i).RawBright); %6/29/2010 HO
end

Grouped.ImSize = Dots.ImSize; %6/29/2010 HO
Grouped.Num = length(groups); %6/29/2010 HO
Grouped.GroupParam = GroupParam; %1/6/2011 HO
Grouped.ManualAllDotsFlag = ManualAllDotsFlag % 1/8/2011 HO
Grouped.GroupedFlag = GroupedFlag % 1/8/2011 HO
Grouped.ManualGroupedFlag = ManualGroupedFlag % 1/8/2011 HO
Grouped.ManualNonGroupedFlag = ManualNonGroupedFlag % 1/8/2011 HO
if ManualAllDotsFlag == 2; % 1/8/2011 HO
    Grouped.GroupGrayZoneParam = GroupGrayZoneParam;
end

save([TPN 'Grouped.mat'],'Grouped')

FinalNumDots = Grouped.Num


%% Draw 
%for f = 1:10 %there is no point in repeating the same thing 10 times
%6/29/2010 HO
maxSum=zeros(Grouped.ImSize(1),Grouped.ImSize(2)); %Dots.ImSize to Grouped.ImSize 6/29/2010 HO
maxID1=maxSum; maxID2 = maxSum; maxID3 = maxSum;

%P=361
YXsize=Grouped.ImSize(1)*Grouped.ImSize(2); %Dots.ImSize to Grouped.ImSize 6/29/2010 HO

for i = 1:Grouped.Num %size(Grouped.DPos,1) to Grouped.Num 6/29/2010 HO
    Pos=Grouped.Vox(i).Pos; %Grouped.Vox to Grouped.Vox.Pos 6/29/2010 HO
    PosI=sub2ind(Grouped.ImSize(1:2),Pos(:,1),Pos(:,2)); %just ImSize to Grouped.ImSize 6/29/2010 HO
    maxID1(PosI)=rand*200+50;  
    maxID2(PosI)=rand*200+50;  
    maxID3(PosI)=rand*200+50;   
    maxSum(PosI)=maxSum(PosI)+1;
end

maxPassed2=uint8((maxSum>0)*200);

MaxC=maxID1+(maxSum>1)*1000;
MaxC(:,:,2)=maxID2+(maxSum>1)*1000;
MaxC(:,:,3)=maxID3+(maxSum>1)*1000;
MaxC=uint8(MaxC);
clf;
image(MaxC),pause(.01)
GroupI=MaxC;
%end

%% load Raw
% load([TPN 'images\maxRaw.mat']) 
% GroupI=maxRaw;
% GroupI(:,:,3)=maxSum*200;
imwrite(GroupI,[TPN 'find\GroupI.tif'],'Compression','none')
%image(GroupI),pause(.01)





