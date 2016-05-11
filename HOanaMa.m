%% Comments Errors
% 
% 091009 -47 'added lower boundry cache to SegMids.adam 
% 110509 -109 'added dots in mask if statement to remove unecessary
%   function dist runs'.adam
% 110509 -115 'corrected mod() percent dots done to be accurate
%   percent'.adam
% 120109 -62 'changed crap region to 3 instead of 2'. adam 


% This program removes dendritic segments in crap regions assigned as 3
% duing masking in Amira.
% Then, go through all the PSD95 dots, and calculate distance from the
% position of dots to the mask and distance from the position of dots to
% dendritic segments. If dots are within the mask, distance to the mask is
% assumed to be zero. So, with our new methods, which find dots only within
% the mask, mDotToMaskDist part is no longer necessary. HO 1/5/2010
%
% HO 6/25/2010 Removed Dots.cut and AllSegCut and AllSegBackup. 
% They are not necessary with our new method (we have no junk region in the
% mask which was supposed to cut from the AllSeg and the remaining segments
% saved in AllSegCut, plus cut dots registered in Dots.cut in the old 
% method). The old method was useful for example to remove dots in cell 
% bodies by assigning the cell bodies as crap region.
% The old method was kept as HOanaMaOriginal.mat.
%
% HO 7/30/2010 Save DDm.mat in TPN instead of saving bunch of tif files for
% each z plane under TPN\pics\DDm. pics folder had only this folder and
% files, so was removed (not generated in RunCell any more). Now voxels in
% DDm have 2 for the dendritic nodes (points in AllSeg) and 1 in between
% nodes, 0 for elsewhere.

%% 
function[]=HOanaMa(TPN)
%% Apply Mask to Data

colormap gray(255)

%% Load Dots and Dendrites
load([TPN 'Dots.mat'])
load([TPN 'Settings.mat'])
ImageInfo = Settings.ImInfo;
xyum=ImageInfo.xyum; %changed to reflect structure format of ImInfo HO 1/5/2010
zum=ImageInfo.zum; %changed to reflect structure format of ImInfo HO 1/5/2010
DotPos=Dots.Pos;
DotPos(:,1:2)=DotPos(:,1:2)*xyum;
DotPos(:,3)=DotPos(:,3)*zum;


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'
load([TPN 'D.mat'])
[ys xs zs]=size(D);


%% List mask
clear mask
mask=[0 0 0];
D=logical(D==1);

%Next 3 lines will remove not only isolated voxels but also remove voxels 
%having only 1 neighboring voxel within 3*3*3 space, but this is not
%necessary with our new masking method 1/5/2010 HO
% 'Counting surround'
% D=CountSurround(D);
% 'Done With Surround'
% D=D>1;  %remove small objects 
for i=1:size(D,3)
    [y x z]=ind2sub(size(D),find(D(:,:,i)));
    z=(z>0)*i;
    mask=cat(1,mask,[y x z]);
end
mask=mask(2:size(mask,1),:);
mask(:,1:2)=mask(:,1:2)*xyum;
mask(:,3)=mask(:,3)*zum;
% mask becomes the list of all the masking voxels (n*3 matrix, n is the number
% of voxels) written in um.

%% Check Dots against Mask
%%Find distance between dot and nearest masked pixel
[ys xs zs]=size(D);

%With our new maskign method, all the dots are within mask, and the
%distance to mask is always 0. 1/5/2010 HO
mDotToMaskDist = zeros(1,Dots.Num);

save([TPN 'mDotToMaskDist.mat'],'mDotToMaskDist')

Dots.DistToMask=mDotToMaskDist;
%change dotdata to those puncta within 2um
save([TPN 'Dots.mat'],'Dots')


%% Draw Dend and Calculate Dots.Dist2Dend
'Drawing Dendrites'
xyum = xyum;
zum = zum;
Scxy = (1/xyum);
Scz = (1/zum);
DDm = uint8(zeros(ys,xs,zs));
%%Draw Segments
SkelRes=.1; %subdivide dendritic segments into 0.1um segments
load([TPN 'data\AllSeg.mat'])
AllSegCut = AllSeg; %to take advantage of already-written program, just change the name from AllSeg to AllSegCut HO 6/25/2010
for i=1:size(AllSegCut,1)
        Dist=sqrt((AllSegCut(i,1,1)-AllSegCut(i,1,2))^2 + (AllSegCut(i,2,1)-AllSegCut(i,2,2))^2 + (AllSegCut(i,3,1)-AllSegCut(i,3,2))^2); %find distance
        Length(i)=Dist;
          devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
        for d=1:devs+1
            sy=AllSegCut(i,1,1)+((AllSegCut(i,1,2)-AllSegCut(i,1,1))/devs)*(d-1);
            sx=AllSegCut(i,2,1)+((AllSegCut(i,2,2)-AllSegCut(i,2,1))/devs)*(d-1);
            sz=AllSegCut(i,3,1)+((AllSegCut(i,3,2)-AllSegCut(i,3,1))/devs)*(d-1);
            
            syVox = ceil(sy*Scxy); sxVox = ceil(sx*Scxy); szVox = ceil(sz*Scz); %ceil because um scale starts from 0 but voxel goes 1,2,3,...
            if syVox < 1, syVox = 1; end
            if sxVox < 1, sxVox = 1; end
            if szVox < 1, zVox = 1; end; %somehow some pts still get 0, fix it, HO 10/6/2010
            if syVox > ys, syVox = ys; end
            if sxVox > xs, sxVox = xs; end
            if szVox > zs, szVox = zs; end; %somehow some pts still get more than matrix dimention, fix it, HO 10/15/2011
           
            DDm(syVox, sxVox, szVox)=1; %draw Skel
        end
end
clear Dist

%%DrawNodes
for i=1:size(AllSegCut,1)
    for j=1:2;
        syVox = ceil(AllSegCut(i,1,j)*Scxy); sxVox = ceil(AllSegCut(i,2,j)*Scxy); szVox = ceil(AllSegCut(i,3,j)*Scz); %ceil because um scale starts from 0 but voxel goes 1,2,3,...
        if syVox < 1, syVox = 1; end;
        if sxVox < 1, sxVox = 1; end;
        if szVox < 1, szVox = 1; end; %somehow some pts still get 0, fix it, HO 10/6/2010
        if syVox > ys, syVox = ys; end;
        if sxVox > xs, sxVox = xs; end;
        if szVox > zs, szVox = zs; end; %somehow some pts still get more than matrix dimention, fix it, HO 10/15/2011

        DDm(syVox, sxVox, szVox)=2; %draw nodes
    end
end
   

%%linearize DDm
[yDend,xDend,zDend] = ind2sub(size(DDm),find(DDm>0));
dendMask = cat(2,round(yDend*xyum),round(xDend*xyum),round(zDend*zum));


%%Calculate Distance of Dots to Dend

clear Dist mDot2DendDist
d = 0;

for i = 1:Dots.Num
    if DDm(sub2ind(size(DDm),Dots.Pos(i,1),Dots.Pos(i,2),Dots.Pos(i,3)))>0 % if dot is inside dendrite mask (DDm == 1), do not need to run function 'dist' to know it is 0
        mDot2DendDist(i) = 0; 
    else    
      Dist=dist(dendMask,DotPos(i,:)); 
      mDot2DendDist(i)=min(Dist);
      if mod(i,Dots.Num/10)==0 % corrected to represent real percent values {was mod(i,100) ==0}
      PercentDoneWithDots = double(i)/Dots.Num*100
      end
    end
end

save([TPN 'mDot2DendDist.mat'],'mDot2DendDist')

Dots.Dist2Dend=mDot2DendDist;
%change dotdata to those puncta within 2um
save([TPN 'Dots.mat'],'Dots')

%HO 7/30/2010 No need to do this because DistToMask is zero with a new way.
%%Draw Dots
% for i=1:Dots.Num
%         DDm(round(DotPos(i,1)*Scxy)+1,round(DotPos(i,2)*Scxy)+1,round(DotPos(i,3)*Scz)+1)=Dots.DistToMask(i)*10;
% end


colormap colorcube
image(max(DDm,[],3)*100) %30 to 100 HO 7/30/2010
%HO 7/30/2010 instead of saving bunch of tif files, save DDm.mat. DDm.mat
%will have 1 in dendritic skeleton, 2 at the nodes (points in AllSeg), 0 in
%elsewhere.
%imwriteNp(TPN,DDm,'DDm')
save([TPN 'DDm.mat'],'DDm');








