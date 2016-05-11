function[]=HOanaRdGrouped(TPN,DPN)

% HO 6/25/2010 flipped the division to calculate Round.Long, Round.Oblong
% and Round.Compact so that less spherical or less smooth ugly dots get
% lower values, then you can set threshold to take higher values in anaSG.
% However, Compact doesn't work well with small dots which have only
% several voxels because if you think 7-voxel sphere, 1 in the center and 6
% around the center voxel, the average 6-connectivity face per voxel will 
% be ~4.3, but if you think 2by2by2 square (8 voxels in total), the value
% will be 3, much lower than sphere. So, actually lots of small dots whose
% Compact values come out large with my flipping the division are noise, so
% the thresholding in anaSG using Compact wouldn't work. This is why the
% old definition somehow worked fine. Just don't include this criteria for
% minimum thresholding in anaSG.
% HO 7/6/2010 modified the generation of reference sphere. Since I use
% different xyz voxel dimentions for different types of imagings, I had to
% increase the volume of reference sphere to work with finer images. Then,
% the old code took forever to generate it, so I modified the code to make
% it faster.
%
% 7/6/2010 HO modified the code to be applicable to the grouped dots. It
% simply loads Grouped instead of Dots, convert the name from Grouped to
% Dots to take advantage of already written program, then at the end,
% convert it back to Grouped and save Grouped.


%Instead of loading Dots, load Grouped and change the Grouped to Dots to
%take advantage of already written program. 7/6/2010 HO
%load([TPN 'Dots.mat'])
load([TPN 'Grouped.mat'])
Dots = Grouped;

%% Find appropriate mean faces for perfect sphere as reference

%The initial 'for' loop was completely modified to make it faster. The
%result came out the same as the old way 7/6/2010
TSphere=zeros(31,31,31); %changed from 11*11*11 to 31*31*31 because I do 0.025um xy 0.2um z for the finest image of CtBP2 puncta (so 24 times more possible dot volume compared to 0.103um xy 0.3um z) 6/25/2010
TSphere2=zeros(33,33,33);
TSphere(16,16,16)=1;
Tdists=bwdist(TSphere);
for d=1:160 %chagned from 100 to 160, if you go >160, the sphere tries to get voxels outside the 31*31*31 3D matrix. HO 6/25/2010
   Near=find(Tdists<(d/10)); %d=1:10 will identify only the center point, then d=11 will identify 6 more voxels around the center voxel
   Tvol(d)=size(Near,1); %Tvol(d) will be the number of voxels within the distance of d/10 from the center voxel
   TSphere(Near)=1;
   TSphere2(2:end-1,2:end-1,2:end-1) = TSphere;
   Tperim = bwperim(TSphere,6);
   NearPerim = find(Tperim);
   [NearPerimY,NearPerimX,NearPerimZ] = ind2sub(size(TSphere), NearPerim);
   FaceCount=0;
   for n = 1:length(NearPerim) %go through each voxel in perim
       if TSphere2(NearPerimY(n), NearPerimX(n)+1,NearPerimZ(n)+1) == 0;
           FaceCount = FaceCount+1;
       end
       if TSphere2(NearPerimY(n)+2, NearPerimX(n)+1,NearPerimZ(n)+1) == 0;
           FaceCount = FaceCount+1;
       end
       if TSphere2(NearPerimY(n)+1, NearPerimX(n),NearPerimZ(n)+1) == 0;
           FaceCount = FaceCount+1;
       end
       if TSphere2(NearPerimY(n)+1, NearPerimX(n)+2,NearPerimZ(n)+1) == 0;
           FaceCount = FaceCount+1;
       end
       if TSphere2(NearPerimY(n)+1, NearPerimX(n)+1,NearPerimZ(n)) == 0;
           FaceCount = FaceCount+1;
       end
       if TSphere2(NearPerimY(n)+1, NearPerimX(n)+1,NearPerimZ(n)+2) == 0;
           FaceCount = FaceCount+1;
       end
   end
   meanFaces(d)=FaceCount/Tvol(d); 
end
c=0;
for v=1:max(Tvol)
   if ~isempty(find(Tvol==v))
       c=c+1;
       tvol(c)=v; %tvol will remove redundancy in Tvol, so tvol would be 1, 7, ...
       v2f(c)=meanFaces(find(Tvol==v,1)); %find(Tvol==v,1) will take only the 1st one among all Tvol==v, so again removing redundancy in meanFaces
   end    
end
RoundFaces=interp1(tvol,v2f,1:max(tvol)); %thinking that v2f is a function of tvol, interpolate between points with a step of tvol=1
%plot(1:max(tvol), RoundFaces, 'o') %HO 1/5/2010


%% Run Dots
for i = 1 :Dots.Num
       
   Cent=Dots.Pos(i,:);
   Vox=Dots.Vox(i).Pos;   
   Dist = dist(Vox,Cent);
   MeanD= max(1,mean(Dist));
   DistN=Dist/MeanD; %normalize by MeanD
   CentN=Cent/MeanD; %normalize by MeanD
   VoxN=Vox/MeanD; %normalize by MeanD
   
 %  MedN = median(DistN,1);
%   ExtN = max(DistN)/MedN;
   
   if size(Vox,1)>1, %if more then one voxel
       [Co, Sc,latent]=princomp(VoxN);
       Dots.Round.Var(i,:)=latent; %variance in three component axes.
       %Dots.Round.Long(i)=latent(1)/max(.1,latent(2)); %ratio of variances between the longest and the second longest axes, indication of roundness
       %Dots.Round.Oblong(i)=max(abs(Sc(:,1)))/max(1,abs(max(Sc(:,2)))); %max(abs(Sc(:,1))) and max(abs(Sc(:,2))) will the lengths from the center to be the furthest point in the principal and the second axes, respectively.
       Dots.Round.Long(i)=max(.1,latent(2))/latent(1); %HO 6/25/2010 flipped the following division to adjust to thresholding direction in anaSG.
       Dots.Round.Oblong(i)=max(1,abs(max(Sc(:,2))))/max(abs(Sc(:,1))); %HO 6/25/2010 flipped the following division to adjust to thresholding direction in anaSG.
       Dots.Round.SumVar(i)=sum(latent);
   else
       Dots.Round.Var(i,:)=[0;0;0];
       Dots.Round.Long(i)=0;
       Dots.Round.SumVar(i)=0;
       Dots.Round.Oblong(i)=1;
   end
   
  
   %% find surface area (Faces for a given voxel within a punctum will be 0
   %% to 6, the nubmer of voxels in 6-connectivity neighbors that are
   %% outside the punctum)
   for v = 1:size(Dots.Vox(i).Pos,1)
       Conn=dist(Dots.Vox(i).Pos, Dots.Vox(i).Pos(v,:));
       Dots.Vox(i).Faces(v)=6-sum(Conn==1);
   end
   %Dots.Round.histFaces(i,:)=hist(Dots.Vox(i).Faces,0:1:6); %histogram doesn't do anything HO 7/6/2010
   Dots.Round.meanFaces(i)=mean(Dots.Vox(i).Faces);
   %HO 6/25/2010 flipped the following division so that more spherical or 
   %smooth dots get higher values and ugly dots get lower values.
   %Thresholding in anaSG works well this way. Now Compact is definded as
   %the ratio of the average number of faces (facing outside) per voxel for
   %an ideal spherical dot with the same volume as a given dot over that
   %for the given dot. This value can be >1 or <1, doesn't tell much about
   %how close to sphere if the dots are small, but tells more about how
   %smooth the surface connectivity is. The actually roundness is given
   %better in Round.Long or Round.Oblong.
   %Dots.Round.Compact(i)=Dots.Round.meanFaces(i)/RoundFaces(Dots.Vol(i)); %This was original Compact definition. HO
   Dots.Round.Compact(i)=RoundFaces(Dots.Vol(i))/Dots.Round.meanFaces(i); 
end

%Instead of saving Dots, change the Dots back to Grouped
%and save Grouped. 7/6/2010 HO
%save([TPN 'Dots.mat'],'Dots')
Grouped = Dots;
save([TPN 'Grouped.mat'],'Grouped')

