function [DepthDev] = LDSStratNoCBmedian(TPN)
%
%7/8/2010 HO added comments
%This program calculate the standard deviations of the z position of dots
%and arbors found within the 30um xy radius circle from each 1um*1um xy
%pixel in the territory. So, the number of standard deviations you can get
%is equal to the number of 1um*1um xy pixel in the territory. Then, the
%average starndard deviation is estimated by simply averanging all the
%standard deviations in the territory.
%
% 
% clear all
% 
% TPN = GetMyDir;
% yxum=0.103;
% zum=0.3;
% 


LookRad=30;

   load([TPN 'CA.mat']) 
   Territory=CA.Arbor(1).Territory;
   
   load([TPN 'Use.mat'])
   DPos=Use.DPos;
   
   load([TPN 'data' filesep 'AllSeg.mat']) %7/28/2010 HO
   AllSegCut = AllSeg; clear AllSeg; %7/28/2010 HO
   Mids=mean(AllSegCut,3);
   
   %Remove area near cell body
   Dist=sqrt((Mids(:,1)-Use.Cent(1)).^2 + (Mids(:,2)-Use.Cent(2)).^2);
   DistS=sort(Dist);
   Inner10=DistS(fix(size(DistS,1)/10)); %assuming that the closest 10% of arbors to cell body is pretty much in cell body, but I think this is overestimate.
   CutPrims=max(Inner10,10); %take 10um if the closest 10% is still <10um.
   Mids=Mids(Dist>CutPrims,:);
   Dist=sqrt((DPos(:,1)-Use.Cent(1)).^2 + (DPos(:,2)-Use.Cent(2)).^2);
   DPos=DPos(Dist>CutPrims,:);
   
   
   clear DepthDev
   'Find dots over Area'
   [y x] = find(Territory);
   NumDots=zeros(size(y,1),1);
   SDevDots=NumDots;
   for i = 1: size(y,1)
      Dist=sqrt((DPos(:,1)-y(i)).^2 + (DPos(:,2)-x(i)).^2);
      Near=Dist<LookRad;
      Nears=DPos(Near,:);
      NumDots(i)=size(Nears,1);
      SDevsDots(i)=std(Nears(:,3));
   end
   
   'Find mids over area'
   NumMids=zeros(size(y,1),1);
   SDevsMids=NumMids;
   for i = 1: size(y,1)
      Dist=sqrt((Mids(:,1)-y(i)).^2 + (Mids(:,2)-x(i)).^2);
      Near=Dist<LookRad;
      Nears=Mids(Near,:);
      NumMids(i)=size(Nears,1);
      SDevsMids(i)=std(Nears(:,3));
   end
     
   OKs=(NumDots>10) & (NumMids > 40);
   DepthDev.DotsA=mean(SDevsDots(OKs));
   DepthDev.DendA=mean(SDevsMids(OKs));
   DepthDev;
   
   
   
   'Find Each Dot'
   %Run from Dot perspective
   %Not clear what you can gain from this HO 7/30/2010
   NumDots=zeros(size(DPos,1),1);
   DtoDots=NumDots;
   for i = 1:size(DPos,1) 
      Dist=sqrt((DPos(:,1)-DPos(i,1)).^2 + (DPos(:,2)-DPos(i,2)).^2);
      Near=Dist<LookRad;
      Nears=DPos(Near,:);
      NumDots(i)=size(Nears,1);
      meanDepth=median(Nears(:,3));
      DtoDots(i)=abs(meanDepth-DPos(i,3));
   end
     
   'Find Each Mid'
      %Run from Dot perspective
      %actually this is from Mids perspective! HO 7/30/2010
      %Not clear what you can gain from this HO 7/30/2010
   NumMids=zeros(size(Mids,1),1);
   DtoMids=NumMids;
   for i = 1:size(Mids,1)
      Dist=sqrt((Mids(:,1)-Mids(i,1)).^2 + (Mids(:,2)-Mids(i,2)).^2); %Mids perspective.
      Near=Dist<LookRad;
      Nears=Mids(Near,:);
      NumMids(i)=size(Nears,1);
      meanDepth=median(Nears(:,3));
      DtoMids(i)=abs(meanDepth-Mids(i,3));
   end

   
   

   DepthDev.DtoDots=mean(DtoDots(NumDots>10));
   DepthDev.DtoMids=mean(DtoMids(NumMids>40));
   save([TPN 'data' filesep 'DepthDevNoCBmedian.mat'],'DepthDev')
   DepthDev   
   
   
   %% Draw Data
%  
%    for i = 1: size(Mids,1)
%       DrawMids(fix(Mids(i,1))+1,fix(Mids(i,2))+1,fix(Mids(i,3))+1)=1;
%    end
%    imwriteNp(TPN,DrawMids,'DrawMids')
%    %}
       
