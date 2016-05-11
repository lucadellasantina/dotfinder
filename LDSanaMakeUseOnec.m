function[] = LDSanaMakeUseOnec(TPN)

% 1/15/2010 HO added comments to make it easier to understand each line.
% Plus, changed ImInfo{4} and {5} to my style, which is ImInfo.xyum and
% zum.
% 6/25/2010 HO load AllSeg.mat instead of AllSegCut.mat. We have no longer
% AllSegCut.mat. See HOanaMa.
% 10/18/2011 HO don't use AllSeg anymore. Instead use SkelFiner.mat.


load([TPN 'Settings.mat'])
ImageInfo = evalin('base', 'Settings');
ImageInfo = ImageInfo.ImInfo;
xyum=ImageInfo.xyum; %changed to the structure style 1/15/2010
zum=ImageInfo.zum; %changed to the structure style 1/15/2010


    if exist([TPN 'find\SG.mat'])
        load([TPN 'Dots.mat'])
        load([TPN 'find\SG.mat'])
%         load([TPN 'data\AllSeg.mat']) %HO 6/25/2010
%         %load([TPN 'data\AllSegCut.mat']) %HO 6/25/2010
%         AllSegCut = AllSeg; %take advantage of already-written program HO 6/25/2010
        load([TPN 'Grouped.mat'])
        load([TPN 'SkelFiner.mat']) %HO 6/25/2010
        
        AllSegCut = cat(2, Skel.SegStats.Seg(:,2,:), Skel.SegStats.Seg(:,1,:), Skel.SegStats.Seg(:,3,:)); %HO 10/18/2011
        
        if exist([TPN 'Cell.mat'])
        load([TPN 'Cell.mat'])
        cName=['P' Cell.Age '_' Cell.Type '_' Cell.Name '_' ]
        Use.cName=cName;
        clear Cell cName
        end
        
        DPos=round(Grouped.Pos); %Grouped.DPos to Grouped.Pos 7/7/2010 HO
        DPos(:,1:2)=DPos(:,1:2)*xyum; DPos(:,3)=DPos(:,3)*zum;
        Cent=Dots.Im.CBpos;
        Cent=Cent*xyum;
        Use.DPos=DPos;
        Use.Cent=Cent;
        

        %%Extract Dend positions 
        Mids=mean(AllSegCut,3); %segment xyz position calculated as mean of two node positions
        Length=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 ...
           + (AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 ...
           + (AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);
        
       Use.Mids=Mids;
       Use.Length=Length;
       
       %HO added recording of stratification if you have
       %SratificationIndMap.mat 10/15/2011
       if isfield(Skel.FilStats, 'SkelStratification');
            Edges=Skel.FilStats.aEdges+1;
            SegStrat=[Skel.FilStats.SkelStratification(Edges(:,1)); Skel.FilStats.SkelStratification(Edges(:,2))];
            SegStrat=SegStrat';
            MidsStrat = Skel.FilStats.EdgeStratification;
            Use.SegStrat = SegStrat;
            Use.MidsStrat = MidsStrat;
       end
       
       
       Nearest = zeros(size(DPos, 1),1);
        for i = 1:size(DPos,1)
            Ndist=dist(Mids,DPos(i,:)); %find dist from dot to all nodes
            Near=min(Ndist); %find shortest distance
            Nearest(i)=find(Ndist==Near,1); %get node at that distance
            %DotToNN(i,:)=Near; %record that distance for posterity %not
            %used, skip this line HO
        end
         NN=Mids(Nearest,:); %assign that node to NearestNode list for dots
         Use.NN=NN;
        
        %HO added recording of stratification if you have
        %StratificationIndMap.mat 10/15/2011
        if isfield(Skel.FilStats, 'SkelStratification');
            Use.NNStrat = MidsStrat(Nearest);
        end
       
        clear NN OK DPos Cent Dots SG
        clear Mids Length AllSgeCut
       
       %%Extract Depth restrictions
        if exist([TPN 'data\Results.mat'])
        load([TPN 'data\Results.mat'])
        clear Top Bottom
        for i = 1: size(Results.Arbor,2)
            Top(i)=Results.Arbor(i).Top;
            Bottom(i)=Results.Arbor(i).Bottom;
        end
        clear Results
        
        Use.Top=Top;
        Use.Bottom=Bottom;
        clear Top Bottom
        end
       
       
        
       save([TPN 'Use.mat'],'Use') 
        
    end
        
