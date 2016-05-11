function LDSDotsDD(TPN)

%7/8/2010 HO added comments.
%AllSegCut.mat was changed to AllSeg.mat 6/25/2010 HO
%load Settings to get xyum and zum instead of including in the input
%argument of the program. 7/30/2010 HO

%Added simple interpretation of arbors for bistratified (or more). Instead
%of taking FWHM, take whole arbors divided by valley. 10/13/2011 HO
%Also changed figure style and saving into tif instead of ai file.
%10/13/2011 HO
%Use SkelFiner.mat instead of AllSeg. 10/18/2011

'Starting analysis of Dots and Dendrites'


%% Get Path to Big Centroid

DPNd=[TPN 'data\'];
load([TPN 'Settings']); %7/30/2010 HO
ImageInfo = evalin('base', 'Settings');
xyum = ImageInfo.ImInfo.xyum; %7/30/2010 HO
zum = ImageInfo.ImInfo.zum; %7/30/2010 HO

if exist([TPN 'Dots.mat']) & exist([TPN 'SkelFiner.mat']) %Get results if info available %HO 6/25/2010
    Cname=TPN;
    load([TPN 'Dots.mat'])
    load([TPN 'SkelFiner.mat'])
    AllSegCut = cat(2, Skel.SegStats.Seg(:,2,:), Skel.SegStats.Seg(:,1,:), Skel.SegStats.Seg(:,3,:)); %HO 10/18/2011
    %load([TPN 'data\AllSeg.mat']) %HO 6/25/2010
    %load([TPN 'data\AllSegCut.mat']) %HO 6/25/2010
    %AllSegCut = AllSeg; %take advantage of already-written program HO 6/25/2010
    %load([TPN 'ShiftInfo.mat']) %7/8/2010 HO I don't have this!
    load([TPN 'Use.mat'])

    BC=Use.DPos; %Extract Dot Position (position in um)



    TotalDots=size(Use.DPos,1);
    save([DPNd 'TotalDots.mat'],'TotalDots')

    %% Assign Dend properties
    'Finding Dendrite Segment Properties'

    %%Find Segment Lengths
    SegLength=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 +(AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 +(AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);
    save([DPNd 'SegLength.mat'],'SegLength')
    TotalLength=sum(SegLength);
    save([DPNd 'TotalLength.mat'],'TotalLength')
    AverageDensity=TotalDots/TotalLength;
    save([DPNd 'AverageDensity'], 'AverageDensity')

    %%Find Segment Midpoints
    SegMP=(AllSegCut(:,:,1)+AllSegCut(:,:,2))/2; %seg midpoints in um

    %% Bin with Depth
    'Binning Data with Depth'

    %%Depth
    bin=zum;%zum;
    step=zum;
    bin2=7; %size to average %smoothing with 6 neighboring points
    filt=ones(bin2,1);

    maxDepth=fix(max(max(SegMP(:,3)), max(BC(:,3))))+1;
    minDepth=max(fix(min(min(SegMP(:,3)), min(BC(:,3)))),1); %max(.., 1) is for the case when min(SegMP(:,3)) or min(BC(:,3)) is zero. Put 1 instead of 0.
    minDepth = minDepth-mod(minDepth,zum); %Set min depth to factor of zum. If minDepth=1 and zum=0.3, mod(minDepth,zum) will be 0.1.
    DotBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);
    DendBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);

    b=0;
    for d=minDepth:step:maxDepth
        b=b+1;
        DotBinDepth(b)=sum(BC(:,3)>=(d-bin/2) & BC(:,3) <= (d+bin/2));
        DendBinDepth(b)=sum(SegLength(SegMP(:,3)>=(d-bin/2) & SegMP(:,3) <= (d+bin/2)));
    end %end d, run all depth bins

    %%filter to smooth
    DotBD=filter(filt,bin2,DotBinDepth);
    DendBD=filter(filt,bin2,DendBinDepth);


    DotPerDendDepth=DotBD./DendBD;

    cfigure(35,15);
    subplot(3,2,1);plot(DendBD,'r');
    title('Total dendritic length in each z bin, smoothed with 7 bins');
    subplot(3,2,3);plot(DotBD,'g');
    title('Total dot number in each z bin, smoothed with 7 bins');
    subplot(3,2,5);plot(DotPerDendDepth);
    title('Dot number / dend length in each z bin, not smoothed');
    pause(.01)
    %}

    %% Identify Strata
    PeakFind='Manual'
    BoarderFind='full width half maxima';
    %Have strata been previously identified
    clear P right left
    if exist([DPNd 'Results.mat']),
        load([DPNd 'Results.mat'])
        for a=1:size(Results.Arbor,2)
            left(a)=Results.Arbor(a).Top/zum;
            right(a)=Results.Arbor(a).Bottom/zum;
            P(a,1)=Results.Arbor(a).Peak/zum
        end
        Valley = Results.Arbor(1).Valley; %HO 10/15/2011
    else
        %%manual Identify Strat

%             if exist([DPNd 'Results.mat']),
%                 load([DPNd 'Results.mat'])
%                 Results.Arbor.Peak
%             end
        Valley = 0; %put dummy Valley HO 10/15/2011
        'Select peaks and click return'
        [x y]=ginput %choose two peaks for bistratified GCs.
        P=round(x); %so P is the z stack number closest to the peak of arbor distribution.
        peaks=DendBinDepth*0;
        if size(x,1)>1 %if bistratified, point valley between two peaks.
            'Identify Valley '
            [x y]=ginput;
            Valley=round(x);
        else, Valley=0; %stray value for Valley
        end

        DendBDmin=imregionalmin(DendBD);
        %%find peak edges using full width half maxima and reginal minimus
        for i =1:size(P,1)
            top=DendBD(P(i)); %top is the max arbor density in DendBD.
            left(i)=1; right(i)=size(DendBD,1);
            for l=1:P(i)-1
                if (DendBD(P(i)-l)< top/2) | sum((P(i)-l)==Valley) %pick the left side at the FWHM, or for bistratified, if the valley is higher than half max, pick the point at Valley.
                    left(i)=P(i)-l; break;   end %search for left side
            end %end searching for left side
            for r=1:size(DendBD,1)-P(i)
                if (DendBD(P(i)+r)< top/2) | sum((P(i)+r)==Valley) %pick the right side at the FWHM, or for bistratified, if the valley is higher than half max, pick the point at Valley.
                    right(i)=P(i)+r; break;   end %search for left side
            end %end searching for left side
        end %find peak width
    end % get edges of peaks

    %}

    %%find arbor stats
    APLength=DendBD*0;
    APDots=DendBD*0;
    APDotDend=DendBD*0;
    for i=1:size(P,1)
        ALength(i)=sum(DendBinDepth(left(i):right(i)));
        APLength(left(i):right(i))=ALength(i);
        ADots(i)=sum(DotBinDepth(left(i):right(i)));
        APDots(left(i):right(i))=ADots(i);
        ADotDend(i)=ADots(i)/ALength(i);
        APDotDend(left(i):right(i))=ADotDend(i);
    end %run all arbors

    %HO 10/13/2011
    %So, arbor stats ignores dendrites and dots that did not fall into
    %the FWHM. I added around line 200 simple interpretation of arbors
    %by just dividing into both sides of valley for bi RGC.


    %%Look at Valley


    if size(P,1)>1 %if bi
        vdm=DendBD*0;
        vr=2; %valley width = vr *2 +1
        for i=min(right):max(left) %run valley
            vdm(round(i))=mean(DotPerDendDepth(i-vr:i+vr)); %smoothing with 5 bins
        end
        vcent=find(vdm==min(vdm(min(right):max(left)))); %valley center, find where is the valley for Dot/Dend, not just Dend.

        VPLength=DendBD*0;
        VPDots=DendBD*0;
        VPDotDend=DendBD*0;

        %%Find edges of valley
        vleft=vcent-vr;
        vright=vcent+vr;
        %%Find stats
        VLength=sum(DendBinDepth(vleft:vright));
        VPLength(vleft:vright)=VLength;
        VDots=sum(DotBinDepth(vleft:vright));
        VPDots(vleft:vright)=VDots;
        VDotDend=VDots/VLength;
        VPDotDend(vleft:vright)=VDotDend;
    end %if bi


    %%Draw stats
    subplot(3,2,2);plot(APLength,'r');
    title('Dend length sum within arbors defined by FWHM');
    subplot(3,2,4);plot(APDots,'g');
    title('Dot number sum within arbors defined by FWHM');
    subplot(3,2,6);plot(APDotDend,'b');
    title('Dot/Dend within arbors defined by FWHM');
    %title(Cname)
    if exist('VPDotDend','var')
        hold on
        plot(VPDotDend,'b')
        hold off
    end
    %}


    %HO added simple ON OFF arbor stats (just diving into both sides
    %from the valley 10/13/2011
    if size(P,1)>1 %if bi
        ValleyList = [0, Valley, length(DendBinDepth)];
        APLengthSimple=DendBD*0;
        APDotsSimple=DendBD*0;
        APDotDendSimple=DendBD*0;
        for i=1:size(P,1)
            ALengthSimple(i)=sum(DendBinDepth(ValleyList(i)+1:ValleyList(i+1)));
            APLengthSimple(ValleyList(i)+1:ValleyList(i+1))=ALengthSimple(i);
            ADotsSimple(i)=sum(DotBinDepth(ValleyList(i)+1:ValleyList(i+1)));
            APDotsSimple(ValleyList(i)+1:ValleyList(i+1))=ADotsSimple(i);
            ADotDendSimple(i)=ADots(i)/ALength(i);
            APDotDendSimple(ValleyList(i)+1:ValleyList(i+1))=ADotDendSimple(i);
        end %run all arbors
        %add the simple arbor stats on the figure
        subplot(3,2,2);plot(APLength,'r');hold on;plot(APLengthSimple,'r:');
        title('Dend length sum within arbors defined by FWHM (solid) and by simple division (dotted)');
        subplot(3,2,4);plot(APDots,'g');hold on;plot(APDotsSimple,'r:');
        title('Dot number sum within arbors defined by FWHM (solid) and by simple division (dotted)');
        subplot(3,2,6);plot(APDotDend,'b');hold on;plot(VPDotDend,'b');plot(APDotDendSimple,'r:');
        title('Dot/Dend within arbors and valley defined by FWHM (solid) and by simple division (dotted)');
    end

    %HO for monostratified RGC, simple calculation of dend and dot is
    %to use whole dend and whole dots. 10/15/2011
    if size(P,1)==1 %if mono
        ALengthSimple = TotalLength;
        ADotsSimple = TotalDots;
        ADotDendSimple = ADots/ALength;
    end


    if isdir([TPN 'images'])==0, mkdir([TPN 'images']); end %create directory to store steps
    %saveas(gcf,[TPN 'images\DD'],'ai') %save figure as illustrator file in images
    saveas(gcf,[TPN 'images\DD'],'tif') %HO save figure as tif file in images 10/13/2011
    %save([DPNd 'DepthFigBin.mat'], 'bin2') %save size of figure bin

    pause(2)



    %% Enter Data to Structured Array
    if exist([TPN 'data' filesep 'Results.mat'])
        load([TPN 'data' filesep 'Results.mat'])
    end
    %Get information
    %Results.type=input('What is the Cell type?  ', 's')
    %Results.age=input('What is the Cell age? ')
    %Results.notes=input('Enter notes here ==> ', 's')


    Results.location=Cname;
    Back=find(Cname== filesep);
    Results.name=Cname(Back(size(Back,2)-1)+1:Back(size(Back,2))-1)
    Results.ImageInfo.xyum=xyum;
    Results.ImageInfo.zum=zum;
    Results.CellStats.TotalLength=TotalLength;
    Results.CellStats.TotalDots=TotalDots;
    Results.CellStats.AverageDensity=AverageDensity;
    Results.Depth.DotsBinDepth=DotBD;
    Results.Depth.DendBinDepth=DendBD;
    Results.Depth.DotPerDendDepth=DotPerDendDepth;
    Results.Depth.AverageingBinWidthInMicrons=bin2*zum;

    for i=1:size(P,1)
        Results.Arbor(i).PeakFind=PeakFind;
        Results.Arbor(i).BoarderFind=BoarderFind;
        Results.Arbor(i).Length=ALength(i);
        Results.Arbor(i).Dots=ADots(i);
        Results.Arbor(i).DotDend=ADotDend(i);
        Results.Arbor(i).TotalDendLengthSimple=ALengthSimple(i); %HO 10/13/2011
        Results.Arbor(i).TotalDotNumSimple=ADotsSimple(i); %HO 10/13/2011
        Results.Arbor(i).NumDotOverDendLengthSimple=ADotDendSimple(i); %HO 10/13/2011
        Results.Arbor(i).Top=left(i)*zum; %now top and bottom are converted to um from the first z stack
        Results.Arbor(i).Peak=P(i)*zum; %now top and bottom are converted to um from the first z stack
        Results.Arbor(i).Bottom=right(i)*zum;
        Results.Arbor(i).Valley = Valley; %for monostratified RGC, dummy Valley=0 will be registered. HO 10/15/2011
    end


    %save Results
    save([TPN 'data\Results.mat'],'Results')


    %{
    if exist('./Dat.mat')
        load('./Dat.mat')
        for i= 1:size(Dat,2)
            if strcmp(Cname,Dat(i).name); targ=i; break
            else targ=size(Dat,2)+1; end
        end
        Dat(targ)=Results;

        save('./Dat.mat','Dat')

    end %if Dat exists
        %}



        %% Draw Dot and Dend
        %{
    'Drawing Dots and Dendrites'

    Sc=(1/xyum)/2;
    DD=uint8(zeros(round(max(max(AllSegCut(:,1,:)))*Sc),round(max(max(AllSegCut(:,2,:)))),round(max(max(AllSegCut(:,3,:)))*Sc)));


    %%Draw Segments
    SkelRes=.1;
    for i=1:size(AllSegCut,1)
        Dist=sqrt((AllSegCut(i,1,1)-AllSegCut(i,1,2))^2 + (AllSegCut(i,2,1)-AllSegCut(i,2,2))^2 + (AllSegCut(i,3,1)-AllSegCut(i,3,2))^2); %find distance
        Length(i)=Dist;
          devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
        for d=1:devs+1
            sy=AllSegCut(i,1,1)+((AllSegCut(i,1,2)-AllSegCut(i,1,1))/devs)*(d-1);
            sx=AllSegCut(i,2,1)+((AllSegCut(i,2,2)-AllSegCut(i,2,1))/devs)*(d-1);
            sz=AllSegCut(i,3,1)+((AllSegCut(i,3,2)-AllSegCut(i,3,1))/devs)*(d-1);
            DD(round(sy*Sc)+1,round(sx*Sc)+1,round(sz*Sc)+1)=1; %draw Skel
        end
    end
    clear Dist

    %%DrawNodes
    for i=1:size(AllSegCut,1)
        DD(round(AllSegCut(i,1,1)*Sc)+1,round(AllSegCut(i,2,1)*Sc)+1,round(AllSegCut(i,3,1)*Sc)+1)=2;
        DD(round(AllSegCut(i,1,2)*Sc)+1,round(AllSegCut(i,2,2)*Sc)+1,round(AllSegCut(i,3,2)*Sc)+1)=2;
    end

    %%Draw Dots
    for i=1:size(BC,1)
        DD(round(BC(i,1)*Sc)+1,round(BC(i,2)*Sc)+1,round(BC(i,3)*Sc)+1)=3;
    end

    save([TPN 'pics\DD.mat'])

        %}

        %% Finish

        [TPN(size(TPN,2)-6:size(TPN,2)-1)]
        DotDendAt=uint16(clock)
        save([TPN 'data\DotDendAt.mat'],'DotDendAt') %record the date and time of this analysis
        'Done DotDend'
else
        'Data Not Available'
end

