function [Grad] = LDSGradient(TPN)

%7/8/2010 HO added comments.
%7/30/2010 HO changed naming of all the parameters to make it intuitive.
%Also, put InnerLimit = 10, which is the same as before but more organized,
%and Middle was set between outer limit and inner limit, NOT between outer
%limit and 0 like before. This would not change the values much.

%%Get Numbers for Inner Half, Outer Half, Inner Third, Outer Third,
%%inner half - CB, Outer Half / Inner Half - CB

%HO 10/18/2011 added loop for bistratified RGCs.

% Run as vectors
load([TPN 'Use.mat'])
% load([TPN 'data\AllSeg.mat']) %7/28/2010 HO
% AllSegCut = AllSeg; clear AllSeg; %7/28/2010 HO

%Why repeating this? Mids and Length are under Use.mat. 10/15/2011 HO
Mids=Use.Mids;
Length=Use.Length;
% Mids=mean(AllSegCut,3);
% Length=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 ...
%     + (AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 ...
%     + (AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);

%Find eccentricities (not direct 3D distance from cell body, but distance
%on the z-projected image)
EccDPos=sqrt((Use.DPos(:,1)-Use.Cent(1)).^2 + (Use.DPos(:,2)-Use.Cent(2)).^2);
EccMids=sqrt((Mids(:,1)-Use.Cent(1)).^2 + (Mids(:,2)-Use.Cent(2)).^2);

%HO 10/18/2011 added stratification stuff for bistratified
if isfield(Use, 'MidsStrat');
    MidsStrat = Use.MidsStrat;
    MidsStrat = MidsStrat';
end
if isfield(Use, 'NNStrat');
    NNStrat = Use.NNStrat;
    NNStrat = NNStrat';
end
if exist([TPN, 'Strat.mat']);
    load([TPN, 'Strat.mat']);
end

load([TPN 'CA.mat'])


%% Run entire arbor first HO added this part 10/18/2011
'Running entire arbor'
for a = 1:length(CA.Arbor);
    %%find edges
    Outer=max(EccMids);
    Inner=0;
    Middle = (Outer+Inner)/2;
    GradAll(a).AnalysisOuterLimit=Outer;
    GradAll(a).AnalysisInnerLimit=Inner;
    %This part is the straightforward calculatin of P/D (without inserting area
    %in between).
    
    if a==1; %entire cell
        InnerDendLength = sum(Length((EccMids>Inner) & (EccMids<=Middle)));
        OuterDendLength = sum(Length((EccMids>Middle) & (EccMids<=Outer)));
        InnerDots = sum((EccDPos>Inner) & (EccDPos<=Middle));
        OuterDots = sum((EccDPos>Middle) & (EccDPos<=Outer));
    else %bi
        'bistratified arbor'
        a-1
        
        InnerDendLength = sum(Length((EccMids>Inner) & (EccMids<=Middle) &...
            (MidsStrat>=Strat(a).DendVolOuterLimit) & (MidsStrat<=Strat(a).DendVolInnerLimit)));
        OuterDendLength = sum(Length((EccMids>Middle) & (EccMids<=Outer) &...
            (MidsStrat>=Strat(a).DendVolOuterLimit) & (MidsStrat<=Strat(a).DendVolInnerLimit)));
        InnerDots = sum((EccDPos>Inner) & (EccDPos<=Middle) &...
            (NNStrat>=Strat(a).DendVolOuterLimit) & (NNStrat<=Strat(a).DendVolInnerLimit));
        OuterDots = sum((EccDPos>Middle) & (EccDPos<=Outer) &...
            (NNStrat>=Strat(a).DendVolOuterLimit) & (NNStrat<=Strat(a).DendVolInnerLimit));
    end
    
    InnerPD = InnerDots/InnerDendLength;
    OuterPD = OuterDots/OuterDendLength;
    PDInnerOuterRatioStraightCalc = InnerPD/OuterPD %Unlike P/D in CAsampleUseNN2007, this is a straighforward calculation of P/D and the inner-outer ratio of that.
 
    GradAll(a).InnerDotStraightCalc = InnerDots;
    GradAll(a).OuterDotStraightCalc = OuterDots;
    GradAll(a).InnerDendStraightCalc = InnerDendLength;
    GradAll(a).OuterDendStraightCalc = OuterDendLength;
    GradAll(a).InnerPDStraightCalc = InnerPD;
    GradAll(a).OuterPDStraightCalc = OuterPD;
    GradAll(a).PDInnerOuterRatioStraightCalc=PDInnerOuterRatioStraightCalc; 



    %% Run from areas
    if a==1; %first round
        if size(CA.Arbor,2)>2 %this will be for multi-stratified GC
            Area=0;
            for i=2:length(CA.Arbor);
                Area=Area+CA.Arbor(i).Territory;
            end
        else %monostratified RGC
            Area=CA.Arbor(1).Territory;
        end
    else %second and third rounds for bi RGC
        Area=CA.Arbor(a).Territory;
    end

    %remember Area (CA.Arbor.Territory) has pixels of 1um-by-1um (see CAsampleUseNN).
    DistMap=zeros(size(Area));
    for y = 1:size(Area,1)
        for x = 1: size(Area,2)
            DistMap(y,x)=sqrt((y-Use.Cent(1))^2+(x-Use.Cent(2))^2);
        end
    end


% Chart all
    Bin=10;
    clear cDist cDots cDend cArea CAreaAbs
    if a==1;
        for i = Inner+Bin/2:Outer-Bin/2
            cDist(i)=i;
            cDots(i)=sum(EccDPos>(i-Bin/2) & EccDPos<=(i+Bin/2));
            cDend(i)=sum(Length(EccMids>(i-Bin/2) & EccMids<=(i+Bin/2)));
            cArea(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2)));
            cAreaAbs(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2))>0);
        end
    else
        for i = Inner+Bin/2:Outer-Bin/2
            cDist(i)=i;
            cDots(i)=sum(EccDPos>(i-Bin/2) & EccDPos<=(i+Bin/2) &...
                NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
            cDend(i)=sum(Length(EccMids>(i-Bin/2) & EccMids<=(i+Bin/2) &...
                MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
            cArea(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2)));
            cAreaAbs(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2))>0);
        end
    end

    cPD=cDots./cDend;
    cPA=cDots./cArea;
    cDA=cDend./cArea;
    
    LineColorList = ['bbrgyk'];
    if a<3;
        cfigure(20,20);
        subplot(3,1,1),hold on,plot(cPA, LineColorList(a));
        title('Entire dendrites');
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/A (puncta/um2)');
        subplot(3,1,2),hold on,plot(cDA, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('D/A (um/um2)');
        subplot(3,1,3),hold on,plot(cPD, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/D (puncta/um)');
        %ylim([0 .3]);
        %This part is not working. Omit. HO 10/18/2011
%         hold on
%         %Get Slope
%         [p, S]=polyfit(cDist,cPD,1);
%         Line=p(1)*cDist.^2+p(2);
%         plot(Line, LineColorList(a));
%         %ylim([0 .3]);
%         hold off
        pause(3)
    else
        subplot(3,1,1),hold on,plot(cPA, LineColorList(a));
        title('Entire dendrites, blue=Off arbor, red=On arbor');
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/A (puncta/um2)');
        subplot(3,1,2),hold on,plot(cDA, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('D/A (um/um2)');
        subplot(3,1,3),hold on,plot(cPD, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/D (puncta/um)');
        %ylim([0 .3]);
        %This part is not working. Omit. HO 10/18/2011
%         hold on
%         %Get Slope
%         [p, S]=polyfit(cDist,cPD,1);
%         Line=p(1)*cDist.^2+p(2);
%         plot(Line, LineColorList(a));
%         %ylim([0 .3]);
%         hold off
        pause(3)
    end
    
    if a==1;
        Name=[TPN 'images\PA&DA&PDvsEccentricityAllDend.tif'];
    else
        Name=[TPN 'images\PA&DA&PDvsEccentricityAllDend_Bi.tif'];
    end
    saveas(gcf,Name) %save figure with title

    GradAll(a).EccBin = Bin;
    GradAll(a).Eccentricity = cDist;
    GradAll(a).PvsEcc = cDots;
    GradAll(a).DvsEcc = cDend;
    GradAll(a).AvsEcc = cArea;
    GradAll(a).PAvsEcc = cPA;
    GradAll(a).DAvsEcc = cDA;
    GradAll(a).PDvsEcc = cPD;


    IArea=sum(Area(DistMap>Inner & DistMap<=Middle));
    OArea=sum(Area(DistMap>Middle & DistMap<=Outer));
    if a==1;
        IDot=sum(EccDPos>Inner & EccDPos<=Middle);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer));
    else
        IDot=sum(EccDPos>Inner & EccDPos<=Middle &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
    end
    I2O=OE-IE; %mean distance from inner to outer.
    IPA=IDot/IArea;
    OPA=ODot/OArea;
    IDA=IDend/IArea;
    ODA=ODend/OArea;
    IPD=IDot/IDend;
    OPD=ODot/ODend;

    dPA=(IPA-OPA)/((IPA+OPA)); %~~~~!!!!!!!!!! Scale?
    dDA=(IDA-ODA)/((IDA+ODA));
    dPD=(IPD-OPD)/((IPD+OPD)); 

    GradAll(a).InnerArea=IArea;
    GradAll(a).OuterArea=OArea;
    GradAll(a).InnerDot=IDot;
    GradAll(a).OuterDot=ODot;
    GradAll(a).InnerDend=IDend;
    GradAll(a).OuterDend=ODend;
    GradAll(a).InnerPA=IPA;
    GradAll(a).OuterPA=OPA;
    GradAll(a).InnerDA=IDA;
    GradAll(a).OuterDA=ODA;
    GradAll(a).InnerPD=IPD;
    GradAll(a).OuterPD=OPD;
    GradAll(a).PAInnerOuterRatio = IPA/OPA;
    GradAll(a).DAInnerOuterRatio = IDA/ODA;
    GradAll(a).PDInnerOuterRatio = IPD/OPD;

    PDInnerOuterRatio = IPD/OPD %just for displaying reasion to compare with straightforward calculation

    GradAll(a).dPA=dPA;
    GradAll(a).dDA=dDA;
    GradAll(a).dPD=dPD;

    GradAll(a).InnerMeanDendEcc = IE;
    GradAll(a).OuterMeanDendEcc = OE;
    GradAll(a).InnerOuterEccDiff=I2O;

    save([TPN 'GradAll.mat'],'GradAll')

end
    
%% Josh's way, remove innermost 10um radius and 2% of skels (Mids) that lie outermost (not simply 2% outermost radius)
'Excluding innermost 10um radius and outermost 2% of skels'
for a = 1:length(CA.Arbor);
    %%find edges
    SortMids=sort(EccMids(EccMids>10));
    Outer=SortMids(fix(size(SortMids,1)*.98)); %remove the outer most 2% (probably trying to remove outliers, for example if only one process goes far away, it will have a large impact on the decision of inner vs outer)
    Inner=10; %remove 10um around the cell body center because it is pretty much cell body, and almost no puncta on the cell body.
    Middle = (Outer+Inner)/2;
    Grad(a).AnalysisOuterLimit=Outer;
    Grad(a).AnalysisInnerLimit=Inner;
    %This part is the straightforward calculatin of P/D (without inserting area
    %in between).
    
    if a==1; %entire cell
        InnerDendLength = sum(Length((EccMids>Inner) & (EccMids<=Middle)));
        OuterDendLength = sum(Length((EccMids>Middle) & (EccMids<=Outer)));
        InnerDots = sum((EccDPos>Inner) & (EccDPos<=Middle));
        OuterDots = sum((EccDPos>Middle) & (EccDPos<=Outer));
    else %bi
        'bistratified arbor'
        a-1
        
        InnerDendLength = sum(Length((EccMids>Inner) & (EccMids<=Middle) &...
            (MidsStrat>=Strat(a).DendVolOuterLimit) & (MidsStrat<=Strat(a).DendVolInnerLimit)));
        OuterDendLength = sum(Length((EccMids>Middle) & (EccMids<=Outer) &...
            (MidsStrat>=Strat(a).DendVolOuterLimit) & (MidsStrat<=Strat(a).DendVolInnerLimit)));
        InnerDots = sum((EccDPos>Inner) & (EccDPos<=Middle) &...
            (NNStrat>=Strat(a).DendVolOuterLimit) & (NNStrat<=Strat(a).DendVolInnerLimit));
        OuterDots = sum((EccDPos>Middle) & (EccDPos<=Outer) &...
            (NNStrat>=Strat(a).DendVolOuterLimit) & (NNStrat<=Strat(a).DendVolInnerLimit));
    end
    
    InnerPD = InnerDots/InnerDendLength;
    OuterPD = OuterDots/OuterDendLength;
    PDInnerOuterRatioStraightCalc = InnerPD/OuterPD %Unlike P/D in CAsampleUseNN2007, this is a straighforward calculation of P/D and the inner-outer ratio of that.
 
    Grad(a).InnerDotStraightCalc = InnerDots;
    Grad(a).OuterDotStraightCalc = OuterDots;
    Grad(a).InnerDendStraightCalc = InnerDendLength;
    Grad(a).OuterDendStraightCalc = OuterDendLength;
    Grad(a).InnerPDStraightCalc = InnerPD;
    Grad(a).OuterPDStraightCalc = OuterPD;
    Grad(a).PDInnerOuterRatioStraightCalc=PDInnerOuterRatioStraightCalc; 



    %% Run from areas
    if a==1; %first round
        if size(CA.Arbor,2)>2 %this will be for multi-stratified GC
            Area=0;
            for i=2:length(CA.Arbor);
                Area=Area+CA.Arbor(i).Territory;
            end
        else %monostratified RGC
            Area=CA.Arbor(1).Territory;
        end
    else %second and third rounds for bi RGC
        Area=CA.Arbor(a).Territory;
    end

    %remember Area (CA.Arbor.Territory) has pixels of 1um-by-1um (see CAsampleUseNN).
    DistMap=zeros(size(Area));
    for y = 1:size(Area,1)
        for x = 1: size(Area,2)
            DistMap(y,x)=sqrt((y-Use.Cent(1))^2+(x-Use.Cent(2))^2);
        end
    end


%% Chart all
    Bin=10;
    clear cDist cDots cDend cArea CAreaAbs
    if a==1;
        for i = Inner:Outer
            cDist(i)=i;
            cDots(i)=sum(EccDPos>(i-Bin/2) & EccDPos<=(i+Bin/2));
            cDend(i)=sum(Length(EccMids>(i-Bin/2) & EccMids<=(i+Bin/2)));
            cArea(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2)));
            cAreaAbs(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2))>0);
        end
    else
        for i = Inner:Outer
            cDist(i)=i;
            cDots(i)=sum(EccDPos>(i-Bin/2) & EccDPos<=(i+Bin/2) &...
                NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
            cDend(i)=sum(Length(EccMids>(i-Bin/2) & EccMids<=(i+Bin/2) &...
                MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
            cArea(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2)));
            cAreaAbs(i)=sum(Area(DistMap>(i-Bin/2) & DistMap<=(i+Bin/2))>0);
        end
    end

    cPD=cDots./cDend;
    cPA=cDots./cArea;
    cDA=cDend./cArea;
    
    LineColorList = ['bbrgyk'];
    if a<3;
        cfigure(20,20);
        subplot(3,1,1),hold on,plot(cPA, LineColorList(a));
        title('Innermost 10um radius and outermost 2% dend skel removed');
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/A (puncta/um2)');
        subplot(3,1,2),hold on,plot(cDA, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('D/A (um/um2)');
        subplot(3,1,3),hold on,plot(cPD, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/D (puncta/um)');
        %ylim([0 .3]);
        %This part is not working. Omit. HO 10/18/2011
%         hold on
%         %Get Slope
%         [p, S]=polyfit(cDist,cPD,1);
%         Line=p(1)*cDist.^2+p(2);
%         plot(Line, LineColorList(a));
%         %ylim([0 .3]);
%         hold off
        pause(3)
    else
        subplot(3,1,1),hold on,plot(cPA, LineColorList(a));
        title('Innermost 10um radius and outermost 2% dend skel removed, blue=Off arbor, red=On arbor');
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/A (puncta/um2)');
        subplot(3,1,2),hold on,plot(cDA, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('D/A (um/um2)');
        subplot(3,1,3),hold on,plot(cPD, LineColorList(a));
        xlabel('eccentricity (2D distance from cell body)');ylabel('P/D (puncta/um)');
        %ylim([0 .3]);
        %This part is not working. Omit. HO 10/18/2011
%         hold on
%         %Get Slope
%         [p, S]=polyfit(cDist,cPD,1);
%         Line=p(1)*cDist.^2+p(2);
%         plot(Line, LineColorList(a));
%         %ylim([0 .3]);
%         hold off
        pause(3)
    end
    
    if a==1;
        Name=[TPN 'images\PA&DA&PDvsEccentricity.tif'];
    else
        Name=[TPN 'images\PA&DA&PDvsEccentricity_Bi.tif'];
    end
    saveas(gcf,Name) %save figure with title

    Grad(a).EccBin = Bin;
    Grad(a).Eccentricity = cDist;
    Grad(a).PvsEcc = cDots;
    Grad(a).DvsEcc = cDend;
    Grad(a).AvsEcc = cArea;
    Grad(a).PAvsEcc = cPA;
    Grad(a).DAvsEcc = cDA;
    Grad(a).PDvsEcc = cPD;


    IArea=sum(Area(DistMap>Inner & DistMap<=Middle));
    OArea=sum(Area(DistMap>Middle & DistMap<=Outer));
    if a==1;
        IDot=sum(EccDPos>Inner & EccDPos<=Middle);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer));
    else
        IDot=sum(EccDPos>Inner & EccDPos<=Middle &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
    end
    I2O=OE-IE; %mean distance from inner to outer.
    IPA=IDot/IArea;
    OPA=ODot/OArea;
    IDA=IDend/IArea;
    ODA=ODend/OArea;
    IPD=IDot/IDend;
    OPD=ODot/ODend;

    dPA=(IPA-OPA)/((IPA+OPA)); %~~~~!!!!!!!!!! Scale?
    dDA=(IDA-ODA)/((IDA+ODA));
    dPD=(IPD-OPD)/((IPD+OPD)); 

    Grad(a).InnerArea=IArea;
    Grad(a).OuterArea=OArea;
    Grad(a).InnerDot=IDot;
    Grad(a).OuterDot=ODot;
    Grad(a).InnerDend=IDend;
    Grad(a).OuterDend=ODend;
    Grad(a).InnerPA=IPA;
    Grad(a).OuterPA=OPA;
    Grad(a).InnerDA=IDA;
    Grad(a).OuterDA=ODA;
    Grad(a).InnerPD=IPD;
    Grad(a).OuterPD=OPD;
    Grad(a).PAInnerOuterRatio = IPA/OPA;
    Grad(a).DAInnerOuterRatio = IDA/ODA;
    Grad(a).PDInnerOuterRatio = IPD/OPD;

    PDInnerOuterRatio = IPD/OPD %just for displaying reasion to compare with straightforward calculation

    Grad(a).dPA=dPA;
    Grad(a).dDA=dDA;
    Grad(a).dPD=dPD;

    Grad(a).InnerMeanDendEcc = IE;
    Grad(a).OuterMeanDendEcc = OE;
    Grad(a).InnerOuterEccDiff=I2O;

    save([TPN 'Grad.mat'],'Grad') %used to be saved as GradA.mat 7/30/2010 HO
end

%% Run inner only <40um
for a = 1:length(CA.Arbor);
    Outer=40;
    Inner=10;
    Middle = (Outer+Inner)/2;

    IArea=sum(Area(DistMap>Inner & DistMap<=Middle));
    OArea=sum(Area(DistMap>Middle & DistMap<=Outer));
    if a==1;
        IDot=sum(EccDPos>Inner & EccDPos<=Middle);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer));
    else
        IDot=sum(EccDPos>Inner & EccDPos<=Middle &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        ODot=sum(EccDPos>Middle & EccDPos<=Outer &...
            NNStrat>=Strat(a).DendVolOuterLimit & NNStrat<=Strat(a).DendVolInnerLimit);
        IDend=sum(Length(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        ODend=sum(Length(EccMids>Middle & EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        IE=mean(EccMids(EccMids>Inner & EccMids<=Middle &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
        OE=mean(EccMids(EccMids>Middle& EccMids<=Outer &...
            MidsStrat>=Strat(a).DendVolOuterLimit & MidsStrat<=Strat(a).DendVolInnerLimit));
    end
    I2O=OE-IE; %mean distance from inner to outer.
    IPA=IDot/IArea;
    OPA=ODot/OArea;
    IDA=IDend/IArea;
    ODA=ODend/OArea;
    IPD=IDot/IDend;
    OPD=ODot/ODend;

    dPA=(IPA-OPA)/((IPA+OPA)); %~~~~!!!!!!!!!! Scale?
    dDA=(IDA-ODA)/((IDA+ODA));
    dPD=(IPD-OPD)/((IPD+OPD)); 

    GradI(a).AnalysisOuterLimit = Outer;
    GradI(a).AnalysisInnerLimit = Inner;
    GradI(a).InnerArea=IArea;
    GradI(a).OuterArea=OArea;
    GradI(a).InnerDot=IDot;
    GradI(a).OuterDot=ODot;
    GradI(a).InnerDend=IDend;
    GradI(a).OuterDend=ODend;
    GradI(a).InnerPA=IPA;
    GradI(a).OuterPA=OPA;
    GradI(a).InnerDA=IDA;
    GradI(a).OuterDA=ODA;
    GradI(a).InnerPD=IPD;
    GradI(a).OuterPD=OPD;
    GradI(a).PAInnerOuterRatio = IPA/OPA;
    GradI(a).DAInnerOuterRatio = IDA/ODA;
    GradI(a).PDInnerOuterRatio = IPD/OPD;

    GradI(a).dPA=dPA;
    GradI(a).dDA=dDA;
    GradI(a).dPD=dPD;

    GradI(a).InnerMeanDendEcc = IE;
    GradI(a).OuterMeanDendEcc = OE;
    GradI(a).InnerOuterEccDiff=I2O;
    
    save([TPN 'GradI.mat'],'GradI')

end


%% HO added arbor comparison of these parameters 10/18/2011
%since ectopic OFF arbor of G10 happend in the peripheral, the fair
%comparison or fair estimate of P/D will be to compare it with ON layer
%arbor around the same eccentricity because PSD expression tends to go down
%in the peripheral. So, use dend length of OFF arbor at each eccentricity
%as weighting factor for ON arbor
if length(GradAll)>2; %if bistratified RGC, make effective comparison between different arbors
    OffDendEccInnerLimit = find(GradAll(2).DvsEcc>0, 1, 'first');
    OffDendEccOuterLimit = find(GradAll(2).DvsEcc>0, 1, 'last');
    EccRange = OffDendEccOuterLimit - OffDendEccInnerLimit + 1;
    EccBin = GradAll(2).EccBin;
    EccNotFitting = mod(EccRange, EccBin);
    if EccNotFitting == 0;
        NumBin = (EccRange-EccNotFitting)/EccBin;
        DendAnalysisInnerStartEcc = OffDendEccInnerLimit+EccBin/2;
    else
        NumBin = (EccRange-EccNotFitting)/EccBin + 1;
        DendAnalysisInnerStartEcc = OffDendEccInnerLimit + round(EccNotFitting/2);
    end
    
    Ecc = DendAnalysisInnerStartEcc:EccBin:DendAnalysisInnerStartEcc+(NumBin-1)*EccBin;
    for i=1:NumBin;
        CorrectEccInd = find(GradAll(2).Eccentricity==Ecc(i));
        OffP(i) = GradAll(2).PvsEcc(CorrectEccInd);
        OffD(i) = GradAll(2).DvsEcc(CorrectEccInd);
        OffA(i) = GradAll(2).AvsEcc(CorrectEccInd);
        OnP(i) = GradAll(3).PvsEcc(CorrectEccInd);
        OnD(i) = GradAll(3).DvsEcc(CorrectEccInd);
        OnA(i) = GradAll(3).AvsEcc(CorrectEccInd);
    end
    
    OffPTotal = sum(OffP); OffDTotal = sum(OffD); OffATotal = sum(OffA);
    OnPTotal = sum(OnP); OnDTotal = sum(OnD); OnATotal = sum(OnA);
    OffPoverDTotal = OffPTotal/OffDTotal; OffPoverATotal = OffPTotal/OffATotal; OffDoverATotal = OffDTotal/OffATotal;
    OnPoverDTotal = OnPTotal/OnDTotal; OnPoverATotal = OnPTotal/OnATotal; OnDoverATotal = OnDTotal/OnATotal;
    
    %the other way to calculated this
    if find(OffD == 0 | OnD == 0 | OffA == 0 | OnA == 0) %remove eccentricity where OFF or ON dendrite or area don't exist
        ZeroDorAInd = find(OffD == 0 | OnD == 0 | OffA == 0 | OnA == 0);
        OffD(ZeroDorAInd) = [];
        OffP(ZeroDorAInd) = [];
        OffA(ZeroDorAInd) = [];
        OnD(ZeroDorAInd) = [];
        OnP(ZeroDorAInd) = [];
        OnA(ZeroDorAInd) = [];
    end
        
    OffPD = OffP./OffD; OffPA = OffP./OffA; OffDA = OffD./OffA;
    OnPD = OnP./OnD; OnPA = OnP./OnA; OnDA = OnD./OnA;
    %OffPoverDTotal2 = sum(OffPD.*OffD)/sum(OffD); %this should be the same as OffPoverDTotal
    
    %now replace OffPD with OnPD. So if the On arbor has the same
    %dend length as Off arbor at each eccentricity but different P/D at
    %each eccentricity, what is the total P/D? The number will become much
    %lower than OnPoverDTotal if there are more inner dendrites with higher
    %P/D on the ON arbor.
    OnPoverDTotalWeightedByOffD = sum(OnPD.*OffD)/sum(OffD);
    
    %try same thing for PA and DA
    %OffPoverATotal2 = sum(OffPA.*OffA)/sum(OffA); %this should be the same as OffPoverATotal
    OnPoverATotalWeightedByOffA = sum(OnPA.*OffA)/sum(OffA);
    
    %OffDoverATotal2 = sum(OffDA.*OffA)/sum(OffA); %this should be the same as OffPoverATotal
    OnDoverATotalWeightedByOffA = sum(OnDA.*OffA)/sum(OffA);
    
    GradBi.OffDendEccInnerLimit = OffDendEccInnerLimit;
    GradBi.OffDendEccOuterLimit = OffDendEccOuterLimit;
    GradBi.EccBin = EccBin;
    GradBi.Eccentricity = Ecc;
    GradBi.OffPvsEcc = OffP;
    GradBi.OffDvsEcc = OffD;
    GradBi.OffAvsEcc = OffA;
    GradBi.OffPoverDvsEcc = OffPD;
    GradBi.OffPoverAvsEcc = OffDA;
    GradBi.OffDoverAvsEcc = OffDA;
    GradBi.OnPvsEcc = OnP;
    GradBi.OnDvsEcc = OnD;
    GradBi.OnAvsEcc = OnA;
    GradBi.OnPoverDvsEcc = OnPD;
    GradBi.OnPoverAvsEcc = OnDA;
    GradBi.OnDoverAvsEcc = OnDA;
    GradBi.OffPTotal = OffPTotal;
    GradBi.OffDTotal = OffDTotal;
    GradBi.OffATotal = OffATotal;
    GradBi.OffPoverDTotal = OffPoverDTotal;
    GradBi.OffPoverATotal = OffPoverATotal;
    GradBi.OffPoverATotal = OffPoverATotal;
    GradBi.OnPTotal = OnPTotal;
    GradBi.OnDTotal = OnDTotal;
    GradBi.OnATotal = OnATotal;
    GradBi.OnPoverDTotal = OnPoverDTotal;
    GradBi.OnPoverATotal = OnPoverATotal;
    GradBi.OnPoverATotal = OnPoverATotal;
    GradBi.OnPoverDTotalWeightedByOffD = OnPoverDTotalWeightedByOffD;
    GradBi.OnPoverATotalWeightedByOffA = OnPoverATotalWeightedByOffA;
    GradBi.OnPoverATotalWeightedByOffA = OnPoverATotalWeightedByOffA;
    
    GradBi
    
    save([TPN 'GradBi.mat'], 'GradBi');
end
    
    
    
    
    
    
        




