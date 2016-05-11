%% GET DIRECTORY OF CELL
TPN = GetMyDir;
%yxum=0.103;
%zum=0.3;

%% DEPTHDEV iss based on z, here use real stratification 7/30/2010 HO
%do this before use so that I can add my stratification in Use.
%but skip this if you don't have Strat.mat. Everything will work fine
%without Strat.mat if the cell is monostratified for now 11/10/2011 HO
%[Strat] = HOStratificationAnalysis(TPN);

%% GENERATE USE
%[Use] = HOMakeUse2007(TPN, yxum, zum); %#ok<NASGU> %Exactly the same
close all;
LDSanaMakeUseOnec(TPN)

LDSDotsDD(TPN); %input argument of yxum and zum were removed, instead load Settings to get those within the program.
%Although I don't like the way this program defines stratification-related
%parameters, just keep the format as it is.
%Can be modified in the future for bistratified GCs.

%% GENERATE CA FOR AREA
close all;
[CA] = LDSCAsampleUse(TPN);

%% DRAW OUTPUT
LDSCAsampleCollect(TPN);

%% GENERATE DEPTHDEV
%Not necessary to do this in each arbor for bistratified RGC.
%So the program was not modified for bistratified RGC. 10/15/2011 HO.
close all;
[DepthDev] = LDSStratNoCBmedian(TPN);

%% GENERATE GRAD
[Grad] = LDSGradient(TPN);
