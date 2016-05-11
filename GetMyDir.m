function[TPN] = GetMyDir()
% allows user to assign a variable to a directory, and maintain that directory in a temp folder.


%% Get Folder 

if ~exist(['.' filesep 'temp']), mkdir(['.' filesep 'temp']);end

if exist(['.' filesep 'temp' filesep 'Last.mat'])
     load(['.' filesep 'temp' filesep 'Last.mat']);
     if exist(Last)
        TPN=uigetdir(Last);
     else
         TPN=uigetdir;
     end
else
    TPN=uigetdir;
end

TPN= [TPN filesep];

Last=TPN;
if Last>0
save(['.' filesep 'temp' filesep 'Last.mat'],'Last')
end
