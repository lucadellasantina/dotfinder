%% anaRead 
%   script for reading in tiff files of various formats
%% Comments and error log
% 051409 aab added stand alone TPN line 3
%   added write out 'done writing' line 93
% 1/5/2010 HO changed the tiff format to read because our new method uses
%   Fiji which saves a single multi-image (3-D) tiff file for each channel.
% 1/10/2010 HO changed 'title' to 'name' for at line 70 and 72 because 
%   there is Matlab function called 'title', which somehow caused
%   confliction.

%%
function[] = LDSanaRead(TPN)

if ~exist('TPN')
    TPN = GetMyDir;
end

colormap gray(255)

TPNi=[TPN 'I' filesep]
dTPNi=dir([TPNi '*.tif']); %by adding '*.tif' you can extract only tif files (9/2/09 HO)

%% Get dims from the first z-section
Is=imread([TPNi dTPNi(1).name]);
Iclass=class(Is);
% [ys xs cs]=size(Is); %old method where Amira outputs multiple 2-D multi-channel tiff files
% zs=size(dTPNi, 1); 
[ys xs]=size(Is); %new method where Fiji outputs a single multi-image (3-D) tiff file for each channel
zs = length(imfinfo([TPNi dTPNi(1).name]));
cs = size(dTPNi, 1);
%% read
Iraw=zeros(ys, xs, zs, cs,Iclass);
'reading'
% for i = 1:zs %old method where Amira outputs multiple 2-D multi-channel tiff files
%    Iraw(:,:,i,:)=imread([TPNi dTPNi(i).name]);
% end
for i = 1:cs %new method where Fiji outputs a multi-image (3-D) tiff file for each channel
    for j = 1:zs
        Iraw(:,:,j,i)=imread([TPNi dTPNi(i).name], j);
    end
end
'done reading'


Imax=squeeze(max(Iraw,[],3));
if xs > 2*ys
    cfigure(40, 40*(ys/xs));
else
    cfigure(20*(xs/ys), 20);
end
for i = 1:size(Imax,3)    
    subplot(2,2,i)
    image(Imax(:,:,i)*(500/double(max(max(Imax(:,:,i))))))
    titleword = ['channel ' num2str(i)];
    title(titleword);
end
colormap gray(256);
if size(Imax,3)==3; %fixed the error, now working fine with 2-channel data 2/11/2010 HO
    subplot(2,2,4),image(Imax)
elseif size(Imax,3)==2;
    Imax2=Imax;
    Imax2(:,:,3)=0;
    subplot(2,2,4),image(Imax2)
end


%% Get image info, changed to use structure for ImInfo instead of Cell (9/2/09 HO)

clear v
v.DenCh = 3;
v.PostCh = 2;
v.ColoCh = 1;
v.xyum = .103;
v.zum=.3;
v.MedFilt = 0;
v.MedFiltKern = 0;

'getting image info', pause(.1)
prompt = {'dendrite channel: ', 'post synaptic channel', 'colocalizing channel:  ',...
    'xy resolution : ', 'z resolution :','median filter post? (1=yes, 0=no)','median kernel dimension: '};
name = 'Define image channels. (0 if not applicable)'; %'title' to 'name' because there is Matlab function called 'title', which somehow conflicted at this line HO 1/10/2010

v = getVars(v, name, prompt); %'title' to 'name' because there is Matlab function called 'title', which somehow conflicted at this line HO 1/10/2010

v.xNumVox = xs; %added by HO 12/16/09
v.yNumVox = ys; %added by HO 12/16/09
v.zNumVox = zs; %added by HO 12/16/09

ImInfo = v;

DenCh=ImInfo.DenCh;
PostCh=ImInfo.PostCh;
ColoCh=ImInfo.ColoCh;
xyum=ImInfo.xyum;
zum=ImInfo.zum;
MedFilt=ImInfo.MedFilt;
MedFiltKern=ImInfo.MedFiltKern;

if exist([TPN 'Settings.mat'])
    load([TPN 'Settings.mat']);
end
Settings.ImInfo = ImInfo;
save([TPN 'Settings.mat'], 'Settings') %fixed to save only Settings (9/2/09 HO)

%% Write Channels


if DenCh,  Dend=Iraw(:,:,:,DenCh); save([TPN 'Dend.mat'],'Dend');end
if PostCh, Post=Iraw(:,:,:,PostCh);
    if MedFilt
        for i=1:size(Post,3)
            Post(:,:,i)=medfilt2(Post(:,:,i),[MedFiltKern,MedFiltKern]);
        end
        'Medean filter applied'
    end
    save([TPN 'Post.mat'],'Post');
end
if ColoCh, Colo=Iraw(:,:,:,ColoCh); save([TPN 'Colo.mat'],'Colo'); end
'done writing'

close all;


%% 













