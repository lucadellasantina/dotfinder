function[] = HOanaCBGrouped(TPN)

% This program calculate the position of the center of cell body as mean of
% all voxels assigned as a cell body region in Amira (assigned as 2).
% Then, calculate distance from this cell body center to individual dots.

% New version will assume the center of cell body as a starting point of
% marching during Imaris filament skeletonization whose location is stored
% in Skel.FilStats.SomaPtXYZ. HO 1/10/2010
%
% 7/6/2010 HO modified the code to be applicable to the grouped dots. It
% simply loads Grouped instead of Dots, convert the name from Grouped to
% Dots to take advantage of already written program, then at the end,
% convert it back to Grouped and save Grouped.

%Find cell body

%colormap gray(255)

load([TPN 'Settings.mat'])
ImageInfo = evalin('base', 'Settings');
ImageInfo = ImageInfo.ImInfo;
xyum=ImageInfo.xyum; %changed to reflect structure format of ImInfo HO 1/5/2010
zum=ImageInfo.zum; %changed to reflect structure format of ImInfo HO 1/5/2010


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'

%Instead of loading Dots, load Grouped and change the Grouped to Dots to
%take advantage of already written program. 7/6/2010 HO
%load([TPN 'Dots.mat'])
load([TPN 'Grouped.mat'])
Dots = Grouped;

%%{

% TPNm = [TPN 'mask\'];
% dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));
% 
% if size(dTPNm,1)==1
%     Ir=tiffread2([TPNm dTPNm(1).name]);
%     IM=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
%     for i = 1: size(Ir,2)
%         IM(:,:,i)=Ir(i).data;
%     end
%     clear Ir
% else
% 
%     clear I Ic
%     IM(:,:)=imread([TPNm dTPNm(1).name]); %read
%     IM(1,1,size(dTPNm,1))=0;
%     c=0;
%     for i=1:size(dTPNm,1)
%         nam=dTPNm(i).name;
%         naml=length(nam);
%         if nam(naml-3:naml)=='.tif'
%             c=c+1;
%             IM(:,:,c)=imread([TPNm nam]);
%         end
%         PercentRead=i/size(dTPNm,1)*100
%     end
% 
% 
% end
% 
% [y x z] = ind2sub(size(IM),find(IM==2));
% CBpos=[mean(y) mean(x) mean(z)];
% Dots.Im.CBpos=CBpos;
%}

load([TPN 'Skel.mat'])

CBpos = [ceil(Skel.FilStats.SomaPtXYZ(2)/xyum) ceil(Skel.FilStats.SomaPtXYZ(1)/xyum) ceil(Skel.FilStats.SomaPtXYZ(3)/zum)];
Dots.Im.CBpos=CBpos;

CBpos=Dots.Im.CBpos; %open and scale
CBpos(1:2)=CBpos(1:2)*xyum; CBpos(3)=CBpos(3)*zum;
Dpos=Dots.Pos; Dpos(:,1:2)=Dpos(:,1:2)*xyum; Dpos(:,3)=Dpos(:,3)*zum;

Dist2CB=dist(Dpos,CBpos);
Dots.Dist2CB=Dist2CB;

%Instead of saving Dots, change the Dots back to Grouped
%and save Grouped. 7/6/2010 HO
%save([TPN 'Dots.mat'],'Dots')
Grouped = Dots;
save([TPN 'Grouped.mat'],'Grouped')




