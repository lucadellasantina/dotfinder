%function[] = pic2dat(PicFold)
%reads all the images in a folder, compiles them as a 3D matrix and saves
%the matrix as a mat file in the same directory that the folder is in

%% Read directory
d=dir(PicFold); %read selected directory
d=d(3:size(d,1)); %eliminate . and ..

GetSize=imread([PicFold '\' d(1).name]);
I=zeros(size(GetSize,1),size(GetSize,2),size(d,1),'uint8');
for i = 1 : size(d,1)
   I(:,:,i)=imread([PicFold '\' d(i).name]); 
end

%% Rename I to folder name
Back=find(PicFold=='\');
Back=Back(size(Back,2))+1;
Name=PicFold(Back:size(PicFold,2));
evalin('caller',[Name '=' 'I' '; clear ' 'I']) %Switches file names

%% Save 
save([PicFold '.mat'],Name)
%clear(Name)