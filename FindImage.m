function[DPN]=FindImage(TPN)

Check=0;
    if isdir(TPN) %first condition Dat(i) must be folder
        
        %%look for image folder
         %create target path name
        d=dir(TPN); %map target folder
        d=d(3:size(d,1)); %clip .  ..
        Check=zeros(size(d,1),1);  %how many did you find
        for c=1:size(d,1)  % Check all files in folder
            if  isdir([TPN  d(c).name]) && ...
                ~strcmp(d(c).name,'images') && ...
                ~strcmp(d(c).name,'pics') && ...
                ~strcmp(d(c).name,'data') && ...
                ~strcmp(d(c).name,'dataFix') && ...
                ~strcmp(d(c).name,'mask') && ...
                ~strcmp(d(c).name,'temp') && ...
                ~strcmp(d(c).name,'other'), %% if not another folder type
                
                %Check files in directory
                Files=dir([TPN  d(c).name]);
                Files=Files(3:size(Files,1));
                if size(Files,1)>2 %at least three files
                for f = 2:size(Files,1)-1 %run all but last file 
                    Tag=Files(f).name;
                    if strcmp('.tif',Tag(1,max(1,size(Tag,2)-3):size(Tag,2))) %check if tif
                        OK=1; %check if all the same size
                        if ~Files(f-1).bytes==Files(f).bytes, 
                            OK=0; %not the same size
                        end
                        Check(c)=OK;  
                    end %check all files
                  end %run all middle folders
                end% end if folder big enought
            end% End if directory
        end% End check all files in folder
    end %if Data file is folder
    if sum(Check)==1
        DPN=[TPN d(find(Check,1)).name '\'];
    else
        DPN='?';
    end