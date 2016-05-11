function I = LDStif2mat(varargin) %first input argument is title but not necessary

%load multi-tif (1 color) file
if length(varargin)>0; %varargin is cell format, so if there is any, the length should be >0.
    if isstr(varargin{1}); %varargin is cell format, so to access the content, bracket has to be used instead of parenthesis.
        [FileName PathName] = uigetfile('*.tif', varargin{1});
        ImInfo = imfinfo([PathName FileName]);
        ImSize = [ImInfo(1).Height ImInfo(1).Width length(ImInfo)]; %[y x z]
        if ImSize(3) == 1;
            I=imread([PathName FileName]);
            I = I(:,:,1); %kill the second and the third z planes if created.
        else
            for j = 1:ImSize(3)
                I(:,:,j)=imread([PathName FileName], j);
            end
        end
    else
        display('ERROR: Input argument has to be string for title or do not enter any input argument!');
        I = [];
        return;
    end
else %no input argument, repeat the above without title input.
    [FileName PathName] = uigetfile('*.tif');
    ImInfo = imfinfo([PathName FileName]);
    ImSize = [ImInfo(1).Height ImInfo(1).Width length(ImInfo)]; %[y x z]
    if ImSize(3) == 1;
        I=imread([PathName FileName]);
        I = I(:,:,1); %kill the second and the third z planes if created.
    else
        for j = 1:ImSize(3)
            I(:,:,j)=imread([PathName FileName], j);
        end
    end
end