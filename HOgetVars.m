function[structIn] = HOgetVars(structIn,title,prompt)
%%Input dialouge for changing Structure array
%%Enter the Structure first.  optionally enter title, and then prompts
%%class (num or string) should be preserved

%I made the window size adjustable HO 2/10/2011

fnames = fieldnames(structIn);
if nargin < 3
    prompt = fnames;
end
if nargin < 2
    title = 'Define Variables'
end


nLines = 1;
for i = 1:length(fnames)
    var=getfield(structIn,fnames{i});
    if isstr(var)
        notstr(i)=0;
        gVars{i} = var;
    else
        gVars{i} = num2str(var);
        notstr(i) = 1;
    end
end
options.Resize='on'; %I made the window size adjustable HO 2/10/2011
gVars= inputdlg(prompt,title,nLines,gVars, 'on'); %I made the window size adjustable HO 2/10/2011
pause(.1)
for i = 1:length(fnames)
    if notstr(i)
        structIn=setfield(structIn,fnames{i},str2num(gVars{i}));
    else
        structIn=setfield(structIn,fnames{i},gVars{i});
    end
end