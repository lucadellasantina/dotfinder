function[]=imwriteN(I,name)
%% write matrix to series of tiffs. function(3dmatrix,'filename')
'writing...'

if isdir('pics')==0, mkdir 'pics'; end %create directory to store steps
folder=['pics/' name];
if isdir(folder)==0, mkdir(folder); else rmdir(folder,'s'), mkdir(folder); end %create directory to store steps

%{
%convert to 255
top=max(I(:));
I=I*(255/top);
%}

I=uint8(I);
%}

%% find planes size and number of digits
 sI=size(I,3);
if sI<10,pdiddy=1; elseif sI<100, pdiddy=2; elseif sI<1000 pdiddy=3; end

for i=1:sI
    %PercentWritten=i/sI*100
    %image(I(:,:,i)),pause(.01)
    name2=[folder '/' name num3str(i,pdiddy) '.tif'];
    imwrite(I(:,:,i),name2,'tif','Compression','none')
end
'done writing'
