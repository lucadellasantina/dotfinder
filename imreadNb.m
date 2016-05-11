function[I1Raw,I2Raw,I3Raw]=imreadN(name,planes);
%% read a stack of tiffs

channels=3;
filetype=2; % 1=multitif, 2=single tifs

%%find number of digits
if planes<10,pdiddy=1; 
elseif planes<100, pdiddy=2; 
elseif planes<1000, pdiddy=3; end


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image...'

if filetype==1 %if multitif
    for i = 1:planes*channels      % read entire stack inot I
    I(:,:,i)=imread(name,i);
    end
    I=single(I(200:630,:,:));  %% Sub sample
    
    %Color seperate
    [ys,xs,zs]=size(I);
    Ig=I(:,:,1:planes);     % Green Channel
    if channels>1
        Ir=I(:,:,planes+1:planes*2);    %Red Channel
    end
    if channels>2
        Ib=I(:,:,planes*2+1:planes*3); %Blue Channel
    end

elseif filetype==2 %if sequence of individual tifs
    ns=size(name,2);
    for i=1:planes
        PercentRead=uint8((i/planes)*100) 
        name2=[name(1:ns-4) num3str(i,pdiddy) name(ns-3:ns)];
        I(:,:,:)=imread(name2);
        I=single(I);
        chan=size(I,3);
        
        %% assign channels
        I1Raw(:,:,i)=I(:,:,1);
        if chan>1        
            I2Raw(:,:,i)=I(:,:,2);
            if chan>2, 
                I3Raw(:,:,i)=I(:,:,3);            
            end %end run channel 3
        end %End run other channells

        %% make max 255
        if max(I1Raw(:))>260 
            I1Raw=I1Raw*(255/max(I1Raw(:))); 
            if chan>1
            I2Raw=I2Raw*(255/max(I2Raw(:))); 
            if chan>2
            I3Raw=I3Raw*(255/max(I3Raw(:)));
            end,end
        end %make 8 bit
        
    end %run planes
end

'done reading'
