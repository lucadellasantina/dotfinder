%% Bin Dots and Dendrites relative to distance from center
%{
this program reads a specified xls file in the form 
(:,1)=length, (:,2)=x1,(:,3)=y1, (:,3)=z1 ect.
moves along line in 0.01um incriments assigning length to appropriate bin
%}

'start'
if ~exist('xyum'),xyum=0.1035;end %pix size in um if not assigned


clear DendDis DotDis
bin=5; %assign bin size in um.
filename='./ManualData/allSeg.xls';
xc=538; %cell center x in pixels
yc=580; %cell center y in pixels
zc=1;
zum=.3;
xc=xc*xyum; %position pixel center in um center
yc=yc*xyum;
zc=zc*zum;
estimate=.1; %incriments of estimating line identity


%%          BIN DENDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%read data 
segs=xlsread(filename,'Sheet1'); %read in segment data
DendDis=zeros(10000,1);

num=max(find(segs(:,3)>0));  %find number of segments
for i=1:num; %run all segments
    
    move=0;
    x1=segs(i,2);y1=segs(i,3);z1=segs(i,4); %extract point1 xyz
    x2=segs(i,5);y2=segs(i,6);z2=segs(i,7); %extract point2 xyz
    xd=x2-x1; yd=y2-y1; zd=z2-z1; %find xyz distance between points for slope
    L=sqrt(xd^2+yd^2+zd^2);
	a=L/estimate;
	%a=estimate^2/(xd^2+yd^2+zd^2); %find xy muliplier to = estimated change
    x3=x1;y3=y1;z3=z1;
    L2=0;
    for move=1:a-1
    	x3=x1+(xd/a)*move;y3=y3+(yd/a)*move;z3=z1+(zd/a)*move;
        d=sqrt((x3-xc)^2+(y3-yc)^2); %find distance of new segment
        DendDis(fix(d/bin)+1)=DendDis(fix(d/bin)+1)+estimate;
        %L2=sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2);
        %trackd(move)=L2;
        
    end
end

last=max(find(DendDis>0));
DendDis=DendDis(1:last);
subplot(3,1,1)
plot(DendDis),pause(.01)

%%      BIN DOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DotDis=zeros(10000,1);

[y,x,z]=find3(I>0,1); %list dot positions
y=y*xyum; x=x*xyum; z=z*.3; %convert to microns
num=size(y,1);

for i=1:num; %run all dots
    d=sqrt((x(i)-xc)^2+(y(i)-yc)^2) %find distance to first point
    DotDis(fix(d/bin)+1)=DotDis(fix(d/bin)+1)+1; %record another dot in bin
end

last=max(find(DotDis>0));
DotDis=DotDis(1:last);
subplot(3,1,2)
plot(DotDis),pause(.01)



%%          COMBINE BOTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sDend=size(DendDis,1); sDot=size(DotDis,1); sb=max(sDot,sDend);
DendDis2=zeros(sb,1);DendDis2(1:sDend)=DendDis;
DotDis2=zeros(sb,1);DotDis2(1:sDot)=DotDis;
DotDendDis=DotDis2./DendDis2
subplot(3,1,1),plot(DendDis2)
subplot(3,1,2),plot(DotDis2)
subplot(3,1,3),plot(DotDendDis),pause(.01)


'stop'



















'stop'


















