
%function[DisBins,DistCorr]=ShiftCorr(Ig, Ib)
% shift center of green image relative to blue and calculate correlation
% Corrolation in 3D using correlation from corr2
colormap gray(255)

%{
%% Blur image
dic=ones(10,10,10);
Ig=convn(Ig,dic,'same');
%}

Ig=Igm; Ib=Ibm;

[yend,xend,planes]=size(Ig);

%{
%tests maricies
Ig=ones(yend,xend,planes)*100;
for y=1:yend, for x=1:xend, for z=1:planes
            Ib(y,x,z)=rand*10;
end,end,end
%}





s=5; %shift distance
c=0; %start counter
for y=1:2*s+1, for x=1:2*s+1
yshift=y-s-1; xshift=x-s-1;

%Igs=Ig*0;
Igs=Ig(s+1:yend-s,s+1:xend-s,:);
Ibs=Ib(s+1+yshift:yend-s+yshift,s+1+xshift:xend-s+xshift,:);

subplot(2,1,1)
image(max(Igs,[],3)*(255/max(Igs(:)))),pause(.01)
subplot(2,1,2)
image(max(Ibs,[],3)*(255/max(Ibs(:)))),pause(.01)

A=Igs; %define first matrix
B=Ibs; %define second matrix
Bm=mean(B(:)); Am=mean(A(:));  %find mean of matrix
Ad=A-Am;  Bd=B-Bm;  %Make matrix that is difference from mean
ABd=Ad.*Bd; AAd=Ad.^2; BBd=Bd.^2;  %matrix multiplication
r=sum(ABd(:))/sqrt((sum(AAd(:)))*(sum(BBd(:)))); %find r
c=c+1; %advance counter
DistCorr(c,1)=sqrt(yshift^2+xshift^2); %find distance to origin
DistCorr(c,2)=r; %record correlation value
DistProd(c,1)=sqrt(yshift^2+xshift^2);
DistProd(c,2)=sum(sum(sum(A.*B)));
clear A B Bm Am Ad Bd ABd AAd BBd

end,PercentShifted=y/(s*2+1),end %run x and y block for shift, plus y counter
%}

for i=1:fix(sqrt(s^2+s^2))+2
    DisBin=find(DistCorr(:,1)<=i-1 & DistCorr(:,1)>i-2); %bin correlation by distance
    if isempty(DisBin), DisBins(i)=0; else
        DisBins(i)=(sum(DistCorr(DisBin,2)))/size(DisBin,1);
    end
end


plot(DisBins)

%{ 
%% strait product
for i=1:fix(sqrt(s^2+s^2))+2
    DisBinP=find(DistProd(:,1)<=i-1 & DistProd(:,1)>i-2); %bin correlation by distance
    if isempty(DisBinP), DisBinsP(i)=0; else
        DisBinsP(i)=mean(DistProd(DisBinP,2)); end
end

subplot(3,2,6)
plot(DisBinsP)
%}


pause(.01)








