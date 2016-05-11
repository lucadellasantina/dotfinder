function[DistBins]=shift(Ig, Ib, Ir)
% shift center of green image relative to blue and calculate correlation
% Corrolation in 3D using correlation from corr2

[yend,xend,planes]=size(Ig);

s=5; %shift distance
c=0; %start counter
for y=1:2*s+1, for x=1:2*s+1
yshift=y-s-1; xshift=x-s-1;
Igs=Ig*0;
Igs(s+1+yshift:yend-s+yshift,s+1+xshift:xend-s+xshift,:)=Ig(s+1:yend-s,s+1:xend-s,:);
subplot(1,3,2:3)
image(max(Igs,[],3)*50),pause(.01)


A=Igs; %define first matrix
B=Ib; %define second matrix
Bm=mean(B(:)); Am=mean(A(:));  %find mean of matrix
Ad=A-Am;  Bd=B-Bm;  %Make matrix that is difference from mean
ABd=Ad.*Bd; AAd=Ad.^2; BBd=Bd.^2;  %matrix multiplication
r=sum(ABd(:))/sqrt((sum(AAd(:)))*(sum(BBd(:)))); %find r
clear A B Bm Am Ad Bd ABd AAd BBd
c=c+1; %advance counter
DistCorr(c,1)=sqrt(yshift^2+xshift^2); %find distance to origin
DistCorr(c,2)=r; %record correlation value
DistProd(c,1)=DistCorr(c,1);
DistProd(c,2)=sum(sum(sum(Igs.*Ib)));

end,y,end %run x and y block for shift, plus y counter
%}

for i=1:fix(sqrt(s^2+s^2))+2
    DisBin=find(DistCorr(:,1)<=i-1 & DistCorr(:,1)>i-2); %bin correlation by distance
    if isempty(DisBin), DisBins(i)=0; else
        DisBins(i)=mean(DistCorr(DisBin,2)); end
end
subplot(3,2,4)
plot(DisBins)

for i=1:fix(sqrt(s^2+s^2))+2
    DisBinP=find(DistProd(:,1)<=i-1 & DistProd(:,1)>i-2); %bin correlation by distance
    if isempty(DisBinP), DisBinsP(i)=0; else
        DisBinsP(i)=mean(DistProd(DisBinP,2)); end
end

subplot(3,2,6)
plot(DisBinsP)










