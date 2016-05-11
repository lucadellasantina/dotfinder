function[Pmean,Pmedian]=jmboot(Dat1,Dat2,reps);

if ~exist('reps','var'), reps=10000;end

if size(Dat1,2)>1, Dat1=Dat1'; end
if size(Dat2,2)>1, Dat2=Dat2'; end

meanDat1=mean(Dat1);
Size1=size(Dat1,1);
meanDat2=mean(Dat2);
Size2=size(Dat2,1);

RealDif=abs(meanDat1-meanDat2);
RealMedDif=abs(median(Dat1)-median(Dat2));


%%flawed in that dat size will occasionally be off by one
poolDat=cat(1,Dat1,Dat2);
SizeP=size(poolDat,1);
clear tDif
for i = 1:reps
    List=rand(SizeP,1);
    Ord=sort(List);
    tDat1=poolDat(List<=Ord(Size1));
    tDat2=poolDat(List>Ord(Size1));
    mDat1=mean(tDat1);
    mDat2=mean(tDat2);
    medDat1=median(tDat1);
    medDat2=median(tDat2);
    medDif(i,1)=abs(medDat1-medDat2);
    tDif(i,1)=abs(mDat1-mDat2);
    %showAns=tDif(i,1),pause
end

hist(tDif),pause(.1)
reps;
Pmean=sum(tDif>=RealDif)/size(tDif,1);
Pmedian=sum(medDif>=RealMedDif)/size(medDif,1);
