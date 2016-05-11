function[I]=showvec(Va,siz);

V=Va(:,1:2); %only use first two dimensions of vector

mins=min(V,[],1);
maxs=max(V,[],1);
difs=maxs-mins;

iV(:,1)=V(:,1)-mins(1);
iV(:,2)=V(:,2)-mins(2);
mids=mean(iV,1);

if exist('siz')
    I=zeros(siz,'uint8');
    difi=siz/2-mids;
    iV(:,1)=iV(:,1)+difi(1);
    iV(:,2)=iV(:,2)+difi(2);
    
else
    I =zeros(fix(difs)+1,'uint8')
end

iVf=fix(iV)+1;
iVf(iVf(:,1)>siz(1),1)=siz(1);
iVf(iVf(:,2)>siz(2),2)=siz(2);
iVf(iVf(:,1)<1,1)=1;
iVf(iVf(:,2)<1,2)=1;

for i = 1:size(iVf,1)
    I(iVf(i,1),iVf(i,2))=I(iVf(i,1),iVf(i,2))+1;
end