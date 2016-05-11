function[I]=showvec(Va,siz);
%%Show standard vector centered at 0


V=Va(:,1:2); %only use first two dimensions of vector
V=round(V);

    I=zeros(siz,'uint8');
    difi=siz/2;
    V(:,1)=V(:,1)+difi(1);
    V(:,2)=V(:,2)+difi(2);
    

Out=   V(:,1)>siz(1) | V(:,2)>siz(2) | V(:,1)<1 | V(:,2)<1;

iVf=V(~Out,:);

for i = 1:size(iVf,1)
    I(iVf(i,1),iVf(i,2))=I(iVf(i,1),iVf(i,2))+1;
end