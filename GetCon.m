function[C] = GetCon(T,con)
%%Counts connectivity of matrix 'T' by connectivity 'con'

if ~exist('con'),con=26;end %set default connectivity to 26

%T = T>0;

[ys xs zs] = size(T);


C=zeros(ys+2,xs+2,zs+2,'uint8');
Shift=C;


c=0;
for y = -1 : 1,    for x = -1 : 1,         for z = -1 : 1
            c=c+1;
            Lshift(c,:)=[y x z];
end,end,end
sumL=sum(Lshift.^2,2);

if con==6
    L=Lshift(sumL==1,:);
elseif con == 18
    L=Lshift(sumL==1 | sumL == 2,:)
elseif con == 27
    L = Lshift(sumL>0,:)
end


for i = 1: size(L,1)   
    Shift=Shift*0;
    Shift(2+L(i,1):ys+1+L(i,1),2+L(i,2):xs+1+L(i,2),2+L(i,3):zs+1+L(i,3))=T;
    C=C+Shift;
end

C=C(2:ys+1,2:xs+1,2:zs+1);