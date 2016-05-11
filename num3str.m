function[s]=num3str(num,d)
%converts 'num' to a 'd' didgit string


n(d)=mod(num,10);
if d>1, n(d-1)=(mod(num,100)-n(1))/10; end
if d>2, n(d-2)=(mod(num,1000)-n(2)*10-n(1))/100; end
if d>3, n(d-3)=(mod(num,10000)-n(3)*100-n(2)*10-n(1))/1000;end

for i=1:d
    s(i)=char(n(i)+48);
end
