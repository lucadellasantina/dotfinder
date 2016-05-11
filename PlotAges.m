function[] = PlotAges(Age,Data)



%%Data
scatter(Age,Data)
hold on
%%Find standard Error as StandardDeviation / sqrt(n)
E(1)=std(Data(Age==5))/sqrt(sum(Age==5));
E(2)=std(Data(Age==7))/sqrt(sum(Age==7));
E(3)=std(Data(Age==12))/sqrt(sum(Age==12));
E(4)=std(Data(Age>30))/sqrt(sum(Age>30));


errorbar([5 7 12 35],[mean(Data(Age==5)) mean(Data(Age==7)) ...
    mean(Data(Age==12)) mean(Data(Age>30))],E,'r')


hold off

ys=YLim;
ylim([min(0,ys(1)) ys(2)]);