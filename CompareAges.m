function[] = PlotAges(Age,Data,Age2,Data2)


%%Data
scatter(Age,Data,'b')
hold on
%%Find standard Error as StandardDeviation / sqrt(n)
E(1)=std(Data(Age==5))/sqrt(sum(Age==5));
E(2)=std(Data(Age==7))/sqrt(sum(Age==7));
E(3)=std(Data(Age==12))/sqrt(sum(Age==12));
E(4)=std(Data(Age>30))/sqrt(sum(Age>30));


errorbar([5 7 12 35],[mean(Data(Age==5)) mean(Data(Age==7)) ...
    mean(Data(Age==12)) mean(Data(Age>30))],E,'c')




scatter(Age2,Data2,'r')
%%Find standard Error as StandardDeviation / sqrt(n)
E(1)=std(Data2(Age2==5))/sqrt(sum(Age2==5));
E(2)=std(Data2(Age2==7))/sqrt(sum(Age2==7));
E(3)=std(Data2(Age2==12))/sqrt(sum(Age2==12));
E(4)=std(Data2(Age2>30))/sqrt(sum(Age2>30));


errorbar([5 7 12 35],[mean(Data2(Age2==5)) mean(Data2(Age2==7)) ...
    mean(Data2(Age2==12)) mean(Data2(Age2>30))],E,'m')



ys=YLim;
ylim([0 ys(2)]);

hold off
