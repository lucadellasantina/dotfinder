function[ME] =MyPlot(Data,Lab)

[ys xs]=size(Data);

if nargin == 1
    Lab= 1:xs
end

%%Data


Holding = ishold;
plot(0)
hold on
for a = 1:xs
    scatter(ones(ys,1)*Lab(a),Data(:,a)) 
    %%Find standard Error as StandardDeviation / sqrt(n)
    E(a)=std(Data(:,a))/sqrt(ys);
    M(a)=mean(Data(:,a));
end
errorbar(Lab,M,E,'b')
ME=[M'  E'];

scatter(Lab,M,'r')

hold off


if Holding
    hold on
end
