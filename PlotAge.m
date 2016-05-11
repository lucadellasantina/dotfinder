function[ME] = PlotAge(Age,Data)

%%Scatter plot Data by Age, then plot mean and standard error for each age

%% Find Ages
if size(Age,2)>1; Age=Age';end %make Age a column
Ages=[];
for i = 1:size(Age,1)
    if isempty(find(Ages==Age(i)))
        Ages=[Ages;Age(i)];
    end
end
Ages=sort(Ages);

%%Data
scatter(Age,Data)
hold on
for a = 1:size(Ages,1)
    %%Find standard Error as StandardDeviation / sqrt(n)
    E(a)=std(Data(Age==Ages(a)))/sqrt(sum(Age==Ages(a)));
    M(a)=mean(Data(Age==Ages(a)));
end
errorbar([Ages],M,E,'r')
ME=[M'  E'];
hold off

ys=YLim;
ylim([0 ys(2)]);