function[] = Stats(In);
%Colums equal data sets, Rows equal Ns
In = In(:);
[ys xs zs]=size(In);
Ns=ys
Mean=mean(In,1)
Median=median(In,1)
SE=std(In)./sqrt(Ns)