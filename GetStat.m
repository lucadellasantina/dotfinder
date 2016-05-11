function = Stats(In);

[ys xs zs]=size(In);
Ns=ys*xs*zs
Mean=mean(In)
Median=median(In)
SE=std(In)/sqrt(Ns)