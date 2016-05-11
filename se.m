function[SE]=se(Dat)
%% Gives standard error as standard deviation divided by 
SE=std(Dat(:))/sqrt(numel(Dat));