function[show] = BootCorr(x,y)
%%Generates a nonparametric P value for a correlation coefficient.

x=x(:);y=y(:);
xs=size(x,1);
ys=size(y,1);

[R P]=corrcoef(x,y);
R=R(1,2);

reps =10000;
randR=zeros(reps,1);
for r = 1: reps
    randx=fix(rand(xs,1)*xs)+1;
    randy=fix(rand(ys,1)*ys)+1;
    sampx=x(randx);
    sampy=y(randy);
    rR=corrcoef(sampx,sampy);
    randR(r)=rR(1,2);
    
end
    
p=sum(abs(randR)>=abs(R))/reps;
    
show.R=R;
show.studentTP=P(1,2);
show.bootP=p;
