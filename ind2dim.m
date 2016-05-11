%% Turns linear index yxz into y,x,z index for matrix of size ys,xs,zs

function[y,x,z]=ind2dim(yxz,ys,xs,zs)


y=(mod((mod(yxz-1,ys*xs)+1)-1,ys)+1)';
x=(fix(((mod(yxz-1,ys*xs)+1)-1)./ys)+1)';
z=(fix((yxz-1)./(ys*xs))+1)';