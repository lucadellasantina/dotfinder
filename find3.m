function[y,x,z]=find3(matrix)
% finds a number (second argument) in a 3D matrix (first argument)


[ys,xs,zs]=size(matrix);
[y,xz]=find(matrix);
x=mod(xz-1,xs)+1;
z=fix((xz-1)./xs)+1;