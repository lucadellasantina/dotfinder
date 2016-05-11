function[d]=dist2(A,B)

%finds distance between two vectors in form A (n,3,n) and B (1,3) 
A=double(A); B=double(B);
d=sqrt((A(:,1,:)-B(1)).^2 + (A(:,2,:)-B(2)).^2 );