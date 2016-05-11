function[Out]=CoRot(In,Co)
%%Rotates 3D vector (n,3) matrix by 3 X 3 princom Coefficient matrix 
%%In the form CoRot(Matrix, Coefficients)

Co=(Co);
Out=In*0;
for i = 1:size(In,1)
    I=In(i,:);
    O = Co.*[I; I ;I];
    Ou=sum(O,2)';
    Out(i,:)=Ou;
    
end

