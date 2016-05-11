function[n e]=findv(List,Targ)

%Finds the index in an [n X 3 X n] dimension 
%list of vectors of a particular target vector (1 X 3)
List=double(List);
Targ=double(Targ);
[n,e]=find(List(:,1,:)==Targ(1) & List(:,2,:)==Targ(2) & List(:,3,:)==Targ(3));
  