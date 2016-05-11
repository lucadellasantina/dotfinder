function[Grid]=hexes(siz,Diam)
%% Generate list of yx positions for centers of hexagonal grid covering
%% field of size siz, with an inter hex distance of Diam. 


c=0;
for y = 0-Diam:Diam:siz(1)+Diam;
    Shift=Diam/2*mod(c,2);
    for x = 0-Diam:Diam:siz(2)+Diam;
        c=c+1; 
        Grid(c,1)=y; 
        Grid(c,2)=x+Shift;
    end
end
