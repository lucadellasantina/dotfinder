function[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);%clear all


c=0;

NumCells=size(UseCells,1);
Age=zeros(NumCells,1);
ONs=Age; OFFs=Age; Others=Age; BIs=Age;
Type=cell(NumCells,1);

for k = 1:NumCells
     TPN = char(UseCells(k,:)); 
     
    if exist([TPN 'CA.mat'])
        c=c+1;
        load([TPN 'CA.mat'])
        load([TPN 'Cell.mat'])        
        
        Name(c)={Cell.Name};
        Age(c)=str2num(Cell.Age);
        Type{c}=Cell.Type;
        
                  %% read Cell type and mark appropriate matrix
    switch Cell.Type
        case 'ON'
            ONs(c)=1;
        case 'OFF'
            OFFs(c)=1;
        case 'BI'
            BIs(c)=1;
        case 'Other'
            Others(c)=1;
        otherwise
            Unknowns(c)=1;
    end %End Cell type switch

    end
end
ONs=ONs>0; OFFs=OFFs>0; BIs=BIs>0; Others=Others>0;
Age(Age>30)=35;

