function[Num] = str2numCell(Str)

for x = 1:size(Str,1)
    for y = 1:size(Str,2)
        Num{x,y}=str2num(Str{x,y})
    end
end