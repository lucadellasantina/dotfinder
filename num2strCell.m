function[Str] = num2strCell(Num)


for y = 1: size(Num,1)
for x = 1: size(Num,2)
    Str{y,x}=num2str(Num{y,x})
end
end