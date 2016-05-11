function[IM] = CountSurround(IM)

%%very slow way to count the nonzero voxels in a 3x3x3 cube at each nonzero
%%voxel.   Has the benefit of not taking additional memory
%%output is a matrix where every non zero value is replaced with a number
%%1:27 representing the number of nonzero voxels in the 3x3x3 cube

%IM=IM2;
IM=IM>0;
IM=uint8(IM);
[ys xs zs]=size(IM);
%% Check core
for y = 2 : ys-1
    for x = 2 : xs-1
        for z= 2 : zs-1
            surround=IM(y-1:y+1,x-1:x+1,z-1:z+1)>0;
            IM(y,x,z)= sum(surround(:))* IM(y,x,z) ;
        end
    end
end

%% Check Edges
for y = [1, ys]
    for x= [1 , xs]
        for z = [1, zs]
            surround=IM(max(1,y-1):min(y+1,ys),max(x-1,1):min(x+1,xs),max(1,z-1):min(z+1,zs))>0;
            IM(y,x,z)= IM(y,x,z)*sum(surround(:));
        end
    end
end
            
           
            
            
            