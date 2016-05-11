function[IM] = CountSurround(IM)
% replaces binary matrix with matrix with number from 1:7 representing 
% how many 6 connected voxels are present

%IM=IM2;
IM=IM>0;
IM=uint8(IM);

[ys xs zs]=size(IM);

%% collect xy
xySur = [0 1 0; 1 1 1; 0 1 0];
for i = 1 : zs
   fPlane = IM(:,:,i);  
   fPlane = imfilter(fPlane,xySur,'same');
   fPlane = fPlane .* IM(:,:,i);
   IM(:,:,i) = fPlane;
end

%% collect z
for i = 1: zs + 1

    zPlaneNew=uint8(sum(IM(:,:,max(1,i-1):min(zs,i+1))>0,3));
    if i>1
        IM(:,:,i-1) = IM(:,:,i-1) + (zPlaneOld-1).* uint8(IM(:,:,i-1)>0);
    end
    zPlaneOld = zPlaneNew;
end


            
            
            