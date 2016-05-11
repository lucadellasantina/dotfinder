function[D]=skel3(D,num);

for i=1:num
D=~D;
Dst = bwdist(D,'chessboard'); % Distance Transform (Euclidean) ref. p.6
DstPos = Dst>0;
H = repmat(1/26,[3 3 3]); H(2,2,2)=0; % mask 
DFlt = imfilter(Dst,H); % filtered distance, ref. p.8 Note: 3d mean flt not available in matlab
tp=.52 % Thinness parameter control tp ref. p.8 
% An high (a low) tp value results in a thinner (thick) skeleton
D = DFlt<(Dst-tp); % Skeleton , ref. p.9 
image(max(D,[],3)*1000)
pause(1)
image(D(:,:,30)*1000)
pause(.01)
end

