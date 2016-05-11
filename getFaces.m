
nearby = 50;  %Distance at which faces are looked at

Dists = zeros(length(IDs),length(IDs));
for i = 1: length(IDs)
   Dists(i,:) = dist(Dots.Pos(IDs,:),Dots.Pos(IDs(i),:));
  
end
[x y] = find((Dists > 0) & (Dists < nearby));


for i = 1: length(x)
   cDotX = Dots.Vox(x(i)).Pos;
   cDotY = Dots.Vox(y(i)).Pos;
   face = zeros(size(cDotX,1),1);
   for v = 1: size(cDotX,1)
    diff = [cDotY(:,1)-cDotX(v,1)  cDotY(:,1)-cDotX(v,1) cDotY(:,1)-cDotX(v,1)];
    face(v) = sum(sum(abs(diff),2)<2);
   end  
    faces(i) = sum(face);
    
end

x(faces>0)
y(faces>0)
