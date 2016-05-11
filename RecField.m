g=fspecial('gaussian',50,10);
g=single(g*(255/max(g(:))));
image(g)
Rec=convn(Ic,g);