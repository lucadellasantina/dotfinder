function[R P] = jmscat(x,y);

fsize=500;

colormap gray(255)

maxX=max(x); maxY=max(y);
minY=0; minX=0;
scx=fsize/(maxX-minX);
scy=fsize/(maxY-minY);

Field=zeros(fsize,fsize);


for i = 1 :size(x,1)
    Field(max(round(y(i)*scy),1),max(round(x(i)*scx),1))=Field(max(round(y(i)*scy),1),max(round(x(i)*scx),1))+1;
end

Field=flipdim(Field,1);

%Field=Field*200/max(Field(:));
%Field(Field>0)=Field(Field>0)+50;
Field(Field>0)=Field(Field>0)*50+50;
image(Field),pause(.1)