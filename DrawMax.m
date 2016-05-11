function[Imax]=DrawMax(DPN)

Idir=dir(DPN); Idir=Idir(3:size(Idir,1));

Imax=imread([DPN Idir(1).name]);
for i = 2:size(Idir,1)
   Ip=imread([DPN Idir(i).name]);
   Imax=max(Imax,Ip);
end
