function LDSCAsampleCollect(TPN)

%7/8/2010 HO added comments.
%This program doesn't do any analysis, but just generates different CA4,
%which is basically the same as the previous CA4, but adding +2 to the all
%the pixels within the territory so that outside the territory remains to
%be black but the pixes with zero value (no dots arround) inside the
%territory gets at least 2 so that it becomes blue and differentiated from
%outside the territoty. It also creates figures with slightly different
%range of P/A,D/A,P/D spread to 8-bit scale, which was now noted in the
%file name as well as in the title of the figure. 7/30/2010 HO
%3/11/2011 AB and HO added colorbar for the final saving image.

close all;
cfigure(40,10);
colormap(jet(256))
cmap=colormap;
cmap(1,:)=0;
colormap(cmap)

% load('UseCells.mat')
% KPN = '\\128.208.64.36\wonglab\Josh\Analyzed\Output2\'

if exist([TPN 'CA.mat'])
    load([TPN 'CA.mat'])
    
    for a = 1:length(CA.Arbor); %HO added the loop 10/15/2011 HO

        DotMap=CA.Arbor(a).DotMap;
        DendMap=CA.Arbor(a).DendMap;
        DotDist=CA.Arbor(a).DotDist;
        DendDist=CA.Arbor(a).DendDist;
        Territory =CA.Arbor(a).Territory;

        [ys xs]=size(DotMap);

        DotDend=DendMap;
        DotDend(DotDend>0)=245;
        DotDend(DotMap>0)=100;
        image(DotDend)

        CA4=zeros(ys,961);
        CA4(1:ys,481:481+xs-1)=DotDist*1500+Territory*2;
        CA4(1:ys,241:241+xs-1)=DendDist*300+Territory*2;
        CA4(1:ys,1:xs)=DotDend;
        DDi=DotDist./DendDist*256+Territory*2; %(256*desidered max)
        DDi(DendDist==0)=Territory(DendDist==0)*2;
        CA4(1:ys,721:721+xs-1)=DDi;
        image(CA4)
        title('From left: Dot+DendMap, D/A (0-0.85um/um2), P/A (0-0.17puncta/um2), P/D (0-1puncta/um)');
        pause(5)


        Perimeter=bwperim(Territory);

        if a==1;
            Name=[TPN 'images\DA_0-pt85UmPerUm2&PA_0-pt17PunctaPerUm2&PD_0-pt85PunctaPerUm.tif'];
        else
            Name=[TPN 'images\DA_0-pt85UmPerUm2&PA_0-pt17PunctaPerUm2&PD_0-pt85PunctaPerUm_Arbor' num2str(a-1) '.tif'];
        end
        %Name2=[TPN 'images\Territories_' num2str(a) 'tif'];

        %no point of saving these HO
    %     if isempty(find(Name=='?'));
    %         imwrite(CA4,cmap,Name,'Compression','none')
    %         imwrite(Perimeter,Name2,'Compression','none')
    %     end

        %Adam added color bar 3/11/2011
        colorbar('Ytick',[]);
        axis image; %adding color bar changes image x y dimention, so resize the image.
        axis off %remove the axis

        saveas(gcf,Name) %save figure with title
    end
end



