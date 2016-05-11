function aviStack = Avi2Img()
[filename pathname]= uigetfile('*.avi');
 movieFullFileName = [pathname filename];
mov = aviread(movieFullFileName); 
frame=mov(1).cdata;
grayImage = rgb2gray(frame);

aviStack = uint8(zeros(size(mov,2), size(grayImage,1), size(grayImage,2)));
for i=1:size(mov,2)
    frame=mov(i).cdata;
    grayImage = rgb2gray(frame);
    %imshow(grayImage);
    %drawnow;
    aviStack(i,:,:) = grayImage(:,:);
end
% To show individual frames
%imshow(squeeze(I(30,:,:)))

clear filename pathname mov movieFullFileName i frame grayImage;
end