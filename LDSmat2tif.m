function LDSmat2tif(VarName);

% This part is made for the case where you start without input argument,
% but didn't work when I called this function from outside because the
% variable you want to save is outside of this function

% whos %list all the working variables
% 
% ImName = input('Type the name of your working variable that you want to save as tif. \n', 's');
% 
% [FileName,PathName,FilterIndex] = uiputfile('*.tif','Name&Save the tif file.');
% 
% EvalStr = ['NumZ = size(' ImName ',3);'];
% eval(EvalStr); %check number of z planes
% 
% EvalStr = ['imwrite(' ImName '(:,:,1), [PathName FileName], ''tif'', ''compression'', ''none'');'];
% eval(EvalStr); %write first z plane
% 
% if NumZ > 1; %write the rest of the z planes
%     EvalStr = ['imwrite(' ImName '(:,:,i), [PathName FileName], ''tif'', ''compression'', ''none'', ''WriteMode'', ''append'');'];
%     for i=2:NumZ
%         eval(EvalStr);
%     end
% end

NumZ = size(VarName,3);
[FileName,PathName,FilterIndex] = uiputfile('*.tif','Name&Save the tif file.');
imwrite(VarName(:,:,1), [PathName FileName], 'tif', 'compression', 'none'); %write first z plane
if NumZ > 1; %write the rest of the z planes
    for i=2:NumZ
        imwrite(VarName(:,:,i), [PathName FileName], 'tif', 'compression', 'none', 'WriteMode', 'append');
    end
end

   

