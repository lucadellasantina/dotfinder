%% Dots count
numel(find(SG.passF))

%% All Dots visualization

A = zeros(Settings.ImInfo.yNumVox, Settings.ImInfo.xNumVox, Settings.ImInfo.zNumVox, 'uint8');

for i=1:Dots.Num;
    A(Dots.Vox(i).Ind)=1;
end

LDSmat2tif(A);
clear A
%% Passing Dots visualization

A = zeros(Settings.ImInfo.yNumVox, Settings.ImInfo.xNumVox, Settings.ImInfo.zNumVox, 'uint8');

for i=1:Dots.Num;
    if SG.pass1(i)
        A(Dots.Vox(i).Ind)=1;
    end
end

LDSmat2tif(A);
clear A
%% Rejected Dots visualization

A = zeros(Settings.ImInfo.yNumVox, Settings.ImInfo.xNumVox, Settings.ImInfo.zNumVox, 'uint8');

for i=1:Dots.Num;
    if ~SG.pass1(i)
        A(Dots.Vox(i).Ind)=1;
    end
end

LDSmat2tif(A);
clear A