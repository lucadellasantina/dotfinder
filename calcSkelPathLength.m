function Skel = calcSkelPathLength(Skel)
%This program was modified from HOSkelFinerGenerator.mat.
%This program takes in Skel, march from soma to connected skels while
%registering the path length from soma for each skel point. Then put the
%calculated path length from soma under Skel.FilStats and return.

debug = 0;

SkelIDsPool = 1:1:size(Skel.FilStats.aXYZ,1);
SkelIDsConnectivityPool = Skel.FilStats.aEdges+1; %+1 because Imaris ID starts from zero.
SkelSegLengthsPool = Skel.SegStats.Lengths;

SomaPtID = Skel.FilStats.SomaPtID+1; %+1 because Imaris ID starts from zero.
SomaPtXYZ = Skel.FilStats.SomaPtXYZ;
SourceSkelIDs = SomaPtID; %start from soma.
SkelIDsPool(SomaPtID) = []; %deplete the soma point.

SkelPathLength2Soma = zeros(1,size(Skel.FilStats.aXYZ,1));

% March from the soma
txtBar('Marching through current skeleton to calculate lengths ... ');
while ~isempty(SkelIDsPool)
    txtBar(100 - 100 * numel(SkelIDsPool)/size(Skel.FilStats.aXYZ,1) +1);
    %length(SkelIDsPool) %how many more skel to deplete?
    if isempty(SourceSkelIDs) %if the previous skel was the dead end, resume from the first entry within the remaining pool.
        SourceSkelIDs = SkelIDsPool(1); 
        SkelIDsPool(1) = []; %deplete the newly grabbed skel.
    end

    NextSkelIDs = [];
    for j=1:length(SourceSkelIDs)
        SourceSkelID = SourceSkelIDs(j);
        [SourceSkelRows, SourceSkelCols] = find(SkelIDsConnectivityPool == SourceSkelID);
        if ~isempty(SourceSkelCols) %if the source skel is not dead end, find the connected partner skel
            NextSkelCols = SourceSkelCols*(-1)+3; %to reverse 1 and 2 to get the partner skel.
            for i=1:length(NextSkelCols)
                SourceSkelRow = SourceSkelRows(i);
                SourceSkelCol = SourceSkelCols(i);
                NextSkelCol = NextSkelCols(i);
                NextSkelID = SkelIDsConnectivityPool(SourceSkelRow, NextSkelCol);
                Lengths2NextSkel = SkelSegLengthsPool(SourceSkelRow);
                
                SkelPathLength2Soma(NextSkelID) = SkelPathLength2Soma(SourceSkelID) + Lengths2NextSkel;
                
                NextSkelIDs = [NextSkelIDs, NextSkelID];
            end
            SkelIDsConnectivityPool(SourceSkelRows,:) = []; %deplete the found connectivity from the pool.
            SkelSegLengthsPool(SourceSkelRows) = []; %also deplete seg lengths
        end
    end

    if ~isempty(NextSkelIDs)
        for id=1:length(NextSkelIDs)
            SkelIDsPool(SkelIDsPool==NextSkelIDs(id)) = []; %deplete the next grabbed skel.
        end
    end
    SourceSkelIDs = NextSkelIDs; %switch next to source for the next loop.
end
txtBar('DONE');

if debug
    disp(['Farthest skel path distance is: ' max(SkelPathLength2Soma)]);
    disp(['Soma point ID is (Imaris soma pt ID + 1): ' num2str(SomaPtID)]);
    disp(['Skel IDs of zero path distance are (only soma point should have zero path distance): ' num2str(find(SkelPathLength2Soma==0))]);
end;

% Calculated path length of edges by taking the mean of skels
EdgePathLength2Soma = zeros(1,size(Skel.FilStats.aEdges,1));
for i=1:length(EdgePathLength2Soma)
    EdgePathLength2Soma(i) = mean([SkelPathLength2Soma(Skel.FilStats.aEdges(i,1)+1), SkelPathLength2Soma(Skel.FilStats.aEdges(i,2)+1)]);
end

% Store the results under Skel.FilStats
Skel.FilStats.SkelPathLength2Soma = SkelPathLength2Soma;
Skel.FilStats.EdgePathLength2Soma = EdgePathLength2Soma;
end