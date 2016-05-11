function SkelFiner = HOSkelFinerGenerator(Skel, maxFinerEdgeLength);
%This program takes in Skel.mat, and make finer Skel to fill the voxel gaps
%with the skel to skel distance at least maxFinerEdgeLength. None of the 
%edges will be longer than maxFinerEdgeLength. The result is returned as
%SkelFiner.mat with its format same as Skel.mat.
% HO 7/1/2011.

SkelFinerXYZ = Skel.FilStats.aXYZ; %keep the original skel IDs
SkelFinerEdges = [];
%SkelFinerSeg = [];
%SkelFinerLengths = [];

SkelIDsPool = 1:1:size(Skel.FilStats.aXYZ,1);
SkelIDsConnectivityPool = Skel.FilStats.aEdges+1; %+1 because Imaris ID starts from zero.
SkelSegPool = Skel.SegStats.Seg;
SkelSegLengthsPool = Skel.SegStats.Lengths;

if isfield(Skel.FilStats, 'SomaPtID') %if soma pt is set in Imaris
    SomaPtID = Skel.FilStats.SomaPtID+1; %+1 because Imaris ID starts from zero.
    SomaPtXYZ = Skel.FilStats.SomaPtXYZ;
    SourceSkelIDs = SomaPtID; %start from soma.
    SkelIDsPool(SomaPtID) = []; %deplete the soma point.
else %soma pt not set in Imaris, just grab ID=1 to start
    SourceSkelIDs = 1;
    SkelIDsPool(1) = [];
end


while ~isempty(SkelIDsPool);
    %length(SkelIDsPool) %how many more skel to deplete?
    if isempty(SourceSkelIDs);%if the previous skel was the dead end, resume from the first entry within the remaining pool.
        SourceSkelIDs = SkelIDsPool(1); 
        SkelIDsPool(1) = []; %deplete the newly grabbed skel.
    end

    NextSkelIDs = [];
    Lengths2NextSkels = [];
    for j=1:length(SourceSkelIDs);
        SourceSkelID = SourceSkelIDs(j);
        [SourceSkelRows, SourceSkelCols] = find(SkelIDsConnectivityPool == SourceSkelID);
        if ~isempty(SourceSkelCols) %if the source skel is not dead end, find the connected partner skel
            NextSkelCols = SourceSkelCols*(-1)+3; %to reverse 1 and 2 to get the partner skel.
            for i=1:length(NextSkelCols);
                SourceSkelRow = SourceSkelRows(i);
                SourceSkelCol = SourceSkelCols(i);
                NextSkelCol = NextSkelCols(i);
                NextSkelID = SkelIDsConnectivityPool(SourceSkelRow, NextSkelCol);
                Lengths2NextSkel = SkelSegLengthsPool(SourceSkelRow);
                %this part was modified from anaMa.
                devs=max(1,ceil(Lengths2NextSkel/maxFinerEdgeLength)); %Find number of subdivisions
                NumNewSkels = devs-1;
                NumEdges = NumNewSkels+1;
                if NumNewSkels==0; %if no need to add extra skels
                    SkelFinerEdges(end+1,:) = [SourceSkelID, NextSkelID];
                else %need to add extra skels
                    NewSkelIDs = size(SkelFinerXYZ,1)+1:1:size(SkelFinerXYZ,1)+NumNewSkels;
                    NewSkelIDs = uint32(NewSkelIDs); %need to use the same class
                    clear sy sx sz;
                    for d=1:NumNewSkels;
                        sx(d)=SkelSegPool(SourceSkelRow,1,SourceSkelCol)+((SkelSegPool(SourceSkelRow,1,NextSkelCol)-SkelSegPool(SourceSkelRow,1,SourceSkelCol))/NumEdges)*d;
                        sy(d)=SkelSegPool(SourceSkelRow,2,SourceSkelCol)+((SkelSegPool(SourceSkelRow,2,NextSkelCol)-SkelSegPool(SourceSkelRow,2,SourceSkelCol))/NumEdges)*d;
                        sz(d)=SkelSegPool(SourceSkelRow,3,SourceSkelCol)+((SkelSegPool(SourceSkelRow,3,NextSkelCol)-SkelSegPool(SourceSkelRow,3,SourceSkelCol))/NumEdges)*d;
                    end
                    sy=sy';sx=sx';sz=sz';
                    NewSkelFinerXYZ = [sx, sy, sz];
                    SkelFinerXYZ = [SkelFinerXYZ; NewSkelFinerXYZ];

                    NewSkelIDs = NewSkelIDs';
                    SkelIDSequence = [SourceSkelID; NewSkelIDs; NextSkelID];
                    NewSkelFinerEdges = cat(2, SkelIDSequence(1:end-1), SkelIDSequence(2:end));
                    SkelFinerEdges = [SkelFinerEdges; NewSkelFinerEdges];

                end

                NextSkelIDs = [NextSkelIDs, NextSkelID];
            end
            SkelIDsConnectivityPool(SourceSkelRows,:) = []; %deplete the found connectivity from the pool.
            SkelSegLengthsPool(SourceSkelRows) = []; %also deplete seg lengths
            SkelSegPool(SourceSkelRows,:,:) = []; %also deplete seg
        end
    end

    if ~isempty(NextSkelIDs);
        for id=1:length(NextSkelIDs);
            SkelIDsPool(SkelIDsPool==NextSkelIDs(id)) = []; %deplete the next grabbed skel.
        end
    end
    SourceSkelIDs = NextSkelIDs; %switch next to source for the next loop.
end

%put the finer skel into Skel format.
disp('Total coarse skel length is');
sum(Skel.SegStats.Lengths) %total seg length must be the same between coarse and fine.
disp('Total coarse skel node number is');
size(Skel.FilStats.aXYZ,1)
disp('Total coarse skel edge number is');
size(Skel.FilStats.aEdges,1)
disp('Average coarse skel edge length is');
sum(Skel.SegStats.Lengths)/size(Skel.FilStats.aEdges,1)

SkelFiner.FilStats.aXYZ = SkelFinerXYZ;
SkelFiner.FilStats.aEdges = SkelFinerEdges-1; %-1 to bring it back to Imaris format.

if isfield(Skel.FilStats, 'SomaPtID') %if soma pt is set in Imaris
    SkelFiner.FilStats.SomaPtID = SomaPtID-1; %-1 to bring it back to Imaris format.
    SkelFiner.FilStats.SomaPtXYZ = SomaPtXYZ;
end

%clear Skel;

%calculate SegStats
SkelFiner.SegStats.Seg=cat(3, SkelFinerXYZ(SkelFinerEdges(:,1),:), SkelFinerXYZ(SkelFinerEdges(:,2),:));
SkelFiner.SegStats.Lengths =sqrt((SkelFiner.SegStats.Seg(:,1,1)-SkelFiner.SegStats.Seg(:,1,2)).^2 + ...
    (SkelFiner.SegStats.Seg(:,2,1)-SkelFiner.SegStats.Seg(:,2,2)).^2 + ...
    (SkelFiner.SegStats.Seg(:,3,1)-SkelFiner.SegStats.Seg(:,3,2)).^2);
disp('Total fine skel length is');
sum(SkelFiner.SegStats.Lengths) %total seg length must be the same between coarse and fine.
disp('Total fine skel node number is');
size(SkelFiner.FilStats.aXYZ,1)
disp('Total fine skel edge number is');
size(SkelFiner.FilStats.aEdges,1)
disp('Average fine skel edge length is');
sum(SkelFiner.SegStats.Lengths)/size(SkelFiner.FilStats.aEdges,1)


