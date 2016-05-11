function[Status] = Stat(TPN)
%% Checks and retrives cell Status when given folder name

  %%Get or make Status
    clear Status
    if ~strcmp(TPN(1,size(TPN,2)),'\')
    TPN=[TPN '\']; end

    if exist([TPN 'Status.mat'])
       load([TPN 'Status.mat']);
    else 
       Status.Fix=0;
       Status.Ra=0;
       Status.Ma=0;
       Status.Dr=0;
       Status.Running=0;
       Status.FD=0;
       Status.Sk=0;
       Status.Results=0;
       save([TPN 'Status.mat'],'Status')
    end
   
    
        %%Check Status
       if exist([TPN 'data\BigFilled.mat']), Status.FD=1;
       else Status.FD=0; end
       
       if exist([TPN 'data\AllSeg.mat']), Status.Sk=1;
       else Status.Sk=0; end       
       
       if exist([TPN 'data\MaskedAt.mat']), Status.Ma=1;
       else Status.Ma=0; end   
       %% and because Im an idiot
       if exist([TPN 'data\DotStatsBackup.mat']),Status.Ma2=1;
       else Status.Ma2=0; end
                    
       if exist([TPN 'data\DotStats.mat']), Status.Ra=1;
       else Status.Ra=0; end  
       
       if exist([TPN 'data\Results.mat']), Status.Results=1;
       else Status.DD=0; end   

       if exist([TPN 'dataFix']), Status.Fix =1;
       else Status.Fix=0; end
   
       if exist([TPN 'images\Combo.tif']), Status.Dr =1;
       else Status.Dr =0; end
       
       Status.Fields=0;
       Status.Fields=fieldnames(Status);
       save([TPN 'Status.mat'],'Status')
        
          
end %run all data
