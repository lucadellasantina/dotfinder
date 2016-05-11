function[Cell,Results,Status,Combo]=LoadData(DPN)
        Cell=0; Results=0; Combo=0; Status=0;
        [DPN '\images\Combo.tif']
        if exist([DPN '\Cell.mat']),load([DPN '\Cell.mat']),end
        if exist([DPN '\data\Results.mat']),load([DPN '\data\Results.mat']),end
        Status=Stat(DPN); 
        
        if exist([DPN '\images\Combo.tif'])
            Combo=imread([DPN '\images\Combo.tif']);
        end
        
       
         