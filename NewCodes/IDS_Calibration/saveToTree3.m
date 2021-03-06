function saveToTree3(saveToShots, stt, CAL_LAMBDA, LAMBDA, VOLTAGE, MASS, PEAKS, REL_INT, PIX_SP, IMPACTS, hitsi3)
% SAVE TO TREE 3, UPDATED TO OBJECT ORIENTED JAVA FOR MDSPLUS 
%% GET THE MDSPLUS DIRECTORY STUFF TOGETHER, SO WE CAN IMPORT IT
import MDSplus.*
for n = 1:length(saveToShots)
    
    shot = saveToShots(n);
    
    disp(['Editing ANALYSIS Tree for shot ' num2str(shot)]);
    
    % Open Tree for Edit JAVA
    if(~hitsi3)
        HitTree = Tree('analysis',shot,'EDIT');
    else
        HitTree = Tree('analysis3',shot,'EDIT');
    end
    
    % Navigate to SPECTROSCOPY node JAVA
    %Try is there in the even that you are looking at shots
    %14062501 - 14070917
    try
        HitTree.setDefault(HitTree.getNode('.spectroscopy.ids')) % I ASSUME THIS IS WHAT THE IDS NODE IS CALLED
    catch
        HitTree.addNode('.SPECTROSCOPY', 'STRUCTURE')
        HitTree.addNode('.SPECTROSCOPY.IDS', 'STRUCTURE')
        HitTree.setDefault(HitTree.getNode('.spectroscopy.ids'))
    end
    
    
    % Delete any Legacy Nodes JAVA
%     try HitTree.deleteNode('IDS_CENTER');end;
%     try HitTree.deleteNode('IDS_FWHM');end;
%     try HitTree.deleteNode('IDS_BOUNDS');end;
%     try HitTree.deleteNode('IDS_STDEV');end;
%     try HitTree.deleteNode('ION_MASS');end;
%     try HitTree.deleteNode('LINE_LAMDA');end;
%     try HitTree.deleteNode('CALC_IDS');end;
   
   
    if stt.CAL_LAMBDA
        CleanNode(HitTree,'CAL_LAMBDA');
    end
    
    if stt.LAMBDA
        CleanNode(HitTree,'IDS_LAMBDA');
    end
    
    if stt.VOLTAGE
        CleanNode(HitTree,'IDS_VOLTAGE');
    end
    
    if stt.MASS
        CleanNode(HitTree,'IDS_MASS');
    end

    if stt.PEAKS
        CleanNode(HitTree,'IDS_PEAKS');
    end
    
    if stt.REL_INT
        CleanNode(HitTree,'IDS_REL_INT');
    end

    if stt.PIX_SP
        CleanNode(HitTree,'IDS_PIX_SP');
    end

    if stt.IMPACTS
        CleanNode(HitTree,'IDS_IMPACTS');
    end
    %JAVA WRITE,CLOSE, AND OPEN AGAIN JUST CAUSE.
    HitTree.write()
    try % no idea why this sometimes doesnt work
        HitTree.close()
    catch
        mdsclose;
    end
    

      
    %%%%%%%%%%%%%%
    % Fill in data ( reopen to prevent further tree modification 
    if(~hitsi3)
        HitTree = Tree('analysis', shot);
    else
        HitTree = Tree('analysis3', shot);
    end
    % spectroscopy will for sure be here by now.
    HitTree.setDefault(HitTree.getNode('.spectroscopy.ids'));
    
    if stt.CAL_LAMBDA
        AddData(HitTree,'CAL_LAMBDA',CAL_LAMBDA);
    end
    
    if stt.LAMBDA
        AddData(HitTree,'IDS_LAMBDA',LAMBDA);
    end
    
    if stt.VOLTAGE
        AddData(HitTree,'IDS_VOLTAGE',VOLTAGE);
    end
    
    if stt.MASS
        AddData(HitTree,'IDS_MASS',MASS);
    end
    
    if stt.PEAKS
        AddData(HitTree,'IDS_PEAKS',PEAKS);
    end
    
    if stt.REL_INT
        AddData(HitTree,'IDS_REL_INT',REL_INT);
    end
    
    if stt.PIX_SP
        AddData(HitTree,'IDS_PIX_SP',PIX_SP);
    end
    
    if stt.IMPACTS
        AddData(HitTree,'IDS_IMPACTS',IMPACTS);
    end
    
    try
        HitTree.close()
    catch
        mdsclose;
    end
    
end

% mdsdisconnect;
end
    
function CleanNode(HitTree, Tag)
    disp(['Clearing ' Tag]);
    try 
        current_node = HitTree.getNode(Tag);
%         HitTree.deleteNode(Tag);
    catch
%         disp(['NO NODE: ' Tag]);
        HitTree.addNode(Tag, 'NUMERIC')
        current_node = HitTree.getNode(Tag);
        current_node.addTag(Tag)
    end
%     HitTree.addNode(Tag, 'NUMERIC');
%     current_node = HitTree.getNode(Tag);
%     current_node.addTag(Tag);
end
    


function AddData(HitTree, Tag, Data)
        disp(['Modifying ' Tag]);
        current_node = HitTree.getNode(Tag);
        % Make REAL Sure we dont write NaN to tree
        if any(isnan(Data))
            error(['WARNING: ATTEMPT TO WRITE NaN TO ' Tag '. ABORTING']);
        end
        current_node.putData(MDSarg(Data))
        % MDSarg required for multidimensional arrays
        % to read: NATIVEvalue(node.getData())
end
    
    
    
    
    
    
    
    
    