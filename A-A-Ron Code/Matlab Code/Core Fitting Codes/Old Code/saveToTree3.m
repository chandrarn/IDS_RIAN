function saveToTree3(saveToShots, stt, CAL_LAMBDA, LAMBDA, VOLTAGE, MASS, PEAKS, REL_INT, PIX_SP, IMPACTS)
% SAVE TO TREE 3, UPDATED TO OBJECT ORIENTED JAVA FOR MDSPLUS APR 2014 


import MDSplus.*
global HitTree

for n = 1:length(saveToShots)
    
    shot = saveToShots(n);
    
    disp(['Editing ANALYSIS Tree for shot ' num2str(shot)]);
    
    %Open Tree for Edit JAVA 
    %PICK FEB 8, 2014, TO DECIDE ANALYSIS OR ANALYSIS3
    if(shot>140208000)
        HitTree = Tree('analysis3',shot,'EDIT');
    else
        HitTree = Tree('analysis',shot,'EDIT');
    end
    
     
    % Navigate to SPECTROSCOPY node 
    HitTree.setDefault(HitTree.getNode('.spectroscopy.ids')); % I ASSUME THIS IS WHAT THE IDS NODE IS CALLED
   
    
    % Delete any Legacy Nodes JAVA (Try: otherwise MDSplus faults)
    try HitTree.deleteNode('IDS_CENTER'); end
    try HitTree.deleteNode('IDS_FWHM'); end
    try HitTree.deleteNode('IDS_BOUNDS'); end
    try HitTree.deleteNode('IDS_STDEV'); end
    try HitTree.deleteNode('ION_MASS'); end
    try HitTree.deleteNode('LINE_LAMBDA'); end
    try HitTree.deleteNode('CALC_IDS'); end
   
    
    if stt.CAL_LAMBDA
         EditNode('CAL_LAMBDA');
    end  
    
    if stt.LAMBDA
         EditNode('IDS_LAMBDA');
    end
     
     if stt.VOLTAGE
         EditNode('IDS_VOLTAGE');
     end
     
     if stt.MASS
          EditNode('IDS_MASS');
     end
     
     if stt.PEAKS
          EditNode('IDS_PEAKS');
     end
     
     if stt.REL_INT
         EditNode('IDS_REL_INT');
     end

     if stt.PIX_SP
         EditNode('IDS_PIX_SP');

     end
 
     if stt.IMPACTS
         EditNode('IDS_IMPACTS');

     end
     
    % WRITE AND THEN REOPEN: ONLY ONE PERSON CAN EDIT AT TIME
    HitTree.write();  
    HitTree = Tree('analysis',shot,'NORMAL');
    HitTree.setDefault(HitTree.getNode('.spectroscopy.ids')); 
    
    %%%%%%%%%%%%%%
    % Fill in data
    
    if stt.CAL_LAMBDA
        PutData('CAL_LAMBDA',CAL_LAMBDA);
        % May need: current_node.putData(MDSargs(CAL_LAMBDA));
        % or current_node.putData(Float64(CAL_LAMBDA));
    end
    
     if stt.LAMBDA
         PutData('IDS_LAMBDA',LAMBDA);
     end
    
     if stt.VOLTAGE
         PutData('IDS_VOLTAGE',VOLTAGE);
     end
     
     if stt.MASS
         PutData('IDS_MASS',MASS);
     end
     
     if stt.PEAKS
         PutData('IDS_PEAKS',PEAKS);
     end
     
     if stt.REL_INT
         PutData('IDS_REL_INT',REL_INT);
     end
     
     if stt.PIX_SP
         PutData('IDS_PIX_SP',PIX_SP);
     end
     
     if stt.IMPACTS
         PutData('IDS_IMPACTS',IMPACTS);
     end
     
    
end

end
    

function EditNode(Node)
        global HitTree
        disp(['Modifying ',Node]);
         try HitTree.deleteNode(Node); end
         HitTree.addNode(Node,'NUMERIC')
         current_node=HitTree.getNode(Node);
         current_node.addTag(Node);
end

function PutData(Node,Data)
        global HitTree 
        disp(['Saving data to: ',Node]);
        current_node=HitTree.getNode(Node);
        current_node.putData(MDSarg(Data));
end
    
    
    
    
    
    
    
    
    
    