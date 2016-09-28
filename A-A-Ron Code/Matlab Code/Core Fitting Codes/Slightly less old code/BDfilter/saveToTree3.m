function saveToTree3(saveToShots, stt, CAL_LAMBDA, LAMBDA, VOLTAGE, MASS, PEAKS, REL_INT, PIX_SP, IMPACTS)
% SAVE TO TREE 3, UPDATED TO OBJECT ORIENTED JAVA FOR MDSPLUS 
%% GET THE MDSPLUS DIRECTORY STUFF TOGETHER, SO WE CAN IMPORT IT
%%
import MDSplus.*
for n = 1:length(saveToShots)
    
    shot = saveToShots(n);
    
%     mdsconnect('landau.hit')
%     output = [];
    disp(['Editing ANALYSIS Tree for shot ' num2str(shot)]);
    
    %Open Tree for Edit JAVA
    HitTree = Tree('analysis',shot,'EDIT')
%     % Open Tree for Edit
     %output = [output, mdsvalue(['tcl("edit ANALYSIS/shot=' num2str(shot) '")'])];
%     
    % Navigate to SPECTROSCOPY node JAVA
    HitTree.setDefault(HitTree.getNode('.spectroscopy.ids')); % I ASSUME THIS IS WHAT THE IDS NODE IS CALLED
%     % Navigate to SPECTROSCOPY node
%     output = [output, mdsvalue('tcl("set def .spectroscopy.ids")')];
    
    
    % Delete any Legacy Nodes JAVA
    HitTree.deleteNode('IDS_CENTER');
    HitTree.deleteNode('IDS_FWHM');
    HitTree.deleteNode('IDS_BOUNDS');
    HitTree.deleteNode('IDS_STDEV');
    HitTree.deleteNode('ION_MASS');
    HitTree.deleteNode('LINE_LAMDA');
    HitTree.deleteNode('CALC_IDS');
    
    % Delete any Legacy Nodes
%     output = [output, mdsvalue('tcl("delete node IDS_CENTER /noconfirm")')];
%     output = [output, mdsvalue('tcl("delete node IDS_FWHM /noconfirm")')];
%     output = [output, mdsvalue('tcl("delete node IDS_BOUNDS /noconfirm")')];
%     output = [output, mdsvalue('tcl("delete node IDS_STDEV /noconfirm")')];
%     output = [output, mdsvalue('tcl("delete node ION_MASS /noconfirm")')];
%     output = [output, mdsvalue('tcl("delete node LINE_LAMBDA /noconfirm")')];
%     
    
    
    %JAVA VERSION OF CAL LAMBDA
    if stt.CAL_LAMBDA
            disp('Modifying CAL_LAMBDA');
            HitTree.deleteNode('CAL_LAMBDA');
            current_node=HitTree.addNode('CAL_LAMBDA','NUMERIC');
            current_node.addTag('CAL_LAMBDA');
            %THE CURRENT_MODE THING MIGHT NOT WORK
    end     
%     if stt.CAL_LAMBDA
%         disp('Modifying CAL_LAMBDA');
%         output = [output, mdsvalue('tcl("delete node CAL_LAMBDA /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric CAL_LAMBDA")')];
%         output = [output, mdsvalue('tcl("add tag CAL_LAMBDA CAL_LAMBDA")')];
%     end
%     
%     if stt.LAMBDA
%         disp('Modifying IDS_LAMBDA');
%         output = [output, mdsvalue('tcl("delete node IDS_LAMBDA /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_LAMBDA")')];
%         output = [output, mdsvalue('tcl("add tag IDS_LAMBDA IDS_LAMBDA")')];
%     end
%     
%     if stt.VOLTAGE
%         disp('Modifying IDS_VOLTAGE');
%         output = [output, mdsvalue('tcl("delete node IDS_VOLTAGE /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_VOLTAGE")')];
%         output = [output, mdsvalue('tcl("add tag IDS_VOLTAGE IDS_VOLTAGE")')];
%     end
%     
%     if stt.MASS
%         disp('Modifying IDS_MASS');
%         output = [output, mdsvalue('tcl("delete node IDS_MASS /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_MASS")')];
%         output = [output, mdsvalue('tcl("add tag IDS_MASS IDS_MASS")')];
%     end
%     disp('Done Modifying IDS_MASS');
%     if stt.PEAKS
%         disp('Modifying IDS_PEAKS');
%         output = [output, mdsvalue('tcl("delete node IDS_PEAKS /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_PEAKS")')];
%         output = [output, mdsvalue('tcl("add tag IDS_PEAKS IDS_PEAKS")')];
%     end
%     
%     if stt.REL_INT
%         disp('Modifying IDS_REL_INT');
%         output = [output, mdsvalue('tcl("delete node IDS_REL_INT /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_REL_INT")')];
%         output = [output, mdsvalue('tcl("add tag IDS_REL_INT IDS_REL_INT")')];
%     end
% 
%     if stt.PIX_SP
%         disp('Modifying IDS_PIX_SP');
%         output = [output, mdsvalue('tcl("delete node IDS_PIX_SP /noconfirm")')]
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_PIX_SP")')]
%         output = [output, mdsvalue('tcl("add tag IDS_PIX_SP IDS_PIX_SP")')]
%     end
% 
%     if stt.IMPACTS
%         disp('Modifying IDS_IMPACTS');
%         output = [output, mdsvalue('tcl("delete node IDS_IMPACTS /noconfirm")')];
%         output = [output, mdsvalue('tcl("add node/usage=numeric IDS_IMPACTS")')];
%         output = [output, mdsvalue('tcl("add tag IDS_IMPACTS IDS_IMPACTS")')];
%     end
    %JAVA WRITE,CLOSE, AND OPEN AGAIN JUST CAUSE.
    HitTree.write();
    mdsclose;
    
    HitTree = Tree('analysis',shot);
    HitTree.setDefault(HitTree.getNode('.spectroscopy.ids'));
    
%     % Write and close
%     output = [output, mdsvalue('tcl("write")')];
%     output = [output, mdsvalue('tcl("close")')];
%     
    %%%%%%%%%%%%%%
    % Fill in data
    % I DONT THINK I NEED TO OPEN THE ANALYSIS SUBTREE A SECOND TIME
    %output = [output, mdsopen('analysis', shot)];
    
    if stt.CAL_LAMBDA
        %output = [output, mdsput('\CAL_LAMBDA', '$', CAL_LAMBDA)];
        current_node=HitTree.getNode('CAL_LAMBDA');
        current_node.putData(CAL_LAMBDA);
        % May need: current_node.putData(MDSargs(CAL_LAMBDA));
        % or current_node.putData(Float64(CAL_LAMBDA));
    end
    
%     if stt.LAMBDA
%         output = [output, mdsput('\IDS_LAMBDA', '$', LAMBDA)];
%     end
%     
%     if stt.VOLTAGE
%         output = [output, mdsput('\IDS_VOLTAGE', '$', VOLTAGE)];
%     end
%     
%     if stt.MASS
%         output = [output, mdsput('\IDS_MASS', '$', MASS)];
%     end
%     
%     if stt.PEAKS
%         output = [output, mdsput('\IDS_PEAKS', '$', PEAKS)];
%     end
%     
%     if stt.REL_INT
%         output = [output, mdsput('\IDS_REL_INT', '$', REL_INT)];
%     end
%     
%     if stt.PIX_SP
%         output = [output, mdsput('\IDS_PIX_SP', '$', PIX_SP)];
%     end
%     
%     if stt.IMPACTS
%         output = [output, mdsput('\IDS_IMPACTS', '$', IMPACTS)];
%     end
%     
%     mdsclose();
    
end

mdsdisconnect;
end
    
    
    
    
    
    
    
    
    
    
    
    
    