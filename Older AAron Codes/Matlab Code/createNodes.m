% Create Nodes in Tree

shots = [-1]; % shot(s) where nodes should be created

% For the following variables, 1 = create, 0 = no
ids_voltage = 1;

mdsconnect('landau.hit');
for n = 1:length(shots)

    output = [];
    disp(['Editing ANALYSIS Tree for shot ' num2str(shots(n))]);

    % Open Tree for Edit
    output = [output, mdsvalue(['tcl("edit ANALYSIS/shot=' num2str(shots(n)) '")'])];

    % Navigate to SPECTROSCOPY node
    output = [output, mdsvalue('tcl("set def .spectroscopy.ids")')];
    
    if ids_voltage
        output = [output, mdsvalue('tcl("delete node IDS_VOLTAGE /noconfirm")')];
        output = [output, mdsvalue('tcl("add node/usage=numeric IDS_VOLTAGE")')];
        output = [output, mdsvalue('tcl("add tag IDS_VOLTAGE IDS_VOLTAGE")')];
    end
    
    
    % Write and close
    output = [output, mdsvalue('tcl("write")')];
    output = [output, mdsvalue('tcl("close")')];
    
end