%% Phase Settings for Multiplot_2. These are largely unnecessary at this point
% Phase Changes by 2Pi are legal

function [injParam,dataPhase] = phase_settings(dataPhase,n,in,plotType,injParam)

    if (n==3) && in(1).shot==160728011 && plotType ==1
            dataPhase(:,1) = dataPhase(:,1)+2*pi;
    end
    if in(n).shot == 160728011 && in(n).line==2&& plotType ==1
        dataPhase(11:17,1) = dataPhase(11:17,1)+2*pi;
       dataPhase(3,1) = dataPhase(3,1)+2*pi;
    end

    if in(n).shot == 160728011 && in(n).line==1&& plotType ==1
        dataPhase(12:16,1) = dataPhase(12:16,1)+2*pi;
    end
    if in(n).shot == 160728013 && in(n).line==2 && plotType ==1
        dataPhase(1:10,1) = dataPhase(1:10,1)-2*pi;
        dataPhase(:,1) = dataPhase(:,1)+2*pi;
    end
    if in(n).shot == 160728013 && in(n).line==1 && plotType ==1
        dataPhase(11,1) = dataPhase(11,1)-2*pi;
        dataPhase(12,1) = dataPhase(12,1)+2*pi;
        dataPhase(:,1) = dataPhase(:,1)+2*pi;
        %dataPhase(6,2) = dataPhase(6,2)+2*
    end
    if in(n).shot == 160728013 && in(n).line==2 && plotType ==2
        dataPhase=dataPhase;
    end
    if in(n).shot == 160728011 && in(n).line==2 && plotType ==2
        dataPhase=dataPhase;
    end
    if in(n).shot == 160728012 && in(n).line==2 && plotType ==2
        dataPhase=dataPhase+2*pi;
        dataPhase([3:8,15],2) = dataPhase([3:8,15],2)+2*pi;
    end
     if in(n).shot == 160728013 && in(n).line==1&& plotType ==1
        dataPhase(1:8,1) = dataPhase(1:8,1)-2*pi;
        dataPhase(12,1) = dataPhase(12,1)-2*pi;
    end
    if in(n).shot == 160728021 && in(n).line==1
        dataPhase = dataPhase+pi;
        dataPhase(:,1) = dataPhase(:,1)-2*pi;
    end
    if in(n).shot == 160728024 && in(n).line==1
        dataPhase = dataPhase-pi;

    end
    if in(n).shot == 160728024 && in(n).line==2
        dataPhase(:,2) = dataPhase(:,2)-2*pi;

    end
    if in(n).shot == 160728024 && in(n).line==2
        dataPhase = dataPhase+2*pi;
    end
    if in(n).shot == 160728012 && in(n).line==2
        %dataPhase(:,2) = dataPhase(:,2)+2*pi;
    end
    if in(n).shot == 8160609010
        dataPhase(1:6,2) = dataPhase(1:6,2)-2*pi;
        dataPhase(24:33,2) = dataPhase(24:33,2)-2*pi;
        dataPhase(24:36,1) = dataPhase(24:36,1)-2*pi;
    end
    if in(n).shot == 160525016
        dataPhase(6:7,1) = dataPhase(6:7,1) -2*pi;
    end
    if in(n).shot == 160525017
        dataPhase(1:6,1) = dataPhase(1:6,1) -2*pi;
        dataPhase(1:2,2) = dataPhase(1:2,2) -2*pi;
    end
    if in(n).shot == 160525018
       dataPhase([5,6],2) = dataPhase([5,6],2) -2*pi;
    end
    if in(n).shot == 160525019
       dataPhase(5:6,2) = dataPhase(5:6,2) -2*pi;
    end
    if any([in.shot] == 160728021)
        injParam(3)=injParam(3)-2*pi;
    end

    % Nimrod HITSI
    if in(n).shot == 8129499
       dataPhase(1:19,:) = dataPhase(1:19,:) +2*pi;
       %dataPhase(:,2) = dataPhase(:,2) -2*pi;
       dataPhase(16:19,2) = dataPhase(16:19,2) -2*pi;
    end
    if in(n).shot==151217024
        dataPhase = dataPhase-2*pi;
        dataPhase(15,2)=dataPhase(15,2)-2*pi;
    end
    if in(n).shot==151217025
        dataPhase(15:17,2) = dataPhase(15:17,2)-2*pi;
    end
    % Low Perf HitSI3
    if in(n).shot == 151217026
        dataPhase(4:12,1)=dataPhase(4:12,1)-2*pi;
        dataPhase(7:17,2)=dataPhase(7:17,2)-2*pi;
    end
    
end
