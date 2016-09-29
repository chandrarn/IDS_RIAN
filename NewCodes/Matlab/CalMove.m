% move calibration data from one shot in the tree to another
% Sometimes MDSplus gets upset about NATIVEvalue
%clear all;
% refShot = 151124035;
% newShot = [151217016:151217026];
% refShot = 160525019
% newShot = [ 160601011:160601026];
refShot = 160601011
newShot = [ 160802008:160802013,160802015:160802020, 160803011:160803012];

postHoc = 0; % replace info for already corrected data. 

addpath('T:\IDS\Data Analysis\');
addpath('T:\IDS\Calibration\');
addpath('T:\IDS\Data Repository');


if ~postHoc
    %get data from ref shot in tree
    import MDSplus.*
    Conn = Connection('landau.hit');
    Conn.openTree('hitsi3',refShot);

    param.peaks = double(NATIVEvalue(Conn.get('\IDS_PEAKS').getDoubleArray)); 
    param.IonMass = double((Conn.get('\IDS_MASS')));
    param.LineLam = double((Conn.get('\IDS_LAMBDA')));
    param.PIX_SP = double((Conn.get('\IDS_PIX_SP')));
    param.REL_INT = double((Conn.get('\IDS_REL_INT')));
    param.impacts = double((Conn.get('\IDS_IMPACTS')));
    param.CalLam = double((Conn.get('\CAL_LAMBDA')));
    %[param,options] = loadParams(refShot,[1],1,1);
    param.IonMass;
    stt.CAL_LAMBDA = 1;
    stt.LAMBDA = 1;
    stt.VOLTAGE = 1;
    stt.MASS = 1;
    stt.PEAKS = 1;
    stt.REL_INT = 1;
    stt.PIX_SP = 1;
    stt.IMPACTS = 1;

    %%  make changes
    % load('T:\RChandra\NewCodes\Geometry\impacts5.mat');
    % param.impacts = impacts(param.peaks(:,1));
    % param.LineLam = 1e-9 * [465.147, 465.025 464.742];
    % param.IonMass = [12,16,12];
    % param.CalLam = 1e-9 * 656.28;
    %param.PIX_SP = param.PIX_SP *.7615; % fix pixsp
     %param.LineLam = [464.7418, 464.91348,465.0246].*1e-9;
     % Temp fix, Calibration was railing sigx,y,
     param.peaks=reshape(param.peaks,[],5);
%      param.peaks(:,5) = [1.29566603799762;1.34962188121872;1.35706166437118;1.39363669896754;1.35423209447031;1.33867588925782;1.29612925311774;1.39871160888804;1.41899467674756;1.43458558040884;1.43537682109607;1.39882047718025;1.37020234124237;1.40792642263123;1.42216156848763;1.45328522212848;1.48482393612287;1.45683031920635;1.46505395988887;1.52574976273244;1.51177815664657;1.54729877463436;1.52572731612651;1.53516299558698;1.46471033101584;1.53950680012712;1.54270411041929;1.57405563128865;1.62287522238338;1.56451812110863;1.56439475295522;1.56799758883217;1.52795755859407;1.50254278826412;1.51196798181507;1.56095916512397;1.56248600862580;1.59876413684088;1.59970146929494;1.57804885392595;1.59391292757677;1.55116301465120;1.52571592251769;1.57555323811128;1.61013151576554;1.62145355726744;1.65056133662732;1.61522084747414;1.62131206994450;1.59256933023350;1.60230224134206;1.64409281239272;1.63268112778468;1.61034313507482;1.67903811538727];
%      param.peaks(:,4) = [1.33238096204252;1.29946003910481;1.27407798571707;1.30715792452569;1.41713658845607;1.43541605934074;1.55654860339273;1.28521749550807;1.22707935463328;1.21321471994077;1.22257825870154;1.32319624926938;1.37459462306191;1.25835599827057;1.22053674693480;1.20743551579751;1.18817719457109;1.21715287861811;1.19355862594900;1.16022974311681;1.15684953765502;1.14076374842497;1.15503580471834;1.16675779503446;1.23870235150247;1.14111390619505;1.15252764794329;1.12200840732130;0.942647309435541;1.17495676915540;1.18426443273999;1.19004939917720;1.17522104354076;1.23142419354747;1.23958237268361;1.17899780192186;1.17621409895439;1.12651501077242;1.11256132739080;1.13473150320607;1.11824479504028;1.14874379920468;1.23604718587386;1.12709234215238;1.10530702847405;1.06700912231444;1.05700253041841;1.05478811710535;1.06704363094647;1.09125623578023;1.09797250144703;1.05321614589731;1.04786123400952;1.11499934390288;1.03618031963642];
%       param.peaks=reshape(param.peaks,5*55,1);
     % param.peaks(:,4) =  MAY WISH TO FIX THIS AS WELL;
     % Fix callibration error, sigy was hitting ub, 151217026
     %param.impacts = IMPACTS';
     % Fix CalLamb
     %param.CalLam = 4.65025E-7; % TEMP FIX 151217026, it was set to 656nm
     % For 16051800 and on
%      param.LineLam = [4.6418104e-07; 4.64741788164247e-07;4.64913483710916e-07;4.65024612594789e-07];
    %param.IonMass = [16;12;16;12];
    %param.Inst_Temp(:,1:4) = [param.Inst_Temp(:,2) ,param.Inst_Temp(1:3)]; 
    % fix callibration error for 160518 and 19
%     figure; plot(param.peaks(:,2),param.peaks(:,3));
    % for 160518,19
    %param.peaks(30:end,3) = param.peaks(30:end,3) -3*0.202768609636075 -0.195256879282986 -0.001045211788516;% spectrometer bumped between shots
    % for 151217
    %param.peaks(30:end,3) = param.peaks(30:end,3) +0.241756163349005-0.008248785763696;% spectrometer bumped between shots

    %    figure; plot(param.peaks(:,2),param.peaks(:,3));
    for i = newShot
        saveToTree3(i,stt,param.CalLam,param.LineLam,0,param.IonMass,param.peaks,param.REL_INT,param.PIX_SP,param.impacts);
    end
else
    load(['dat' num2str(refShot) '10.mat']);
    param = dat(1).param;
    for i = newShot
        try
        
        load(['dat' num2str(i) '10.mat']);
        display(['Saving shot: ' num2str(i)]);
        %dat(1).param = param;
        dat(1).param.LineLam = param.LineLam;
        %dat(1).impacts = param.impacts;
        save(['T:\IDS\Data Repository\dat' num2str(i) '10.mat'],'dat');
        catch 
            disp([ 'Unable to open shot ' num2str(i)]);
        end
    end
end