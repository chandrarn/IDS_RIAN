%% Modify calibration data (in param structure) in an existing shot
% Sometimes MDSplus gets upset about NATIVEvalue
%clear all;

import MDSplus.*

% What shot to pull reference parameters from
% refShot = 151124035;
% newShot = [151217016:151217026];
% refShot = 160525019
% newShot = [ 160601011:160601026];
%refShot = [170518020];160728013;

% Shots to apply calibration to
newShot = [-1,170518025,170518028,170518030];

postHoc = 4; 
% 0 All movement within tree
% 1 all movement outside tree
% 2 new data in workspace
% 3 old data outside tree
% 4 move each shot from its own dat into tree

addpath('T:\IDS\Data Analysis\');
addpath('T:\IDS\Calibration\');
addpath('T:\IDS\Data Repository');


if postHoc==0 || postHoc == 3 || postHoc == 2
    
    if postHoc==0 %get data from ref shot in tree
        
        Conn = Connection('landau.hit');
        Conn.openTree('hitsi3',refShot);

        param.peaks = double(NATIVEvalue(Conn.get('\IDS_PEAKS').getDoubleArray)); 
        param.IonMass = double((Conn.get('\IDS_MASS')));
        param.LineLam = double((Conn.get('\IDS_LAMBDA')));
        param.PIX_SP = double((Conn.get('\IDS_PIX_SP')));
        param.REL_INT = double((Conn.get('\IDS_REL_INT')));
        param.impacts = double((Conn.get('\IDS_IMPACTS')));
        param.CalLam = double((Conn.get('\CAL_LAMBDA')));
    elseif postHoc==3 % Get data from an existing dat file
        load(['dat' num2str(refShot) '10.mat']);
        param = dat(1).param;
    elseif postHoc==2 % Get data from a variable in the workspace
        try;param.peaks = PEAKS;catch;param.PEAKS = [0];end
        try;param.IonMass = MASS;catch;param.IonMass = [0];end
        try;param.LineLam = LAMBDA;catch;param.LineLam = [0];end
        try;param.PIX_SP = PIX_SP;catch;param.PIX_SP = [0];end
        try;param.REL_INT = REL_INT;catch;param.REL_INT = [0];end
        try;param.impacts = impacts;catch;param.impacts = [0];end
        try;param.CalLam = CAL_LAMBDA;catch;param.CAL_LAMBDA = [0];end
    end
    %[param,options] = loadParams(refShot,[1],1,1);
    param.IonMass;
    stt.CAL_LAMBDA = 0;
    stt.LAMBDA = 0;
    stt.VOLTAGE = 0;
    stt.MASS = 0;
    stt.PEAKS = 0;
    stt.REL_INT = 0;
    stt.PIX_SP = 1;
    stt.IMPACTS = 0;

    %%  make changes
    % load('T:\RChandra\NewCodes\Geometry\impacts5.mat');
    % param.impacts = impacts(param.peaks(:,1));
    % param.LineLam = 1e-9 * [465.147, 465.025 464.742];
     %param.IonMass = [12,16,12];
    % param.CalLam = 1e-9 * 656.28;
    %param.PIX_SP = param.PIX_SP *.7615; % fix pixsp
     %param.LineLam = [464.7418, 464.91348,465.0246].*1e-9;
     % Temp fix, Calibration was railing sigx,y,
     %param.peaks=reshape(param.peaks,[],5);
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
    %param.IonMass = [12;16;12;16];
    %param.Inst_Temp(:,1:4) = [param.Inst_Temp(:,2) ,param.Inst_Temp(1:3)]; 
    % fix callibration error for 160518 and 19
%     figure; plot(param.peaks(:,2),param.peaks(:,3));
    % for 160518,19
    %param.peaks(30:end,3) = param.peaks(30:end,3) -3*0.202768609636075 -0.195256879282986 -0.001045211788516;% spectrometer bumped between shots
    % for 151217
    %param.peaks(30:end,3) = param.peaks(30:end,3) +0.241756163349005-0.008248785763696;% spectrometer bumped between shots
    % for [ 161018009:161018027]
    %param.impacts=[55.1426895042383,54.7460283605416,54.2042712159711,53.5188539152350,51.7246811561919,50.6206808792197,49.3825182058684,48.0134746975174,46.5171787952431,44.8975962031896,43.1590193780746,41.3060561526868,39.3436175235280,37.2769046349659,35.1113949943933,32.8528279549308,30.5071895041457,28.0806963991051,25.5797796898091,23.0110676746732,20.3813683332339,17.6976512826375,14.9670293057320,12.1967394997207,9.39412409533987,6.56661099739543,3.72169409823381,0.866913416323006,-1.99016488741679,-4.84196856256908,-7.68093933831615,-10.4995529554835,-13.2903391084416,-16.0459012440117,-18.7589361649027,-55.1426895042383,-54.7460283605416,-54.2042712159711,-53.5188539152350,-52.6915930523813,-51.7246811561919,-50.6206808792197,-49.3825182058684,-48.0134746975174,-46.5171787952431,-44.8975962031896,-43.1590193780746,-41.3060561526868,-39.3436175235280,-37.2769046349659,-35.1113949943933,-32.8528279549308,-28.0806963991051,-25.5797796898091,-23.0110676746732,-20.3813683332339,-17.6976512826375,-14.9670293057320,-12.1967394997207,-9.39412409533987,-6.56661099739543,-3.72169409823381,-0.866913416323006,1.99016488741679,4.84196856256908,7.68093933831615,10.4995529554835,13.2903391084416,16.0459012440117,18.7589361649027];
    %param.peaks=reshape(param.peaks,55,5);
    %param.peaks(:,3)=[64.5716336319537;64.7152980165642;64.2565331013773;64.3445402214184;64.1250179998928;64.0890338259467;64.1463385146612;64.2174714123468;63.9810086118457;64.0595592833065;64.0454113041528;63.9832502318573;63.9919644797281;63.8377753625362;63.8359635374366;63.7285859978438;63.7272664521733;63.6572336478220;63.6182377553779;63.6114257110770;63.4878595184938;63.4472873595804;63.3180277773227;63.3407453050264;63.4449343679779;63.2957371643466;63.3186866939331;63.1711350645435;63.2771451665117;63.5250198380493;63.2142286840405;63.3895954313549;63.5569488414813;63.4186655089664;64.5098465240965;64.2108206158928;62.5493763726233;63.0240634814728;63.3577060082475;62.8930852872995;62.3332787562951;62.4663542826716;62.4115543338672;62.4592504373161;62.3544579913494;62.2606300852417;62.4598056504478;62.4600286788039;62.1960477622280;62.0840092227450;62.2965089895570;62.2411779184083;62.1388031211641;62.1092184966401;62.1560378112728;62.0951687329900;62.0269058427853;62.0108595031086;61.9417243765801;62.0495196576660;62.0493435396346;61.9393956855534;61.9340556575973;61.8560231970375;61.8179316995842;61.6073533323180;61.7263135768224;61.7206944689361;61.7976898117878;62.0337779527818];
    %param.PIX_SP=param.PIX_SP.*.97; %slight adjustment
    %    figure; plot(param.peaks(:,2),param.peaks(:,3));
    for i = newShot
        saveToTree3(i,stt,param.CalLam,param.LineLam,0,param.IonMass,param.peaks,param.REL_INT,param.PIX_SP,param.impacts);
    end
elseif postHoc == 2
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
elseif postHoc == 4
    stt.CAL_LAMBDA = 0;
    stt.LAMBDA = 0;
    stt.VOLTAGE = 0;
    stt.MASS = 0;
    stt.PEAKS = 0;
    stt.REL_INT = 0;
    stt.PIX_SP = 1;
    stt.IMPACTS = 0;
    
    files=dir('T:\IDS\Data Repository\dat16*.mat');
    for i = 1:length(files)
        shot=files(i).name;
        if length(shot)==18
            load(['T:\IDS\Data Repository\' shot ]);
            if mod(i,10)==0;display(['File ' num2str(i) '/' num2str(length(files))]);end
            param=dat(1).param;
            saveToTree3(str2num(shot(4:12)),stt,param.CalLam,param.LineLam,0,...
                param.IonMass,param.peaks,param.REL_INT,param.PIX_SP,param.impacts);
        end
    end
    
end