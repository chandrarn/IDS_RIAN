function sys = ss(varargin)
%SS  Creates state-space model or converts model to state space.
%       
%   SYS = SS(A,B,C,D) creates a SS object SYS representing the 
%   continuous-time state-space model
%         dx/dt = Ax(t) + Bu(t)
%          y(t) = Cx(t) + Du(t)
%   You can set D=0 to mean the zero matrix of appropriate dimensions.
%   If one or more of the matrices A,B,C,D have uncertainty, SS returns 
%   an uncertain state-space (USS) model (Robust Control Toolbox only).
%
%   SYS = SS(A,B,C,D,Ts) creates a discrete-time state-space model with  
%   sample time Ts (set Ts=-1 if the sample time is undetermined).
%
%   SYS = SS creates an empty SS object.
%   SYS = SS(D) specifies a static gain matrix D.
%
%   In all cases above, the input list can be followed by pairs
%      'PropertyName1', PropertyValue1, ...
%   that set various model properties (see LTIPROPS for details).
%
%   You can create arrays of state-space models by using ND arrays for
%   A,B,C,D.  The first two dimensions of A,B,C,D define the number
%   of states, inputs, and outputs, while the remaining dimensions 
%   specify the array sizes.  For example, 
%      sys = ss(rand(2,2,3,4),[2;1],[1 1],0)
%   creates a 3x4 array of SISO state-space models.  You can also use 
%   indexed assignment and STACK to build SS arrays:
%      sys = ss(zeros(1,1,2))     % create 2x1 array of SISO models
%      sys(:,:,1) = rss(2)        % assign 1st model
%      sys(:,:,2) = ss(-1)        % assign 2nd model
%      sys = stack(1,sys,rss(5))  % add 3rd model to array
%
%   SYS = SS(SYS) converts a system SYS of arbitrary class to state 
%   space, that is, computes a state-space realization of SYS.
%
%   SYS = SS(SYS,'min') computes a minimal realization of SYS.
%
%   See also LTIMODELS, LTIPROPS, DSS, DELAYSS, RSS, DRSS, LTI/SSDATA, 
%            TF, ZPK, FRD.

%   Author(s): A. Potvin, 3-94, P. Gahinet, 4-1-96
%   Copyright 1986-2006 The MathWorks, Inc.
%   $Revision: 1.28.4.8 $  $Date: 2006/12/27 20:32:28 $

% Add static method to be included for compiler
%#function ltipack.utValidateTs

ni = nargin;
if ni>0 && isa(varargin{1},'ss'),
   % Quick exit for SS(SYS) with SYS of class SS
   if ni==1
      sys = varargin{1};
   elseif ni==2 && strncmpi(varargin{2},'min',min(3,length(varargin{2}))),
      sys = minreal(varargin{1},[],false);
   else
      error('Use SET to modify the properties of SS objects.');
   end
   return
end

% Define class-specific properties
% RE: a,b,c,d,e,InternalDelay,StateName all virtual
% superiorto('tf','zpk','double')
% inferiorto('frd')
sys = struct;

% Dissect input list
DataInputs = 0;
LtiInput = 0;
PVStart = ni+1;
for ct=1:ni
   nextarg = varargin{ct};
   if isa(nextarg,'struct') || isa(nextarg,'lti')
      % LTI settings inherited from other model
      LtiInput = ct;   PVStart = ct+1;   break
   elseif ischar(nextarg)
      PVStart = ct;   break
   else
      DataInputs = DataInputs+1;
   end
end

% Handle bad calls
if PVStart==1,
   if ni==1,
     % Bad conversion
      error('Conversion from string to ss is not possible.')
  elseif ni>0
      error('First input must contain numerical data.');
   end
elseif DataInputs>5 || (DataInputs==5 && LtiInput),
   error('Too many numerical inputs.');
end

% Create object
try
   % Process numerical data 
   switch DataInputs,
      case 0   
         if ni, 
            error('Too many LTI arguments or missing numerical data.'); 
         else
            % Empty model
            a = [];  b = [];  c = [];  d = [];
         end
         
      case 1
         % Gain matrix
         a = [];  b = [];  c = [];  
         d = checkMatrixData(varargin{1},'D');
         
      case 2
         error('Undefined C and D matrices.');
         
      case 3
         error('Undefined D matrix.');
         
      otherwise
         % A,B,C,D specified: validate data
         a = checkMatrixData(varargin{1},'A');
         b = checkMatrixData(varargin{2},'B');
         c = checkMatrixData(varargin{3},'C');
         d = checkMatrixData(varargin{4},'D');
   end
   
   % Sample time
   if DataInputs==5
      % Discrete SS
      try 
         Ts = ltipack.utValidateTs(varargin{5});
      catch
         rethrow(lasterror);
      end
   else
      Ts = 0;
   end
   
   % Delays and array size
   if ni>0
      [d,Ny,Nu] = getIOSize(d,b,c);
      ArraySize = getLTIArraySize(a,b,c,d);
      if isempty(ArraySize)
         error('Arrays A,B,C,D,E have incompatible dimensions.')
      end
   else
      Ny = 0;  Nu = 0;  ArraySize = [1 1];
   end
   Delay = struct(...
      'Input',zeros(Nu,1),...
      'Output',zeros(Ny,1),...
      'IO',[],...
      'Internal',zeros(0,1));
   Nsys = prod(ArraySize);
   
   % Create @ssdata object array
   % RE: Inlined for optimal speed
   sn(1:size(a,1),1) = {''};
   if Nsys==1
      % Sparse does not like ND indexing...
      Data = ltipack.ssdata;
      Data.a = a;
      Data.b = b;
      Data.c = c;
      Data.d = d;
      Data.Delay = Delay;
      Data.Ts = Ts;
      Data.StateName = sn;
   else
      Data = handle(zeros(ArraySize));
      for ct=1:Nsys
         D = ltipack.ssdata;
         D.a = a(:,:,min(ct,end));
         D.b = b(:,:,min(ct,end));
         D.c = c(:,:,min(ct,end));
         D.d = d(:,:,min(ct,end));
         D.Ts = Ts;
         D.Delay = Delay;
         D.StateName = sn;
         Data(ct) = D;
      end
      Data = reshape(Data,ArraySize); % UDDREVISIT: handle(zeros(...)) skips singleton dims
   end
   
   % Create LTI parent
   L = lti(Ny,Nu);
   L = setPrivateData(L,Data);
   
   % SS is a subclass of LTI
   sys = class(sys,'ss',L);
   
   if ni>0 
      % RE: Skip when just constructing empty instance for efficiency
      % LTI properties inherited from other system
      if LtiInput,
         Settings = varargin{LtiInput};
         if ~isa(Settings,'struct')
            Settings = utGetSettings(Settings);
         end
         sys = utCopySettings(sys,Settings);
      end
            
      % Property/Value pairs
      [pvpairs,iodsettings] = LocalCheckDelaySettings(varargin(:,PVStart:ni));
      if ~isempty(pvpairs)
         sys = utFastSet(sys,pvpairs{:});
      end
      
      % Consistency checks 
      sys = utCheckSystem(sys);
      
      % I/O delay settings. Must be done after data checks to prevent errors in 
      % setIODelay when A,B,C are not properly formatted (e.g., A=B=C=[] and D=1)
      if ~isempty(iodsettings)
         sys = utFastSet(sys,iodsettings{:});
      end     
   end

catch
   rethrow(lasterror)
end


%--------------------- Local Functions --------------------------------

function [pvp,iodpvp] = LocalCheckDelaySettings(pvp)
% Pulls out ioDelay settings, throws error when setting InternalDelays
if any(strncmpi(pvp(1:2:end),'int',3))
   error('%s\n%s','Internal delays cannot be specified directly. Use DELAYSS or interconnection',...
   'functions like PARALLEL or FEEDBACK to build models with internal delays.')
end
idx = find(strncmpi(pvp(1:2:end),'io',2));
iodpvp = cell(1,0);
for ct=length(idx):-1:1
   k = 2*idx(ct)-1;
   iodpvp = [iodpvp , pvp(k:min(k+1:end))];
   pvp = [pvp(1:k-1) , pvp(k+2:end)];
end
