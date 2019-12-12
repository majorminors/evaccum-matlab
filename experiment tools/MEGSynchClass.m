% Interface for National Instruments PCI 6503 card
% Version 1.3
%
% DESCRIPTION
%
% N.B.: It does not monitor pulses in the background, so you have to make sure that you wait for any pulse before it comes!
%
% Properties (internal variables):
% 	IsValid						= device set and operational
%   Keys                        = set a cell of key names for emulation mode. For key names look for KbName.m.
%                                   N.B.: Requires PTB.
%                                   N.B.: Suppress passing keypresses to MATLAB during the whole experiment
%
% 	Clock						= internal clock (seconds past since the first scanner pulse or clock reset)
%
% 	Buttons						= current state of the button box
% 	LastButtonPress				= index/indices of the last button(s) pressed
% 	TimeOfLastButtonPress		= time (according to the internal clock) of the last button press (any)
%   BBoxTimeout                 = set a non-Inf value (in seconds) to wait for button press only for a limited time
%                               = set a negative value (in seconds) to wait even in case of response
%
% Methods (internal functions):
% 	MEGSynchClass   			= constructor
% 	delete						= destructor
%	ResetClock					= reset internal clock
%
%   SendTrigger(v)              = sends 'v' (integer between 0-255) as a trigger pulse
%
%	SetButtonReadoutTime(t) 	= blocks individual button readout after a button press for 't' seconds (detection of other buttons is still possible)
%	SetButtonBoxReadoutTime(t)	= blocks the whole button box readout after a button press for 't' seconds (detection of other buttons is also not possible)
%	WaitForButtonPress			= wait until a button is pressed
% 	WaitForButtonRelease		= wait until a button is released
%
% USAGE
%
% Initialise:
% 	MEG = MEGSynchClass;
% 	% MEG = MEGSynchClass(1); % for full emulation (no DAQ available)
%   % MEG = MEGSynchClass({'B1','S4';'B2','S5'}); % specify button names and maps them to S4 and S5
% Close:
% 	MEG.delete;
%
% Example for Trigger:
%   for v = 0:255
%       pause(0.1)
%       MEG.SendTrigger(v)
%   end
%
% Example for buttons:
% 	MEG.SetButtonReadoutTime(0.5);      % block individual buttons
%	% MEG.SetButtonBoxReadoutTime(0.5); % block the whole buttonbox
%   % MEG.Keys = {'f1','f2','f3','f4'}; % emulation Buttons #1-#4 with F1-F4
%	n = 0;
%   % MEG.BBoxTimeout = 1.5;            % Wait for button press for 1.5s
%   % MEG.BBoxTimeout = -1.5;           % Wait for button press for 1.5s even in case of response
%	MEG.ResetClock;
%	while n ~= 10                       % polls 10 button presses
%   	MEG.WaitForButtonPress;         % Wait for any button to be pressed
%   	% MEG.WaitForButtonRelease;     % Wait for any button to be	released
%   	% MEG.WaitForButtonPress([],2); % Wait for Button #2
%   	% MEG.WaitForButtonPress(2);    % Wait for any button for 2s (overrides MEG.BBoxTimeout only for this event)
%   	% MEG.WaitForButtonPress(-2);   % Wait for any button for 2s even in case of response (overrides MEG.BBoxTimeout only for this event)
%   	% MEG.WaitForButtonPress(2,2);  % Wait for Button #2 for 2s (overrides MEG.BBoxTimeout only for this event)
%   	% MEG.WaitForButtonPress(-2,2); % Wait for Button #2 for 2s even in case of response (overrides MEG.BBoxTimeout only for this event)
%    	n = n + 1;
%    	fprintf('At %2.3fs, ',MEG.Clock);
%    	fprintf('Button %d ',MEG.LastButtonPress);
%    	fprintf('pressed: %2.3fs\n',MEG.TimeOfLastButtonPress);
%	end
%_______________________________________________________________________
% Copyright (C) 2017 MRC CBSU Cambridge
%
% Original by Tibor Auer.
% Updates by Johan Carlin: johan.carlin@mrc-cbu.cam.ac.uk
%_______________________________________________________________________

classdef MEGSynchClass < handle
    
    
    properties
        Keys = {}
        BBoxTimeout = Inf % second (timeout for WaitForButtonPress)
    end
    
    properties (SetAccess = private)
        LastButtonPress
    end
    
    properties (Access = private)
        Map = {...
            'S3' [0 1];... % S3 - port0/line1
            'S4' [0 2];... % S4 - port0/line2
            'S5' [0 3];... % S5 - port0/line3
            'S6' [0 4];... % S6 - port0/line4
            'S7' [0 0];... % S7 - port0/line0
            };
        DAQ
        nChannels
        
        tID % internal timer
        
        Data % current data
        Datap % previous data
        TOA % Time of access 1*n
        ReadoutTime = 0 % sec to store data before refresh 1*n
        BBoxReadout = false
        BBoxWaitForRealease = false % wait for release instead of press
        
        BTable
        
        isDAQ
        isPTB
    end
    
    properties (Dependent)
        IsValid
        
        Buttons
        Clock
        
        TimeOfLastButtonPress
    end
    
    methods
        
        %% Contructor and destructor
        function obj = MEGSynchClass(varargin)
            fprintf('Initialising MEG Synch...');
            % Create session
            
            buttons = {...
                'LY' 'S3';...
                'RB' 'S4';...
                'RY' 'S5';...
                'RG' 'S6';...
                'RR' 'S7';...
                };
            noDAQ = false;
            for arg = varargin
                if iscell(arg{1}), buttons = arg{1}; end
                if islogical(arg{1}) || isnumeric(arg{1}), noDAQ = logical(arg{1}); end
            end
            obj.BTable = table(buttons(:,1));
            obj.BTable.Properties.VariableNames = {'Buttons'};
            obj.BTable.Properties.RowNames = buttons(:,2);
            
            if ~noDAQ
                warning off daq:Session:onDemandOnlyChannelsAdded
                obj.DAQ = daq.createSession('ni');
                % Add channels for buttons
                for b = obj.BTable.Properties.RowNames'
                    patch = obj.Map{strcmp(obj.Map(:,1),b{1}),2};
                    obj.DAQ.addDigitalChannel('Dev1', sprintf('port%d/line%d',patch(1),patch(2)), 'InputOnly');
                end
                obj.nChannels = numel(obj.DAQ.Channels);
                
                % Add channels for trigger
                for l = 0:7
                    obj.DAQ.addDigitalChannel('Dev1', sprintf('port2/line%d',l), 'OutputOnly');
                end
                obj.isDAQ = 1;
            else
                obj.isDAQ = 0;
                obj.DAQ.isvalid = 1;
                obj.DAQ.Channels = 1:5;
                obj.nChannels = numel(obj.DAQ.Channels);
                fprintf('\n');
                fprintf('WARNING: DAQ card is not in use!\n');
                fprintf('Emulation: ButtonBox is not available           --> ');
                fprintf('You may need to set Keys!\n');
            end
            
            if ~obj.IsValid
                warning('WARNING: MEG Synch is not open!');
            end
            
            obj.isPTB = exist('KbCheck','file') == 2;
            
            obj.Data = zeros(1,obj.nChannels);
            obj.Datap = zeros(1,obj.nChannels);
            obj.ReadoutTime = obj.ReadoutTime * ones(1,obj.nChannels);
            obj.ResetClock;
            fprintf('Done\n');
        end
        
        function delete(obj)
            fprintf('MEG Synch is closing...');
            if obj.isDAQ
                obj.DAQ.release();
                delete(obj.DAQ);
                warning on daq:Session:onDemandOnlyChannelsAdded
            end
            if obj.isPTB, ListenChar(0); end
            fprintf('Done\n');
        end
        
        %% Utils
        function val = get.IsValid(obj)
            val = ~isempty(obj.DAQ) && obj.DAQ.isvalid;
        end
        
        function ResetClock(obj)
            obj.tID = tic;
            obj.TOA = zeros(1,obj.nChannels);
        end
        
        function val = get.Clock(obj)
            val = toc(obj.tID);
        end
        
        function set.Keys(obj,val)
            obj.Keys = val;
            if obj.isPTB, ListenChar(2); end % suppress passing keypresses to MATLAB
        end
        
        %% Trigger
        function SendTrigger(obj,val)
            if obj.isDAQ
                outputSingleScan(obj.DAQ,dec2binvec(val,8));
            end
        end
        
        %% Buttons
        function SetButtonReadoutTime(obj,t)
            obj.ReadoutTime(:) = t;
            obj.BBoxReadout = false;
        end
        
        function SetButtonBoxReadoutTime(obj,t)
            obj.ReadoutTime(:) = t;
            obj.BBoxReadout = true;
        end
        
        function WaitForButtonPress(obj,timeout,ind)
            BBoxQuery = obj.Clock;
            
            if ~exist('timeout','var') || isempty(timeout)
                timeout = Inf;
            end
            
            if ~exist('ind') || isempty(ind)
                ind = 1:5;
            end
            
            % find button label
            buttonlabel = obj.BTable(ind,:).Buttons;
            
            % Reset indicator
            obj.LastButtonPress = [];
            
            % timeout
            wait = timeout < 0; % wait until timeout even in case of response
            timeout = abs(timeout);
            
            while (~obj.Buttons ||... % button pressed
                    wait || ...
                    isempty(intersect(obj.LastButtonPress,buttonlabel))) && ... % correct button pressed
                    (obj.Clock - BBoxQuery) < timeout % timeout
            end
        end
        
        function WaitForButtonRelease(obj,varargin)
            % backup settings
            rot = obj.ReadoutTime;
            bbro = obj.BBoxReadout;
            
            % config for release
            obj.BBoxWaitForRealease = true;
            obj.SetButtonBoxReadoutTime(0);
            
            WaitForButtonPress(obj,varargin{:});
            
            % restore settings
            obj.BBoxWaitForRealease = false;
            obj.ReadoutTime = rot;
            obj.BBoxReadout = bbro;
        end
        
        function val = get.TimeOfLastButtonPress(obj)
            val = max(obj.TOA) * ~isempty(obj.LastButtonPress);
        end
     
        %% Low level access
		function val = get.Buttons(obj)
            val = 0;
            obj.Refresh;
            if obj.BBoxWaitForRealease
                if any(obj.Datap) && all(~(obj.Data.*obj.Datap))
%                     obj.LastButtonPress = find(obj.Datap);
                    obj.LastButtonPress = obj.BTable(find(obj.Datap),:).Buttons;
                    obj.Datap(:) = 0;
                    val = 1;
                end
            else
                if any(obj.Data)
                    %                 obj.LastButtonPress = find(obj.Data);
                    obj.LastButtonPress = obj.BTable(find(obj.Data),:).Buttons;
                    obj.Data(:) = 0;
                    val = 1;
                end
            end
        end
	end
	
	methods (Access = private)	
        function Refresh(obj)
			t = obj.Clock;
			
            % get data
			if obj.isDAQ
                data = inputSingleScan(obj.DAQ); % inverted
            else
                data = zeros(1,obj.nChannels);
                % button press emulation (keyboard) via PTB
                nKeys = numel(obj.Keys);
                if obj.isPTB && nKeys
                    [ ~, ~, keyCode ] = KbCheck;
                    data(1:nKeys) = keyCode(KbName(obj.Keys));
                end
            end
            
            if obj.BBoxReadout, obj.TOA = max(obj.TOA); end
            ind = obj.ReadoutTime < (t-obj.TOA);
            obj.Datap = obj.Data;
            obj.Data(ind) = data(ind);
            obj.TOA(logical(obj.Data)) = t;
        end
        
    end
    
end
