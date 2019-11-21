function [dots_onset_time, pressed, firstPress] = moving_dots(p,dots,MEG,exp_trial)
% [dots_onset_time, pressed, firstPress] = moving_dots(p,dots,MEG,exp_trial)
%
% creates a cloud of moving dots, then creates a fixation by putting a small 
% black square in a bigger white square, then flips the screen
%
% will display for time provided in 'p' structure (see below) otherwise
% return immediately upon first valid keypress
% 
% multiple dot clouds can be displayed by making 'dots' into an array 'dots(n)' (e.g. dots(1), dots(2))
% 
% order that dots are drawn in uses randperm so clouds won't occlude each other
%
% requires KBQueueStart() called prior to function (so may want to
% call KBQueueFlush() prior to function), as well as information from 'dots'
% structure and 'p' (i.e. experiment parameter) structure (see below)
%
% requires MEG, which is a class of functions for the MEG interface with the National Instruments PCI 6503 card (MRC CBU)
%
% returns:
%    dots_onset_time           time of first flip
%    pressed                   assuming 'KBQueueStart()' has been called prior
%                              to moving_dots, will return whether a key has been
%                              pressed (1) or not (2)
%
%    firstPress                assuming 'KBQueueStart()' has been called
%                              prior to moving_dots, will return an array 
%                              with time and key info of first key pressed
%
% requires from 'dots':
%    dots.num_dots             number of dots in the cloud
%    dots.speed                speed of dots in degrees per second
%    dots.direction            angle in degrees (0 is up)
%    dots.lifetime             number of frames dots live for
%    dots.aperture_size        [width,height] in degrees of the rectangular aperture dots are displayed in
%    dots.centre               [x,y] centre of the dot cloud
%    dots.colour               colour of the dots in [r,g,b]
%    dots.visual_angle         visual angle of the dots expressed as a decimal - determines size
%    dots.coherence            percentage of coherent dots from 0-1 (1=100%)
%
%
% requires from 'p':
%    p.MEG_enabled            a boolean (0 or 1) telling the script whether we should record responses with KbQueue or with MEG functions
%    p.stim_mat               needs this to determine the unique MEG trigger to send for each trial condition
%    p.fix_trial_time         if 0 then trial will end on keypress, if 1 will go for duration of p.dots_duration
%    p.dots_duration          seconds for the dot cloud to be displayed
%    p.screen_width           width of the screen in cm
%    p.screen_distance        distance from the screen in cm
%    p.win                    window pointer for the dots
%    p.frame_rate             frame rate for the window
%    p.resolution             pixel resolution for the window
%
% can also specify with 'p':
%    p.fixation               moving_dots will build a default fixation but these parameters can be changed externally
% 
%    if otherwise unspecified, it will create 'p.fixation' structure with the following default values
%       p.fixation.size       size in degrees for fixation square (0.5)
%       p.fixation.mask       size in degrees for mask around fixation (1)
%       p.fixation.colour     array that specifies colour for outer and inner elements of fixation square (white/black - {[255,255,255],[0,0,0]})
%    will also use:
%       p.bg_colour           background colour - will use for 'p.fixation.mask' around the fixation

%% subclasses required
%% MEGSynchClass
% a class of functions for the MEG interface with the National Instruments PCI 6503 card (MRC CBU)
%% subfunctions (at bottom of moving_dots):
%% pix = angle2pix(p,ang)
% calculates pixel size from visual angles, assuming isotropic (square) pixels

% requires:
% p.screen_distance           distance from screen in cm
% p.screen_width              width of screen in cm
% p.resolution                number of pixels of p in horizontal direction - this needs to be calculated after the window is opened for the dots

%% frames=secs2frames(p,secs)
% converts time in seconds to frames using secs

% requires:
% p.frame_rate                this needs to be calculated after the window is opened for the dots

%% adapted from G.M. Boynton (University of Washington) 
%% last edit D. Minors 18 November 2019
%% start function

% calculate total number of dots across clouds
num_dots = sum([dots.num_dots]);

% zero out the colour and size vectors
colours = zeros(3,num_dots);
sizes = zeros(1,num_dots);

% generate random order so dot clouds won't occlude each other when drawn
order = randperm(num_dots);


%% define dot positions and some other initial parameters

count = 1;

% preallocate arrays  for the loop for the left, right top and bottom of each aperture (in degrees)
l = zeros(1,length(dots));
r = zeros(1,length(dots));
b = zeros(1,length(dots));
t = zeros(1,length(dots));

for i=1:length(dots) % loop through the clouds

    % calculate the left, right top and bottom of each aperture (in degrees)
    l(i) = dots(i).centre(1)-dots(i).aperture_size(1)/2;
    r(i) = dots(i).centre(1)+dots(i).aperture_size(1)/2;
    b(i) = dots(i).centre(2)-dots(i).aperture_size(2)/2;
    t(i) = dots(i).centre(2)+dots(i).aperture_size(2)/2;

    % generate random starting positions
    dots(i).x = (rand(1,dots(i).num_dots)-.5)*dots(i).aperture_size(1) + dots(i).centre(1);
    dots(i).y = (rand(1,dots(i).num_dots)-.5)*dots(i).aperture_size(2) + dots(i).centre(2);

    % create a direction vector for a given coherence level
    direction = rand(1,dots(i).num_dots)*360;
    nCoherent = ceil(dots(i).coherence*dots(i).num_dots);  % start w/ all random directions
    direction(1:nCoherent) = dots(i).direction;  % set the 'coherent' directions

    % calculate dx and dy vectors in real-world coordinates
    dots(i).dx = dots(i).speed*sin(direction*pi/180)/p.frame_rate;
    dots(i).dy = -dots(i).speed*cos(direction*pi/180)/p.frame_rate;
    dots(i).life =    ceil(rand(1,dots(i).num_dots)*dots(i).lifetime);
    
    dots(i).size = angle2pix(p,dots(i).visual_angle); % convert dot size from visual angle to pixels

    % fill in the 'colours' and 'sizes' vectors for this cloud
    id = count:(count+dots(i).num_dots-1);  % index into the num_dots length vector for this cloud
    colours(:,order(id)) = repmat(dots(1).colour(:),1,dots(i).num_dots);
    sizes(order(id)) = repmat(dots(i).size,1,dots(i).num_dots);
    count = count+dots(i).num_dots;
end

% zero out the screen position vectors and the 'good_dots' vector
pixpos.x = zeros(1,num_dots);
pixpos.y = zeros(1,num_dots);
good_dots = false(size(num_dots));

% calculate total number of temporal frames
total_frames = secs2frames(p,p.dots_duration);

%% loop through frames

% if p.MEG_enabled then set up MEG for response recording
if p.MEG_enabled == 1
    MEG.SendTrigger(p.stim_mat(exp_trial,p.MEGtriggers.onsets)); % send a trigger for trial onset
    MEG.ResetClock; % reset the timer
    button_pressed = 0; % a counter to make sure we catch the first time a button was pressed
    pause(0.005); % quick pause before we reset triggers
%     MEG.WaitForButtonPress(p.dots_duration); % listen for button press
    MEG.SendTrigger(0); % reset triggers
end

for frame_num = 1:total_frames
    count = 1;
    
    for i = 1:length(dots)  % loop through the clouds (if more than one)
        

        % update the dot position's real-world coordinates
        dots(i).x = dots(i).x + dots(i).dx;
        dots(i).y = dots(i).y + dots(i).dy;

        % move the dots that are outside the aperture back one aperture width.
        dots(i).x(dots(i).x<l(i)) = dots(i).x(dots(i).x<l(i)) + dots(i).aperture_size(1);
        dots(i).x(dots(i).x>r(i)) = dots(i).x(dots(i).x>r(i)) - dots(i).aperture_size(1);
        dots(i).y(dots(i).y<b(i)) = dots(i).y(dots(i).y<b(i)) + dots(i).aperture_size(2);
        dots(i).y(dots(i).y>t(i)) = dots(i).y(dots(i).y>t(i)) - dots(i).aperture_size(2);

        % increment the 'life' of each dot
        dots(i).life = dots(i).life+1;

        % find the 'dead' dots
        dead_dots = mod(dots(i).life,dots(i).lifetime)==0;

        % replace the positions of the dead dots to random locations
        dots(i).x(dead_dots) = (rand(1,sum(dead_dots))-.5)*dots(i).aperture_size(1) + dots(i).centre(1);
        dots(i).y(dead_dots) = (rand(1,sum(dead_dots))-.5)*dots(i).aperture_size(2) + dots(i).centre(2);

        % calculate the index for this cloud's dots into the whole list of
        % dots.  Using the vector 'order' means that, for example, the first
        % cloud is represented not in the first n values, but rather is
        % distributed throughout the whole list.
        id = order(count:(count+dots(i).num_dots-1));
        
        % calculate the pixel positions for this cloud using visual angle from the real-world
        % coordinates of the dot posns, plus the offset for the centre of the screen
        pixpos.x(id) = angle2pix(p,dots(i).x)+ p.resolution(1)/2;
        pixpos.y(id) = angle2pix(p,dots(i).y)+ p.resolution(2)/2;

        % determine which of the dots in this cloud are outside this cloud's
        % elliptical aperture
        good_dots(id) = (dots(i).x-dots(i).centre(1)).^2/(dots(i).aperture_size(1)/2)^2 + ...
            (dots(i).y-dots(i).centre(2)).^2/(dots(i).aperture_size(2)/2)^2 < 1;
  
        count = count+dots(i).num_dots;
    end
    
    % draw all clouds at once
    Screen('DrawDots',p.win,[pixpos.x(good_dots);pixpos.y(good_dots)], sizes(good_dots), colours(:,good_dots),[0,0],1);
    
%% draw a fixation for the dot motion
    
    % establish default values - build fixation parameters if not already set
    if ~isfield(p,'fixation')
        p.fixation = [];
    end
    % build size of fixation
    if ~isfield(p.fixation,'size')
        p.fixation.size = 0.5; % degrees
    end
    % build mask around fixation (little bit of space around the fixation
    %   so dots don't run into it
    if ~isfield(p.fixation,'mask')
        p.fixation.mask = 0.7;  % degrees
    end
    % build colour of fixation
    if ~isfield(p.fixation,'colour')
        p.fixation.colour = {[255,255,255],[0,0,0]};
    end
    
    % first establish 3 box sizes (one for each element in the fixation) using screen coords
    sz(1) = angle2pix(p,p.fixation.size/2); % size for outer circle
    sz(2) = angle2pix(p,p.fixation.size/4); % size for inner circle
    sz(3) = angle2pix(p,p.fixation.mask/2); % size for mask
    
    % establish 3 rectangles (one for each element in the fixation) based on the boxes using screen coords [l,t,r,b]
    centre = p.resolution/2;
    for i=1:3
        rect{i}= [-sz(i)+centre(1),-sz(i)+centre(2),sz(i)+centre(1),sz(i)+centre(2)];
    end
    
    % draw the mask
    Screen('FillOval', p.win, p.bg_colour,rect{3});
    % draw the outer circle - defaults to white (this is the circle they'll 'see')
    Screen('FillOval', p.win, p.fixation.colour{1},rect{1});
    % draw the inner circle - defaults to black
    Screen('FillOval', p.win, p.fixation.colour{2},rect{2});
    
    
    %% flip the screen and check for keypress
      
    if frame_num == 1
        dots_onset_time = Screen('Flip',p.win);
    else
        Screen('Flip',p.win);
    end
    
    % keypress stuff within the frame loop
    if p.MEG_enabled == 0
        if p.fix_trial_time == 0
            [pressed,firstPress] = KbQueueCheck(); % check for keypress in the KbQueue every frame
            if any(pressed)
                break % break the frame-loop which will end the function
            end
        end
    elseif p.MEG_enabled == 1
        MEG.WaitForButtonPress(0);
        if ~isempty(MEG.LastButtonPress) && ~button_pressed % check for a keypress in the MEG key wait function every frame, if a key hasn't been pressed yet
            fprintf('response!\n')
            button_pressed = 1; % record that a key has been pressed this trial
            MEG.SendTrigger(p.stim_mat(exp_trial,p.MEGtriggers.responses)); % send a trigger
            firstPress{1} = MEG.LastButtonPress; % record the key pressed
            firstPress{2} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
            pause(0.005); % quick pause before resetting
            MEG.SendTrigger(0); % reset the triggers
            if p.fix_trial_time == 0
                break % break the frame-loop which will end the function
            end
%             fpidx = 0;
%         elseif ~strcmp(MEG.LastButtonPress,firstPress{1}) && button_pressed > 0 % if it's not the first time a key has been pressed
%             button_pressed = 2; % record that another key has been pressed this trial
%             fpidx = fpidx+1;
%             firstPress{3+fpidx} = MEG.LastButtonPress; % record the key pressed
%             firstPress{4+fpidx} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
        end
    end

end

% keypress stuff outside the frame loop
if p.MEG_enabled == 0
    if p.fix_trial_time == 1
        [pressed,firstPress] = KbQueueCheck(); % check for a keypress outside the frame loop before the function ends
    end
elseif p.MEG_enabled == 1 
    pressed = 0; % just put something in here because its used in KbQueue and moving_dots expects an output
    if button_pressed == 0
        firstPress{1} = cellstr('RR');
        firstPress{2} = 0;
    end
    %     [~,checkquit] = KbQueueCheck(); % do a check to see if we quit, and if so error out
    %     if strcmp(KbName(checkquit),p.quitkey)
    %         fclose('all');
    %         error('%s quit by user (p.quitkey pressed)\n', mfilename);
end


return


%% subfunctions

%% convert visual angles in degrees to pixels
function pix = angle2pix(p,ang) 
pixSize = p.screen_width/p.resolution(1);   % cm/pix
sz = 2*p.screen_distance*tan(pi*ang/(2*180));  %cm
pix = round(sz/pixSize);   % pix 
return
% adapted from G.M. Boynton (University of Washington) 

%% convert time in seconds to frames
function frames = secs2frames(p,secs)
frames = round(secs*p.frame_rate);
% adapted from G.M. Boynton (University of Washington) 
