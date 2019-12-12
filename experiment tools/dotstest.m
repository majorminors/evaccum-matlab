dots.direction = 45; % pulls direction in degrees from 'p.stim_mat' to give to moving_dots
dots.coherence = 0.2; % pulls coherence values from 'p.stim_mat' to give to moving_dots
dots.aperture_size = [10,10]; % [width,height] in degrees of the rectangular aperture dots are displayed in
dots.centre = [0,0]; % [x,y] centre of the dot cloud
dots.num_dots = 100; % number of dots
dots.colour = [255,255,255]; % colour of the dots in [r,g,b]
dots.visual_angle = 0.15; % visual angle of the dots expressed as a decimal - determines size
dots.speed = 5; % speed of dots in degrees per second
dots.lifetime = 5; % number of frames dots live for

p = struct(); % est structure for parameter values
d = struct(); % est structure for trial data
t = struct();

p.screen_width = 35;   % Screen width in cm
p.screen_distance = 50;
p.dots_duration = 5;
p.bg_colour = [0 0 0];

try
    
    % open screen
    [p.win,p.rect] = Screen('OpenWindow',0,p.bg_colour);
    % then need some info based on the screen for later
    p.frame_rate = 1/Screen('GetFlipInterval', p.win); % is Hz - for moving_dots
    p.resolution = p.rect([3,4]); % pull resolution info from p.rect - used to scale cue image and is passed to moving_dots to do the same
    HideCursor;
    WaitSecs(0.5); % warm up
    
    KbQueueCreate(); % initialises queue to collect response information from the list we made (not listening for response yet)
    KbQueueStart(); % starts delivering keypress info to the queue
    [d.dots_onset, t.pressed, t.firstPress] = moving_dots(p,dots); % pull time of first flip for dots, as well as information from KBQueueCheck from moving_dots
    
    Screen('Flip', p.win);
    
    % close screen
    ShowCursor;
    KbQueueRelease(); %KbReleaseWait();
    clear ans; % clear specific indexes and stuff
    Screen('Close',p.win);
    
    fprintf('done running %s\n', mfilename);
    
catch err
    save(save_file);
    ShowCursor;
    KbQueueRelease(); %KbReleaseWait();
    sca; %Screen('Close',p.win);
    rethrow(err);
end