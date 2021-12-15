%%

p.fullscreen_enabled = 1;
p.screen_num = 0; % screen to display experiment on (0 unless multiple screens)
p.bg_colour = [0 0 0]; % needs to be the same as the cue stimuli background colour (unless transparent)
p.text_colour = [255 255 255]; % colour of instructional text
p.window_size = [0 0 1200 800]; % size of window when ~p.fullscreen_enabled
p.text_size = 40; % size of text


p.usePhotodiode = 1;
p.photodiodeOnColour = [255 255 255]; % colour when photo diode is on
p.photodiodeOffColour = [0 0 0]; % colour when photo diode is off
p.photodiodeRectHeight = 70; % pixel height of rectangle to display
p.photodiodeRectWidth = 70; % pixel width of rectangle to display
p.photodiodeXshift = 120; % shift rectangle x pixels off edge of screen
p.photodiodeYshift = 1080-45; % shift rectangle y pixels off edge of screen
% so we thought we'd use a rectangle, but we'll actually do a circle to minimise surface area
% using drawdots
p.photodiodeDiameter = 50; % I think you can do 50 on the MEG computer, but with drawdots this will vary based on graphics hardware

addpath(genpath(fullfile(pwd, 'tools_exp'))); % add tools folder to path (includes moving_dots function which is required for dot motion, as well as an external copy of subfunctions for backwards compatibility with MATLAB)

try

% psychtoolbox setup
AssertOpenGL; % check Psychtoolbox (on OpenGL) and Screen() is working
KbName('UnifyKeyNames'); % makes key mappings compatible (mac/win)
if p.fullscreen_enabled % zero out p.window_size if p.fullscreen_enabled = 1
    p.window_size=[];
end
[p.win,p.rect] = Screen('OpenWindow',p.screen_num,p.bg_colour,p.window_size);
Screen('BlendFunction',p.win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); % allows transparency in .png images
Screen('TextSize', p.win, p.text_size); % set the text size
% then need some info based on the screen for later
p.frame_rate = 1/Screen('GetFlipInterval', p.win); % is Hz - for moving_dots
p.resolution = p.rect([3,4]); % pull resolution info from p.rect - used to scale cue image and is passed to moving_dots to do the same
HideCursor;
WaitSecs(0.5); % warm up

doPhotodiode(p,'on');

DrawFormattedText(p.win, sprintf('\nadjust screen to photodiode\npress y to to continue'), 'center', 'center', p.text_colour);

Screen('Flip', p.win);

while 1 % loop forever until y or n
    ListenChar(2);
    [secs,keyCode] = KbWait; % wait for response
    key_name = KbName(keyCode); % find out name of key that was pressed
    if strcmp(key_name, 'y')
        ListenChar(0);
        clear secs keyCode key_name
        break % break the loop and continue
    elseif strcmp(key_name, 'n')
        ListenChar(0);
        clear secs keyCode key_name
        error('not happy\naborting %s\n', mfilename); % error out
    end
end % end response loop


doPhotodiode(p,'off');

% tell them it's over
DrawFormattedText(p.win,'done!', 'center', 'center', p.text_colour); % tell them it's over!
Screen('Flip', p.win);
WaitSecs(1);
ShowCursor;
Screen('CloseAll');

catch err
    
    ShowCursor;
 
    sca; %Screen('Close',p.win);

    rethrow(err);
    
end