% pix = angle2pix(p,ang)
%
% calculates pixel size from visual angles, assuming isotropic (square) pixels
%
% requires:
% p.screen_distance           distance from screen in cm
% p.screen_width              width of screen in cm
% p.resolution                number of pixels of p in horizontal direction - this needs to be calculated after the window is opened for the dots
%% adapted from G.M. Boynton (University of Washington)
%% last edit D. Minors 18 November 2019
%% start function

function pix = angle2pix(p,ang) 
pixSize = p.screen_width/p.resolution(1);   % cm/pix
sz = 2*p.screen_distance*tan(pi*ang/(2*180));  %cm
pix = round(sz/pixSize);   % pix 
return
end