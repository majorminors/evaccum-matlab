% pix = angle2pix(p,ang)
% calculates pixel size from visual angles, assuming isotropic (square) pixels

function pix = angle2pix(p,ang) 
pixSize = p.screen_width/p.resolution(1);   % cm/pix
sz = 2*p.screen_distance*tan(pi*ang/(2*180));  %cm
pix = round(sz/pixSize);   % pix 
return
end