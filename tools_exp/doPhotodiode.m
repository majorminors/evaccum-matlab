function doPhotodiode(p,direction)
%% puts something on the screen for a photoresistor to pick up
% requires:
%   p.usePhotodiode
%   p.win
%   p.photodiodeOnColour
%   p.photodiodeXshift
%   p.photodiodeYshift
%   p.photodiodeRectwidth
%   p.photodiodeRectHeight
%
% does not flip screen


if p.usePhotodiode
    
    if strcmp(direction,'on')
        
        colourToUse = p.photodiodeOnColour;
        
    elseif strcmp(direction,'off')
        
        colourToUse = p.photodiodeOffColour;
        
    end
      
    Screen('DrawDots',p.win,[p.photodiodeXshift p.photodiodeYshift],p.photodiodeDiameter,colourToUse,[],1); % draw a dot for the photodiode (last two values []=default centre, 1 = circular (not square) dots)
    
end

return
end
