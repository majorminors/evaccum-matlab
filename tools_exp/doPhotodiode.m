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
    
    Screen('FillRect',p.win,colourToUse,[p.photodiodeXshift p.photodiodeYshift p.photodiodeRectwidth p.photodiodeRectHeight]); %show the photoresistor blank rectangle
    
end

return
end