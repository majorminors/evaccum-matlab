% Define the number of colors in the colormap
numColors = 181;

% Create a range of values from 0 to 1, with a power of 0.7 for a smoother transition
x = linspace(0, 1, numColors)'.^0.8;

% Define the RGB values for the purple, orange, and middle colors
purpleColour = [0.6, 0.4, 0.8]; % RGB values for purple
orangeColour = [0.9, 0.6, 0.2]; % RGB values for orange
middleColour = [0.8, 0.8, 0.8]; % RGB values for middle color

% Interpolate the RGB values between purple and middle, and between middle and orange
colours1 = middleColour + (purpleColour - middleColour) .* x;
colours2 = middleColour + (orangeColour - middleColour) .* x;
colours = [flip(colours2); colours1];

% Set the custom colormap
colormap(colours);

% Display a colorbar for visualization
colorbar;

save('bayes_colourmap.mat', 'colours');