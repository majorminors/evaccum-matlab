function points=simulatecurve

[~,hline] = psychcurve;
points = [hline.XData' hline.YData'];
close(gcf);
