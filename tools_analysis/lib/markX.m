function markX(xValuesToMark)
% this might not work, but needs an array of xvalues
% requires a call to hold on prior to calling
for pidx = 1:numel(xValuesToMark)
    plot([xValuesToMark(pidx) xValuesToMark(pidx)],[ylim],'r');
end; clear pidx
end