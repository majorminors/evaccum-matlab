function check_headmov(movecompfiles)
%move compfiles is a structure with maxfilter's output file names
%e.g.

% movecompfiles = {'outputfile1.log', ...
%                  'outputfile2.log', ...
%                  'outputfile3.log'};

% uses check_movecomp to read and display maxfilter -movecomp output
% Files with maxfilter outputs, e.g. for different subjects and conditions
% NOTE: the green line ('goodness of fit', gof) indicates how well the hpi
% coils could be detected and should always be 1 (without accurate
% localisation of these coils, the rest would be meaningless).
%Translation (red) rotation and drift indicate whether the subject moved
%their head. These values should be as low as possible, ideally below 0.3.



nr_files = length(movecompfiles);

for ff = 1:nr_files,
    
    fprintf(1, 'Processing %s\n', movecompfiles{ff});
    [mv_fig, linet, linee, lineg, linev, liner, lined] = check_movecomp(movecompfiles{ff}); % read info from log-file
    figure( mv_fig );
    try
    [a,b,c,d] = fileparts( movecompfiles{ff} );
    tittxt = [a(end-17:end) '/' b]; % just one way to keep the figure title short
    % check_movecomp only creates the figure, but without the title, so you can do this separately
    ht = title(tittxt); set(ht, 'interpreter', 'none'); % add title to figure
    catch
    end
end;
