function [mv_fig, linet, linee, lineg, linev, liner, lined] = check_movecomp(mvcomp_file)
% read and plot maxfilter output from log-file
% [mv_fig, linet, linee, lineg, linev, liner, lined] = check_movecomp(mvcomp_file);
colours = {'b', 'g', 'r', 'm', 'c', 'k', 'y'};
fprintf(1, 'Input file: %s.\n', mvcomp_file);
fid_in = fopen(mvcomp_file);    % open movecomp log-file
now_line = fgetl(fid_in);
cc = 0;
fprintf(1, 'Reading line ');
while feof(fid_in) == 0,
    if ~isempty(deblank(now_line)),
        cc = cc + 1;
        alllines{cc} = now_line;
        if ~mod(cc,100), fprintf(1, '%d ', cc); end;
        if ~mod(cc,1000), fprintf(1, '\n', cc); end;
    end;
    now_line = fgetl(fid_in);
end;    % while...
nr_lines = cc;
fprintf(1, '\n%d lines read\n', nr_lines);
HPIline = 'Hpi fit OK, movements [mm/s]';
lenH = length(HPIline);
Tline = '#t =';
lenT = length(Tline);
c_fit = 0;
t_fit = 0;
for cc = 1:nr_lines,
    thisline = alllines{cc};
    n = length(thisline);

    if n>=lenH && strcmp(thisline(1:lenH),HPIline),
        c_fit = c_fit+1;
        [tmp, restof] = strtok(thisline, '=');
        [tmp, restof] = strtok(restof(2:end), '/');
        cnt=0;
        while ~isempty(deblank(restof)),
            if length(find(findstr(tmp, '-')))==0,
                cnt=cnt+1;
                xyz(cnt,c_fit) = str2num(tmp);
            end;
            [tmp, restof] = strtok(restof(2:end), '/');
        end;
    end;

    if n>=lenT && strcmp(thisline(1:lenT),Tline),
        t_fit = t_fit+1;

        [tmp, restof] = strtok(thisline, 't');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof, ',');
        linet(t_fit) = str2num( tmp(2:end) );

        [tmp, restof] = strtok(thisline, 'e');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof);
        [tmp, restof] = strtok(restof);
        linee(t_fit) = str2num( tmp(2:end) );

        [tmp, restof] = strtok(thisline, 'g');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof);
        [tmp, restof] = strtok(restof, ',');
        lineg(t_fit) = str2num( tmp(2:end) );

        [tmp, restof] = strtok(thisline, 'v');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof);
        [tmp, restof] = strtok(restof);
        linev(t_fit) = str2num( tmp(2:end) );

        [tmp, restof] = strtok(thisline, 'r');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof);
        [tmp, restof] = strtok(restof);
        liner(t_fit) = str2num( tmp(2:end) );

        [tmp, restof] = strtok(thisline, 'd');
        [tmp, restof] = strtok(restof, '=');
        [tmp, restof] = strtok(restof);
        [tmp, restof] = strtok(restof);
        lined(t_fit) = str2num( tmp(2:end) );

    end;

end;
mv_fig(1) = figure;
subplot(2,1,1);
plot(linet, linee, colours{1}, linet, lineg, colours{2}, linet, linev, colours{3}, linet, liner, colours{4}, linet, lined, colours{5});
grid on
%legend( {'#e (cm)', '#g', '#v (cm/s)', '#r (rad/s)', '#d (cm)'} );
legend( {'fitting error (cm)', 'gof (0-1)', 'translation (cm/s)', 'rotation (rad/s)', 'drift (cm)'} );
if exist('xyz', 'var') && ~isempty(xyz),
    if size(xyz,2)~=length(linet),
        linet = 1:size(xyz,2);
        fprintf(1, '\n\nHPI problems!!! %s\n\n', mvcomp_file);
    end;
    subplot(2,1,2);
    plot(linet, xyz(1,:), colours{1}, linet, xyz(2,:), colours{2}, linet, xyz(3,:), colours{3});
    legend( {'x (mm/s)', 'y (mm/s)', 'z (mm/s)'} );
    grid on
else,
    fprintf(1, '\n\nHPI problems!!! %s\n\n', mvcomp_file);
end;