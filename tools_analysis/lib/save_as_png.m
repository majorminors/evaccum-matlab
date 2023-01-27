function save_as_png(handle, filename, dpi, width, height)
set(handle, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height] / dpi);
print(handle, '-dpng', ['-r' num2str(dpi)], filename);
end
