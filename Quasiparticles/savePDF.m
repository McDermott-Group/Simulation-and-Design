function savePDF(h, filename)
%SAVEPDF Saves a figure specified by handle h to a pdf file named filename.

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize',...
    [pos(3), pos(4)])
% print(h, filename, '-dpdf', '-r0', '-bestfit')
% print(h, filename, '-dpdf', '-r0', '-fillpage')
print(h, filename, '-dpdf', '-r0')

end