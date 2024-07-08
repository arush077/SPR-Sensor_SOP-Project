thickness=[35 40 45 50 55 60];
sensitivity=[178.5 178.5 178.5 178 178 178.5]
fom=[346.601 451.8987 549.2308 635.7143 684.6154 714]

fwhm = [0.515 , 0.395 , 0.325 , 0.28 , 0.26 , 0.25]
rmin = [0.012608 , 0.0075425 , 0.07887 , 0.20868 , 0.35842 ,0.50423]












%30 92    274.6269
%30 0.335 0.00047719


%35 1385.5 2690.2913

% Plotting using plotyy
figure;
[ax, h1, h2] = plotyy(thickness, fwhm, thickness, rmin, 'plot');

% Customizing plot
title('FWHM and Rmin vs. Thickness of Ag');
xlabel('thickness of Ag');
ylabel(ax(1), 'FWHM');
ylabel(ax(2), 'Rmin');

% Adding legend
%legend([h1, h2], {'Rmin', 'FWHM'}, 'Location', 'northwest');