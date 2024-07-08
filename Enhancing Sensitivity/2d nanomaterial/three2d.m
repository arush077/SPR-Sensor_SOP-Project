%2d nanomaterials


no_layers=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];

%gr
sensitivity_gr=[173 173.5 174 174.5 174.5 174.5 175 175.5 176 176.5 176.5 176.5 177 177.5 177.5 178 178.5 179 179] ;
fwhm_gr=[0.24 0.255 0.27 0.285 0.3 0.315 0.33 0.345 0.355 0.37 0.385 0.405 0.42 0.44 0.455 0.475 0.49 0.505 0.525];
fom_gr=[720.8333 680.3922 644.444 612.2807 581.6667 553.9683 530.303 508.6957 495.7746 477.027 458.4416 435.8025 421.4286 403.4091 390.1099 374.7368 364.2857 354.4554 340.9524];
rmin_gr=[0.017364 0.031512 0.048132 0.066437 0.085641 0.10539 0.12555 0.14584 0.16608 0.18612 0.20586 0.22521 0.24404 0.2642 0.28033 0.29775 0.31467 0.33109 0.347];

%bp
sensitivity_bp=[173.5 174 174.5 175 175.5 176 176.5 177 177.5 178 178.5 179 180 180.5 181 181 181.5 182 183];
fwhm_bp=[0.25 0.265 0.285 0.305 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.465 0.485 0.51 0.525 0.55 0.575 0.6 0.62];
fom_bp=[694 656.6038 612.2807 573.7706 548.4375 517.6471 490.2778 465.7891 443.75 423.8095 405.6818 384.9462 371.134 353.9216 344.76 329.0909 315.6522 303.333 295.1613];
rmin_bp=[0.024 0.043 0.064963 0.0883 0.11256 0.13701 0.16135 0.18534 0.2088 0.23162 0.25375 0.27513 0.29575 0.31555 0.3341 0.35295 0.37058 0.38753 0.40382];


%fg
sensitivity_fg=[173 173.5 173.5 174 174 174.5 174.5 175 175 175.5 176 176 176 176.5 177 177 177 177.5 178];
fwhm_fg=[0.215 0.22 0.215 0.215 0.215 0.22 0.22 0.215 0.215 0.22 0.22 0.215 0.215 0.22 0.22 0.215 0.215 0.215 0.22];
fom_fg=[804.65 788.63 806.97 809.30 809.3023 793.181 793.1818 813.9535 813.9535 797.72 800 818.6047 818.6047 802.2727 804.54 823.2558 823.2558 825.58 809.0909];
rmin_fg=[0.000603 0.000802 0.001 0.0006 0.00048 0.00064 0.00099 0.00056 0.0003 0.00048 0.00084 0.000567 0.00035 0.00033 0.000576 0.000673 0.000331 0.000214 0.00031764];










% Define limits
sensitivity_limit = [170, 190];
fom_limit = [160, 860];

% Plotting
figure;
[ax, h1, h2] = plotyy(no_layers, sensitivity_gr, no_layers, fom_gr, 'plot');

% Setting limits
ylim(ax(1), sensitivity_limit);
ylim(ax(2), fom_limit);

% Customizing plot
title('Sensitivity and FOM vs. Number of Layers of Gr');
xlabel('Layers of Gr');
ylabel(ax(1), 'Sensitivity');
ylabel(ax(2), 'FOM');

% Adding legend
legend([h1, h2], {'Sensitivity', 'FOM'}, 'Location', 'northwest');

%-----






% Define limits
sensitivity_limit = [170, 190];
fom_limit = [160, 860];

% Plotting
figure;
[ax, h1, h2] = plotyy(no_layers, sensitivity_bp, no_layers, fom_bp, 'plot');

% Setting limits
ylim(ax(1), sensitivity_limit);
ylim(ax(2), fom_limit);

% Customizing plot
title('Sensitivity and FOM vs. Number of Layers of bp');
xlabel('Layers of bp');
ylabel(ax(1), 'Sensitivity');
ylabel(ax(2), 'FOM');

% Adding legend
legend([h1, h2], {'Sensitivity', 'FOM'}, 'Location', 'northwest');



%-----



% Define limits
sensitivity_limit = [170, 190];
fom_limit = [160, 860];

% Plotting
figure;
[ax, h1, h2] = plotyy(no_layers, sensitivity_fg, no_layers, fom_fg, 'plot');

% Setting limits
ylim(ax(1), sensitivity_limit);
ylim(ax(2), fom_limit);

% Customizing plot
title('Sensitivity and FOM vs. Number of Layers of FG');
xlabel('Layers of FG');
ylabel(ax(1), 'Sensitivity');
ylabel(ax(2), 'FOM');

% Adding legend
legend([h1, h2], {'Sensitivity', 'FOM'}, 'Location', 'northwest');










