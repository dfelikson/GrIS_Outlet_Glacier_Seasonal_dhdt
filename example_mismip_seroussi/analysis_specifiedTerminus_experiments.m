
mdSeasonal = loadmodel('Models_2km_viscous/MISMIP_2km_viscous_Retreat_3yrSeasonal.mat');
mdAverage  = loadmodel('Models_2km_viscous/MISMIP_2km_viscous_Retreat_3yrAvg.mat');

% Figure 1 - seasonal terminus
plot_bed_and_termini(mdSeasonal, 'index', [1 2 3 4 5 6], 'legend', true);
set(gcf, 'Position', [145 237 1198 440]);
export_fig('terminus_seasonal.png', '-r300')

% Figure 2 - average terminus
plot_bed_and_termini(mdAverage, 'legend', true);
set(gcf, 'Position', [145 237 1198 440]);
export_fig('terminus_average.png', '-r300')

% Figure 3 - mass change (average versus seasona)
tSeasonal = [mdSeasonal.results.TransientSolution(:).time];
vSeasonal = [mdSeasonal.results.TransientSolution(:).IceVolumeAboveFloatation] * md.materials.rho_ice * 1e-12; % Gt;
tAverage  = [mdAverage.results.TransientSolution(:).time];
vAverage  = [mdAverage.results.TransientSolution(:).IceVolumeAboveFloatation] * md.materials.rho_ice * 1e-12; % Gt;
figure
hold on
plot(tSeasonal-tSeasonal(1), vSeasonal-vSeasonal(1))
plot(tAverage-tAverage(1), vAverage-vAverage(1))
grid on
title('mass change')
xlabel('time elapsed (years)')
ylabel('mass change (Gt)')
legend('seasonal', 'average')
set(gcf,'color','w');
export_fig('mass_change_seasonal-vs-average.png', '-r300')

