
disp('loading models')
md = struct('name', '', 'path', '', 'model', '');

i = 1;
mds(i).name = 'average';
mds(i).path = 'Models_AIS/MISMIP_AIS_Retreat_specified1.mat';
mds(i).model = loadmodel(mds(i).path);

i = i + 1;
mds(i).name = 'average-seasonal';
mds(i).path = 'Models_AIS/MISMIP_AIS_Retreat_specified2.mat';
mds(i).model = loadmodel(mds(i).path);

%i = i + 1;
%mds(i).name = 'seasonal';
%mds(i).path = 'Models_AIS/MISMIP_AIS_Retreat_specified4.mat';
%mds(i).model = loadmodel(mds(2).path);

i = i + 1;
mds(i).name = 'seasonal';
mds(i).path = 'Models_AIS/MISMIP_AIS_Retreat_specified4.mat';
mds(i).model = loadmodel(mds(i).path);


disp('plotting')
for i = 1:numel(md)
   % Figure - seasonal terminus
   if size(mds(i).model.levelset.spclevelset,2) == 2
      plot_bed_and_termini(mds(i).model, 'legend', true);
   else
      plot_bed_and_termini(mds(i).model, 'index', [1 2 3 4 5 6], 'legend', true);
   end
   set(gcf, 'Position', [145 237 1198 440]);
   export_fig(['terminus_' mds(i).name '.png'], '-r300')
end

% Figure - mass change (average versus seasona)
legend_str = {};
figure
hold on
for i = 1:numel(mds)
   t = [mds(i).model.results.TransientSolution(:).time];
   v = [mds(i).model.results.TransientSolution(:).IceVolumeAboveFloatation] * mds(i).model.materials.rho_ice * 1e-12; % Gt;
   plot(t-t(1), v-v(1))
   legend_str{end+1} = mds(i).name;
end
grid on
title('mass change')
xlabel('time elapsed (years)')
ylabel('mass change (Gt)')
legend(legend_str)
set(gcf,'color','w');
export_fig('mass_change_specifiedTerminus_experiments.png', '-r300')

