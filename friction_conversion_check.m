clear all; close all; clc;

disp('loading')
md_orig = loadmodel('Models/SAtoES_friction_coefficient_MEaSUREsVel/SAtoES_Param.mat');

N = effectivepressure(md_orig);
vb = basal_vel_from_surface_vel_obs(md_conv);

%[~, ~, b1_init] = basalstress(md_orig);
s = averaging(md_orig,1./md_orig.friction.p,0);
r = averaging(md_orig,md_orig.friction.q./md_orig.friction.p,0);
b_orig = N.^r .* (md_orig.friction.coefficient).^2 .* vb.^(s-1) .* vb;%  b( ub==0 & (s-1)<0) = 0;

% Friction conversion
disp('friction conversion')
md_conv = md_orig;
m = 2;
md_conv = friction_coefficient_conversion(md_orig, 'budd', 'budd', 'p', 1/m, 'q', 1/m);

%[~, ~, b2_init] = basalstress(md_conv);
s = averaging(md_conv,1./md_conv.friction.p,0);
r = averaging(md_conv,md_conv.friction.q./md_conv.friction.p,0);
b_conv = N.^r .* (md_conv.friction.coefficient).^2 .* vb.^(s-1) .* vb;%  b( ub==0 & (s-1)<0) = 0;

return

%% ------------------------------------------------------------------------------------------- %%
% Stress balance
md_orig = loadmodel('Models/SAtoES_friction_coefficient_MEaSUREsVel/SAtoES_Inversion.mat');
md_orig.friction.coefficient = md_orig.results.StressbalanceSolution.FrictionCoefficient;
md_orig.inversion.iscontrol = 0;

disp('stress balance: md_orig')
md_orig = solve_sb(md_orig);
md_orig.initialization.vx  = md_orig.results.StressbalanceSolution.Vx;
md_orig.initialization.vy  = md_orig.results.StressbalanceSolution.Vy;
md_orig.initialization.vel = md_orig.results.StressbalanceSolution.Vel;

%[~, ~, b1_sb] = basalstress(md_orig);
s = averaging(md_orig,1./md_orig.friction.p,0);
r = averaging(md_orig,md_orig.friction.q./md_orig.friction.p,0);
b_orig_sb = N.^r .* (md_orig.friction.coefficient).^2 .* md_orig.initialization.vel.^(s-1) .* md_orig.initialization.vel;%  b( ub==0 & (s-1)<0) = 0;

% Stress balance
disp('stress balance: md_conv')
md_conv
md_conv = solve_sb(md_conv);
md_conv.initialization.vx  = md_conv.results.StressbalanceSolution.Vx;
md_conv.initialization.vy  = md_conv.results.StressbalanceSolution.Vy;
md_conv.initialization.vel = md_conv.results.StressbalanceSolution.Vel;

%[~, ~, b2_sb] = basalstress(md_conv);
s = averaging(md_conv,1./md_conv.friction.p,0);
r = averaging(md_conv,md_conv.friction.q./md_conv.friction.p,0);
b_conv_sb = N.^r .* (md_conv.friction.coefficient).^2 .* md_conv.initialization.vel.^(s-1) .* md_conv.initialization.vel;%  b( ub==0 & (s-1)<0) = 0;

% Calculate velocity diff
disp('velocity diff')
ub1 = md_orig.results.StressbalanceSolution.Vel;
ub2 = md_conv.results.StressbalanceSolution.Vel;

ub_diff_rel = (ub2-ub1) ./ ub1;
fprintf('\n');
fprintf('Mean absolute relative difference in velocity = %10.8f m/yr\n', mean(abs(ub_diff_rel), 'omitnan'));
fprintf('\n');

