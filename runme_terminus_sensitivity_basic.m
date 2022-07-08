% Setup
region = 'SAtoES';
start_year = 1985;

%% Mesh sizing
triangleresolution = 1000;

%% Mesh
md = model()
md = triangle(md,['./Exp/' region '.exp'],triangleresolution);

savemodel(['./Models/' region '_mesh.mat'], 'md', md)

%% Parameterize
md = parameterize(md,'Par/Greenland.par');
md = setflowequation(md,'SSA','all');

savemodel(['./Models/' region '_param.mat'], 'md', md)

%% Inversion
md = loadmodel(['./Models.' region '_Param.mat']);
% Control general
md.inversion=m1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.verbose=verbose('solution',false,'control',true);

% Cost functions
md.inversion.cost_functions=[101 103 501]; %Abs, Log, reg
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
md.inversion.cost_functions_coefficients(:,1)=2000;
md.inversion.cost_functions_coefficients(:,2)=40;
md.inversion.cost_functions_coefficients(:,3)=1.6e-06;

%Where vel==0, set coefficients to 0 (i.e., don't try to match this in model
disp(['Removing vel==0 obs from inversion']);
pos = find(md.inversion.vel_obs == 0);
md.inversion.cost_functions_coefficients(pos,1) = 0;
md.inversion.cost_functions_coefficients(pos,2) = 0;

%Controls
md.inversion.control_parameters={'FrictionCoefficient'};
md.inversion.maxsteps=50;
md.inversion.maxiter =50;
md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);
md.inversion.control_scaling_factors=1;

%Set basal friction coefficient initial guess to something low at front
filename = ['Exp/' region '_coeffront.exp'];
  if exist(filename,'file'),
  disp(['Correcting basal friction coefficient initial guess for front inconsistencies']);
  flags = ContourToNodes(md.mesh.x,md.mesh.y,filename,2);
  %flags = md.inversion.vel_obs == 0;
  pos1 = find(flags); pos2 = find(~flags);
  md.friction.coefficient(pos1,:)  = 10;

  md.inversion.max_parameters(pos1) = md.friction.coefficient(pos1,:);
end

%Additional parameters
md.stressbalance.restol=0.01;
md.stressbalance.reltol=0.1;
md.stressbalance.abstol=NaN;
%md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

% Go solve
md.verbose.solution=1;
md.cluster = cluster;
md.settings.waitonlock = waitonlock;
md = solve(md,'Stressbalance');

% Save
savemodel(org,md);


%% Relaxation
%Put results of inversion back into the model for forward runs
md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

% Special post-processing of inverted friction coefficient
filename = ['Exp/' region '_coeffront_after_inversion.exp'];
if exist(filename, 'file')
  pos = find(ContourToNodes(md.mesh.x, md.mesh.y, filename, 1));
  md.friction.coefficient(pos) = 10;
end

md.initialization.pressure = zeros(md.mesh.numberofvertices,1);
md.initialization.temperature = 250*ones(md.mesh.numberofvertices,1);

% Set parameters
md.inversion.iscontrol=0;
md.timestepping.start_time = start_year;
md.timestepping.time_step  = .02;
md.timestepping.final_time = start_year + relaxation_years;
md.settings.output_frequency = (1/md.timestepping.time_step)/5; % 5/yr

% We set the transient parameters
md.transient.ismovingfront=0;
md.transient.isthermal=0;
md.transient.isstressbalance=1;
md.transient.ismasstransport=1;
md.transient.isgroundingline=1;
md.groundingline.migration = 'SubelementMigration';

% We set the calving model (no moving front ... spclevelset is actually ignored)
md.levelset.spclevelset = md.levelset.spclevelset;
md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

% Set the requested outputs
md.stressbalance.requested_outputs={'default'};
md.transient.requested_outputs={'default'};

% Go solve
md.verbose.solution=1;
md.cluster = cluster;
md.settings.waitonlock = waitonlock;
md = solve(md,'transient');