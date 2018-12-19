
steps = [4];
clusterName = ''; % empty for localhost
clusterName = 'oibserve';
%clusterName = 'discover';

%% Setup
region = 'WestGrIS';
start_year = 1985;

% To launch a job remotely and not wait, set the following:
%  md.cluster.interactive = 0;
%  md.settings.waitonlock = nan; or = 0;
switch clusterName %%{{{
   case {'','gs15serac'}
      cluster = generic('name', oshostname(), 'np', 2);
      waitonlock = 10; %nan;

   case 'oibserve'
      cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 12, 'interactive', 0, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
      cluster.interactive = 0; %1;
      waitonlock = nan; %10;

   case 'discover'
      cluster=discover;
      cluster.name='discover.nccs.nasa.gov';
      cluster.login='dfelikso';
      cluster.numnodes=1;
      cluster.cpuspernode=16;
      cluster.time=1*60;
      cluster.interactive=0;
      cluster.processor='sand';
      cluster.queue='allnccs';
      cluster.codepath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/bin';
      cluster.executionpath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/execution';
      cluster.email='denis.felikson@nasa.gov';
      waitonlock = nan;

end
%%}}}

%% Model parameters
org = organizer('repository', ['./Models'], 'prefix', [region '_'], 'steps', steps);

%% Processing {{{
if ~exist(['Exp/' region '.exp'],'file')
   disp(['Domain outline needed: Exp/' region '.exp.'])

   s = input('Is there a shapefile that can be used for the model domain (y/n)?','s');
   if strcmpi(s(1),'y')
      [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/GrIS_Calibrated_SLR/*.shp');
      shp2exp(fullfile(path,file),['Exp/' region '.exp']);
   else
      disp('Exiting ... rerun runme')
      return
   end
end

if perform(org,'Mesh'),% {{{ STEP 1
   % Mesh sizing
   triangleresolution = 1000;
   hmin = 300;
   hmax = 10000;

   fprintf('Creating mesh with the following parameters:\n   %20s = %5d\n   %20s = %5d\n   %20s = %5d\n', 'triangleresolution', triangleresolution, 'hmin', hmin, 'hmax', hmax);
   s = input('-> Continue (y/n)?','s');
   if ~strcmpi(s(1),'y')
      return
   end

   md = triangle(model,['./Exp/' region '.exp'],triangleresolution);

   % Read velocity (for mesh refinement)
   disp('Reading Joughin composite velocities');
   [velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
   vel  = sqrt(velx.^2+vely.^2);

   %Adapt mesh
	disp('Optimizing mesh');

   %Refine mesh beyond terminus (because BAMG refinement above uses velocities that -> 0 where there's still ice in 1985)
   fprintf(['\n\033[' '103;30' 'm   Can I find an automatic way of refining the mesh beyond the terminus of each glacier? \033[0m \n\n']);
   %filename = ['Exp/' glacier '_refineFront.exp'];
   %if ~exist(filename,'file'),
   %   s = input('Is there a shapefile that can be used to refine the front (y/n)?','s');
   %   if strcmpi(s(1),'y')
   %      [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/ModeledInlandThinning/Data/Glacier model domains/*.shp');
   %      shp2exp(fullfile(path,file),['Exp/' glacier '_refineFront.exp']);
   %   else
   %      filename = ['Exp/' glacier '_Front0.5.exp'];
   %   end
   %end

   md = bamg(md,'hmin',hmin,'hmax',hmax,'field',vel,'err',2);

   [md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg = 3413;

	savemodel(org,md);

end %}}}
if perform(org,'Param'),% {{{ STEP 2

	md=loadmodel(org,'Mesh');
	md=parameterize(md,'Par/Greenland.par');
	
   md=setflowequation(md,'SSA','all');
   
   % Weaken shear margins
   filename = ['Exp/' region '_shearmargins.exp'];
   if exist(filename, 'file')
      disp(['Weakening shear margins by 60% using: ' filename]);
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
      md.materials.rheology_B(pos) = 0.60 .* md.materials.rheology_B(pos);
   end

	savemodel(org,md);
end %}}}
if perform(org,'Inversion'),% {{{ STEP 3

	md=loadmodel(org,'Param');

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103 501]; %Abs, Log, reg
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=40;
   switch region
      case 'KAK'
         md.inversion.cost_functions_coefficients(:,3)=.2*50^-3;
      case 'KLG'
         md.inversion.cost_functions_coefficients(:,3)=.2*50^-2;
      otherwise
         md.inversion.cost_functions_coefficients(:,3)=.2*50^-3;
   end

	% %Remove obs where the front from the velocities are upstream of our current front
	% filename = ['Exp/' glacier '_velfront.exp'];
	% if exist(filename,'file'),
	% 	disp(['Correcting cost functions for front inconsistencies']);
	% 	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
	% 	md.friction.coefficient(pos)=min(md.friction.coefficient(pos),100);
	% end
   
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
   %if ~exist(filename,'file'),
   %   plotmodel(md,'data',md.friction.coefficient,'mask',md.mask.ice_levelset<0)
   %   exptool(filename)
   %end
   if exist(filename,'file'),
      disp(['Correcting basal friction coefficient initial guess for front inconsistencies']);
      flags = ContourToNodes(md.mesh.x,md.mesh.y,filename,2);
      %flags = md.inversion.vel_obs == 0;
      pos1 = find(flags); pos2 = find(~flags);
      %md.friction.coefficient(pos1,:) = 50;
      %md.friction.coefficient(pos1,:) = 100;
      md.friction.coefficient(pos1) = 1;
      md.inversion.max_parameters(pos1)= 1;
      %md.friction.coefficient(pos1,:) = 40;
      %md.friction.coefficient(pos1,:) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2,:),md.mesh.x(pos1),md.mesh.y(pos1));

      % % Special case: KLG
      % % Although observed velocities at the front look reasonable, a "ridge" of high friction comes out of the inversion. This is meant to
      % % neglect velocities along this erroneous ridge in the inversion.
      % if glacier == 'KLG'
      %    md.inversion.cost_functions_coefficients(pos1,1) = 0;
      %    md.inversion.cost_functions_coefficients(pos1,2) = 0;
      % end
   end

   % %Fix friction coefficient
   % filename = ['Exp/' glacier '_fixFrictionCoefficient.exp'];
   % if exist(filename,'file'),
   %    disp(['Ignoring manually selected velocities in the inversion using ' filename]);
   %    pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
   %    %md.inversion.min_parameters(pos)=md.friction.coefficient(pos);
   %    md.inversion.max_parameters(pos)=10; %md.friction.coefficient(pos);
   %    %md.inversion.cost_functions_coefficients(pos,1) = 0;
   %    %md.inversion.cost_functions_coefficients(pos,2) = 0;
   %    %md.inversion.cost_functions_coefficients(pos,3) = 0;
   % end

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;
   %md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

   % Go solve
	md.verbose.solution=1;
	md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md=solve(md,'Stressbalance');

   % Save
   savemodel(org,md);
end%}}}
if perform(org,'Relaxation'),% {{{ STEP 4

	md=loadmodel(org,'Inversion');
   
   %Put results of inversion back into the model for forward runs
   md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

   % Special case: UMI
   if strcmp(region, 'UMI')
      md_gimp = loadmodel('models/UMI_Inversion_GIMP_MEaSUREs.mat');
      %pos = md_gimp.mask.ice_levelset<0;
      %pos = ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' glacier '_Front0.exp'],1);
      %inside  = find( pos);
      %outside = find(~pos);
      %md.friction.coefficient(inside)  = md.results.StressbalanceSolution.FrictionCoefficient(inside);
      %md.friction.coefficient(outside) = md_gimp.results.StressbalanceSolution.FrictionCoefficient(outside);
      pos = ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' region '_GIMPfrictioncoefficient.exp'],1);
      inside  = find( pos);
      md.friction.coefficient(inside)  = md_gimp.results.StressbalanceSolution.FrictionCoefficient(inside);
   end
   
   % % Special cases: KLG and HEL
   % % I thought that setting the cost function coefficients to zero along the high friction "ridge" would remove the ridge that comes
   % % from the inversion. It did not. So I manually set the friction coefficient to something low here, AFTER inversion.
   % if glacier == 'KLG' | glacier == 'HEL'
   %    filename = ['Exp/' glacier '_coeffront_after_inversion.exp'];
   %    if exist(filename, 'file')
   %       pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' glacier '_coeffront_after_inversion.exp'],1));
   %       switch glacier
   %          case 'KLG'
   %             md.friction.coefficient(pos) = 10;
   %          case 'HEL'
   %             md.friction.coefficient(pos) = 40;
   %       end
   %    end
   % end

   % Special case: HEL
   % In addition to removing the "ridge" (above), increase friction along the side walls 
   % if glacier == 'HEL'
   %    filename = ['Exp/' glacier '_coefwall.exp'];
   %    if exist(filename, 'file')
   %       pos = find(ContourToNodes(md.mesh.x, md.mesh.y, filename, 1));
   %       md.friction.coefficient(pos) = 2.0 .* md.friction.coefficient(pos);
   %    end
   % end

	md.initialization.pressure = zeros(md.mesh.numberofvertices,1); 
	%md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1); ... WHY WAS THIS HERE?
	md.initialization.temperature = 250*ones(md.mesh.numberofvertices,1);

	% Set parameters
	md.inversion.iscontrol=0;
	md.timestepping.start_time = start_year;
	md.timestepping.time_step  = .02;
	md.timestepping.final_time = start_year + 1;   md.settings.output_frequency = 10; %  1 year relaxation
	%md.timestepping.final_time = start_year + 5;   md.settings.output_frequency = 10; %  5 year relaxation
	%md.timestepping.final_time = start_year + 10;  md.settings.output_frequency = (1/md.timestepping.time_step)/5; % 10 year relaxation
	%md.timestepping.final_time = start_year + 100; md.settings.output_frequency = 50; % 100 year relaxation

	% We set the transient parameters
	md.transient.ismovingfront=1;
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';

	% We set the calving model
	md.levelset.spclevelset = [md.mask.ice_levelset;0];
	md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
	md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

	%no calving
	md.levelset.spclevelset = md.levelset.spclevelset;

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume'};
	md.stressbalance.requested_outputs={'default'};

   % Go solve
	md.verbose.solution=1;
	md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md=solve(md,'transient');

   % Save
   savemodel(org,md);
end%}}}

if perform(org,'TerminusMonthly'),% {{{ STEP X
   % Create polygon to remove ice ... polygon extends from end of fjord picks to terminus
   % NOTE: Must set bufferdistance to 0, otherwise it errors
   fjord = '/Users/dfelikso/Research/Data/GlacierTermini/West_Greenland/Equip/fjord/es_fjord.shp';
   terminus = '/Users/dfelikso/Research/Data/GlacierTermini/West_Greenland/Equip/1985/es_gl_1985_03_29.shp';
   polygon = closedPolygonCreator('', terminus, fjord, 'bufferDistance', 0);
end %}}}

%%}}}

