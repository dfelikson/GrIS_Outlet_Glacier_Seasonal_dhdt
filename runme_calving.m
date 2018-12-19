%% Calving parameters -- examples {{{
if false
   % test540 - PigTranCalvingDevSSA2d
   md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.calving=calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);
   md.transient.ismovingfront = 1;
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   
   % test804 - ValleyGlacierLevelsetCalvingSIA2d
   md.transient.isgroundingline=1;
   md.transient.ismovingfront=1;
   md.calving.calvingrate=1000.*ones(md.mesh.numberofvertices,1);
   md.calving.meltingrate=zeros(md.mesh.numberofvertices,1);
   
   % test805 - ValleyGlacierLevelsetEnthCalvingHO3d
   md.transient.isgroundingline=1;
   md.transient.ismovingfront=1;
   md.calving.calvingrate=1000.*ones(md.mesh.numberofvertices,1);
   md.calving.meltingrate=zeros(md.mesh.numberofvertices,1);
   md.groundingline.melt_interpolation='SubelementMelt1';
   md.levelset.stabilization=2;
   
   % test806
   md.transient.isgroundingline=0; % b/c square ice shelf
   md.transient.isgia=0;
   md.transient.ismovingfront=1;
   md.calving=calvinglevermann();
   md.calving.coeff=4.89e13*ones(md.mesh.numberofvertices,1);
   md.calving.meltingrate=zeros(md.mesh.numberofvertices,1);
   md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
   
   % test807
   md.transient.isgroundingline=0; % b/c square ice shelf
   md.transient.isgia=0;
   md.transient.ismovingfront=1;
   md.calving.calvingrate=zeros(md.mesh.numberofvertices,1);
   md.calving.meltingrate=10000*ones(md.mesh.numberofvertices,1);
   md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
   
   % test808
   md.transient.isgroundingline=0; % b/c square ice shelf
   md.transient.isgia=0;
   md.transient.ismovingfront=1;
   md.calving=calvingminthickness();
   md.calving.min_thickness=400;
   md.calving.meltingrate=zeros(md.mesh.numberofvertices,1);
   md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
   md.levelset.reinit_frequency=1;
end
%%}}}

steps = [8];

%% Setup %%
%glacier = 'KAK'; glacier_epoch = 1985;
glacier = 'KLG'; glacier_epoch = 1999; basinNum = '494';
%glacier = 'HEL'; glacier_epoch = 1981;
region = glacier;

% To launch a job remotely and not wait, set the following:
%  md.cluster.interactive = 0;
%  md.settings.waitonlock = nan; or = 0;
clusterName = ''; % localhost
clusterName = 'oibserve';
%clusterName = 'discover';
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
      cluster.interactive = 1;
      waitonlock = 10; %nan;

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

%% Model parameters %%
org=organizer('repository',['./Models'],'prefix',[glacier '_'],'steps',steps);

%% Processing %%
set(0,'DefaultFigureWindowStyle' , 'normal')

if ~exist(['Exp/' glacier '.exp'],'file'), % {{{
   disp(['Domain outline needed: Exp/' glacier '.exp.'])
   
   s = input('Is there a shapefile that can be used for the model domain (y/n)?','s');
   if strcmpi(s(1),'y')
      [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/ModeledInlandThinning/Data/Glacier model domains/*.shp');
      shp2exp(fullfile(path,file),['Exp/' glacier '.exp']);
   else
      disp('Exiting ... rerun runme')
      return
   end
   %draw_glacier_domain(glacier, 'GIMP', 'MEaSUREs');
end
%}}}

if perform(org,'Mesh'),% {{{ STEP 1
	
   %% Mesh sizing{{{
   switch glacier
      case 'KLG'
         triangleresolution = 1000;
         hmin = 300;
         hmax = 10000;
      case 'KAK'
         triangleresolution = 300;
         hmin = 250;
         hmax = 10000;
      otherwise
         triangleresolution = 300;
         hmin = 150;
         hmax = 10000;
   end
   %%}}}
   fprintf('Creating mesh with the following parameters:\n   %20s = %5d\n   %20s = %5d\n   %20s = %5d\n', 'triangleresolution', triangleresolution, 'hmin', hmin, 'hmax', hmax);
   s = input('-> Continue (y/n)?','s');
   if ~strcmpi(s(1),'y')
      return
   end
	md=triangle(model,['./Exp/' glacier '.exp'],triangleresolution);

   % Read velocity (for mesh refinement)
	%if ~exist('vel','var'),
      disp('Reading Joughin composite velocities');
		[velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
      %disp(['Reading Joughin year ' num2str(2000) ' velocities']);
		%[velx, vely] = interpJoughin(md.mesh.x,md.mesh.y,2000);
		vel  = sqrt(velx.^2+vely.^2);
	%end

	%Adapt mesh
	disp('Optimizing mesh');
   
   %Refine mesh beyond terminus (because BAMG refinement above uses velocities that -> 0 where there's still ice in 1985)
   filename = ['Exp/' glacier '_refineFront.exp'];
   if ~exist(filename,'file'),
      s = input('Is there a shapefile that can be used to refine the front (y/n)?','s');
      if strcmpi(s(1),'y')
         [file,path] = uigetfile('/Users/denisfelikson/Research/Projects/ModeledInlandThinning/Data/Glacier model domains/*.shp');
         shp2exp(fullfile(path,file),['Exp/' glacier '_refineFront.exp']);
      else
         filename = ['Exp/' glacier '_Front0.5.exp'];
      end
   end
   
   if ~exist(filename,'file'),
      fprintf(['\n\033[' '103;30' 'm   WARNING: mesh has not been refined beyond present-day front! \033[0m \n\n']);
      disp(['  -- Refining mesh using velocities'])
      md=bamg(md,'hmin',hmin,'hmax',hmax,'field',vel,'err',2);
   else
      disp(['  -- Refining mesh using velocities and ' filename])
      hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
      in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,filename,'node',1);
      hmaxVertices(find(in))=hmin;
      md=bamg(md,'hmax',hmax,'hmin',hmin,'err',2,'field',vel,'hmaxVertices',hmaxVertices);
      % %Reload velocities (if needed later)
      % [velx, vely] = interpJoughin(md.mesh.x,md.mesh.y,2000);
      % vel  = sqrt(velx.^2+vely.^2);
   end

	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg=3413;

	savemodel(org,md);
end %}}}
if perform(org,'Param'),% {{{ STEP 2

	md=loadmodel(org,'Mesh');
	md=parameterize(md,'Par/Greenland.par');
	
   md=setflowequation(md,'SSA','all');
   
   % Weaken shear margins
   filename = ['Exp/' glacier '_shearmargins.exp'];
   if exist(filename, 'file')
      disp(['Weakening shear margins by 60% using: ' filename]);
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
      md.materials.rheology_B(pos) = 0.60 .* md.materials.rheology_B(pos);
   end

	savemodel(org,md);
end%}}}
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
   switch glacier
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
   filename = ['Exp/' glacier '_coeffront.exp'];
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

	%Go solve
	md.cluster=cluster;
   %md.settings.waitonlock = waitonlock;
   %md=solve(md,'Stressbalance','batch',batch);
   
   md=solve(md,'Stressbalance');

   %if md.clusterwaitonlock == 0 || isnan(waitonlock)
   %   loadandsaveresultsfromcluster;
   %else
      savemodel(org,md);
   %end
end%}}}
if perform(org,'Relaxation'),% {{{ STEP 4

	md=loadmodel(org,'Inversion');
   
   %Put results of inversion back into the model for forward runs
   md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

   % Special case: UMI
   if glacier == 'UMI'
      md_gimp = loadmodel('models/UMI_Inversion_GIMP_MEaSUREs.mat');
      %pos = md_gimp.mask.ice_levelset<0;
      %pos = ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' glacier '_Front0.exp'],1);
      %inside  = find( pos);
      %outside = find(~pos);
      %md.friction.coefficient(inside)  = md.results.StressbalanceSolution.FrictionCoefficient(inside);
      %md.friction.coefficient(outside) = md_gimp.results.StressbalanceSolution.FrictionCoefficient(outside);
      pos = ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' glacier '_GIMPfrictioncoefficient.exp'],1);
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
   md.timestepping.start_time = glacier_epoch;
	%md.timestepping.time_step  = .02;
	md.timestepping.time_step  = .005;
	md.timestepping.final_time = 1;   md.settings.output_frequency = 10; %  1 year relaxation
	%md.timestepping.final_time = 5;   md.settings.output_frequency = 10; %  5 year relaxation
	%md.timestepping.final_time = 10;  md.settings.output_frequency = (1/md.timestepping.time_step)/5; % 10 year relaxation
	%md.timestepping.final_time = 100; md.settings.output_frequency = 50; % 100 year relaxation

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
   %md.settings.waitonlock = waitonlock;
   %md=solve(md,'transient','batch',batch);
   md=solve(md,'transient');

   %if waitonlock == 0 || isnan(waitonlock)
   %   loadandsaveresultsfromcluster;
   %else
      savemodel(org,md);
   %end
end%}}}
if perform(org,'Transient'),% {{{ STEP 5

	md=loadmodel(org,'Relaxation');

	%Transfer resuts to model fields
	md=transientrestart(md);

	% Set parameters
	md.inversion.iscontrol=0;
	md.timestepping.time_step  = .01;
   switch glacier
      case 'KLG'
         md.timestepping.time_step  = .005;
      case 'HEL'
         md.timestepping.time_step  = .005;
   end
   md.timestepping.start_time = glacier_epoch;
	%md.timestepping.final_time = 1;   md.settings.output_frequency = 5;  %    1 year  forward run
	%md.timestepping.final_time = 10;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; %  10 years forward run
	%md.timestepping.final_time = 20;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; %  20 years forward run
	%md.timestepping.final_time = 30;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % 30 years forward run
	%md.timestepping.final_time = 50;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % 30 years forward run
   %md.timestepping.final_time = 100;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % 100 years forward run
	
   % To 2002
   md.timestepping.final_time = 2002; md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2002
   % % To 2015
   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2015
   % % To 2100
   % md.timestepping.final_time = 2100;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100

	% Go through all retreat exp files
   %%{{{
	disp('Setting up levelset time series');
	initlevelset = md.mask.ice_levelset;
	listRetreat=dir(['Exp/' glacier '_FrontRetreat*.exp']);
   if isempty(listRetreat)
      disp 'No terminus retreat exps found. Exiting.'
      return
   end
	for i=1:numel(listRetreat)
		filename = listRetreat(i).name;
		newlevelset = initlevelset;
      time = str2num(filename(numel([glacier '_FrontRetreat'])+1:numel(filename)-4));
		pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
      newlevelset(pos)=+1;
      md.levelset.spclevelset = [md.levelset.spclevelset,[newlevelset;time]];
   end

   % Go through all advance exp files
	listAdvance=dir(['Exp/' glacier '_FrontAdvance*.exp']);
	for i=1:numel(listAdvance)
		filename = listAdvance(i).name;
      time = str2num(filename(numel([glacier '_FrontAdvance'])+1:numel(filename)-4));
      % Find corresponding retreat entry (if one exists)
      idx = find(md.levelset.spclevelset(end,:)==time);
      if ~isempty(idx)
         newlevelset = md.levelset.spclevelset(1:end-1,idx);
      else
         newlevelset = initlevelset;
      end
		pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
      newlevelset(pos)=-1;
      if ~isempty(idx)
         md.levelset.spclevelset(1:end-1,idx) = newlevelset;
      else
         md.levelset.spclevelset = [md.levelset.spclevelset,[newlevelset;time]];
      end
	end

   % Sort by time
   md.levelset.spclevelset = sortrows(md.levelset.spclevelset', size(md.levelset.spclevelset,1))';

   % Check for stranded ice nodes
   [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   if sum(sum(stranded_nodes)) > 0
      fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
                 '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
      s = input('Continue by removing ice from stranded nodes (y/n)?','s');
      if strcmpi(s(1),'y')
         md.levelset.spclevelset = clean_mask_or_levelset;
      else
         return
      end
   end
   
   % Convert to signed distance fields
	if size(md.levelset.spclevelset,2)>1,
		disp('Converting levelsets to signed distance fields');
		for i=1:size(md.levelset.spclevelset,2)
			levelset = md.levelset.spclevelset(1:end-1,i);
			pos      = find(levelset<0);

			if exist('TEMP.exp','file'), delete('TEMP.exp'); end
			expcontourlevelzero(md,levelset,0,'TEMP.exp');
			levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'TEMP.exp'));
			delete('TEMP.exp'); 
			levelset(pos) = - levelset(pos);

			md.levelset.spclevelset(1:end-1,i) = levelset;
		end
	end
   %%}}}

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume'};
	%md.stressbalance.requested_outputs={'default'};
   md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 5*3600;
      fprintf('Check cluster params:\n')
      cluster
      s = input('-> Continue (y/n)?','s');
      if ~strcmpi(s(1),'y')
         return
      end
   end
	md.verbose.solution=1;
	md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md=solve(md,'transient');

   if waitonlock == 0 || isnan(waitonlock)
      loadandsaveresultsfromcluster;
   else
      savemodel(org,md);
   end
end%}}}

if perform(org,'PlotRelaxationAlongflow'),% {{{ STEP 6

	md=loadmodel(org,'Relaxation');

   %  %Create exp file along a flow line
   %  if ~exist(['Exp/' glacier '_flowline.exp'])
   %     plotmodel(md,'data',md.inversion.vel_obs,'log',10,'caxis',[1 1000])
   %     clicktoflowline(md.mesh.elements,md.mesh.x,md.mesh.y,md.inversion.vx_obs,md.inversion.vy_obs,['Exp/' glacier '_flowline.exp'])
   %  end
   [flowline, X] = glacier_flowline(glacier);
   
   % [flowline, bed, surface_matrix] = sample_bed_and_surface(glacier, md);

	%Construct matrix that holds the surface elev for all time steps available
	numt = numel(md.results.TransientSolution);
	surface_matrix = zeros(md.mesh.numberofvertices,numt);
	for i=1:numt,
		surface_matrix(:,i) = md.results.TransientSolution(i).Surface;
	end

   %%{{{ Plot -- along-flow
   figure('Position', [200, 300, 1000, 400])
	clf; set(gcf,'color','w')
   clear title
   title(glacier)

   ax1 = subplot(2,1,1);
	hold on
	surface_matrix=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,surface_matrix,flowline.x,flowline.y,  'default', nan);
	bed           =InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.bed,flowline.x,flowline.y, 'default', nan);
	map = jet(numt);
	h = [];
	for i=1:1:numt;
		%if i==1, set(gca,'FontSize',14,'XDir', 'reverse'); end
      ice_mask = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(i).MaskIceLevelset,flowline.x,flowline.y, 'default', nan);
		data_interp = surface_matrix(:,i);
		hnew=plot(X(ice_mask<0),data_interp(ice_mask<0),'-','Color',map(i,:));
      terminus = min(find(ice_mask<0));
      plot([X(terminus) X(terminus)], [bed(terminus) data_interp(terminus)], '-','Color',map(i,:));
	end
	h=[h;hnew];
	hnew=plot(X,bed,'k-');
	ylabel('s (m)','FontSize',15);
	%xlabel('Distance to front (km)','FontSize',15);
   box on; grid on;

   % Start/end surfaces on ice
   ice_mask = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).MaskIceLevelset,flowline.x,flowline.y);
   surface_start = surface_matrix(:,1);
   surface_end   = surface_matrix(:,end);
   dh = surface_end - surface_start;

   ax2 = subplot(2,1,2);
   plot(X(ice_mask<0), dh(ice_mask<0), 'k.-');
   %set(gca,'FontSize',14,'XDir', 'reverse');
	ylabel('dh (m)','FontSize',15);
	xlabel('Distance to front (km)','FontSize',15);
   box on; grid on;

   linkaxes([ax1, ax2], 'x')
   lastvalid = sum(~isnan(dh)); xlim([-1 X(lastvalid)+1]);
   
   % Save plot to pdf
	filename = ['Graphics/' glacier '_relaxation_flowline.pdf'];
	export_fig(filename,'-nocrop');
	filename = ['Graphics/' glacier '_relaxation_flowline.png'];
	export_fig(filename,'-nocrop');
   %%}}}

   %%{{{ Plot -- map view
   plotmodel(md, 'data', md.results.TransientSolution(end).Thickness - md.results.TransientSolution(1).Thickness, ...
      'mask', md.mask.ice_levelset<0, ...
      'caxis', [-50 50], ...
      'figure', 2)
   expdisp(['Exp/' glacier '_flowline.exp'], 'linestyle', 'k.-')
   axis([min(md.mesh.x) max(md.mesh.x) min(md.mesh.y) max(md.mesh.y)])
   %%}}}
end%}}}
if perform(org,'PlotAlongflow'),% {{{ STEP 7

	md=loadmodel(org,'Transient');

   % Create exp file along a flow line
   %if ~exist(['Exp/' glacier '_flowline.exp'])
   %   create_glacier_flowline_exp(glacier);
   %end
   
   % % Select the manually-picked centerline
   % if ~exist(['Exp/' glacier '_flowline.exp'])
   %    centerline_name = glacier_abbrs(glacier(1:3),'new_abbr','centerline_name');
   %    glacierCenterlineShapefile = ['/Users/denisfelikson/Research/Data/GlacierCenterlines/westGrISmanualCenterlines/' centerline_name '-centerline.shp'];
   %    shp2exp(['Exp/' glacier '_flowline.exp'], glacierCenterlineShapefile);
   % end

   % Sample bed and surface along the flowline
   [flowline, X, bed, surface_matrix] = sample_bed_and_surface(glacier, md);

   % Read observed dynamic dh and Peclet number
   filename = ['ObsMatFiles/' glacier(1:3) '.mat'];
   if exist(filename,'file')
      dh_Pe_obs = load(filename);
      
      % Special case: KLG
      if strcmp(glacier, 'KLG')
         % Shift the distance along flow for obs
         [~, idx] = min( (dh_Pe_obs.x(1) - flowline.x).^2 + (dh_Pe_obs.y(1) - flowline.y).^2 );
         dh_Pe_obs.d = dh_Pe_obs.d + X(idx)*1000;
      end
   end

   if     strcmp(glacier, 'KAK'), plot_times = [2050 2075 2100];
   elseif strcmp(glacier, 'KLG'), plot_times = [2050 2075 2100];
   elseif strcmp(glacier, 'HEL'), plot_times = [2015]; %plot_times = [2015 2050 2075 2100];
   elseif strcmp(glacier, 'UMI'), plot_times = [2005 2020 2075 2100];
   else                         , plot_times = [20 100];
   end
   
   % Setup figure
   %figure('Position', [200, 300, 1000, 600])
   figure('units','normalized','outerposition',[0 0 1 1])
	clf; set(gcf,'color','w')
   clear title;
   title(glacier)

   % Plot modeled surface evolution
   numt = numel(md.results.TransientSolution);
   ax1 = subplot(3,1,1);
	hold on
	cmap = jet(numt);
	h = [];
	%for i = 1 : ceil(1/(md.settings.output_frequency * md.timestepping.time_step)) : numt;
   for i = 1:numel(plot_times)
      data_interp = surface_matrix(:,find(times==plot_times(i)));
		%if i==1, set(gca,'FontSize',14,'XDir','reverse'); end
		if i == 1, set(gca,'FontSize',14); end
      ice_mask = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, ...
         md.results.TransientSolution(find(times==plot_times(i))).MaskIceLevelset, flowline.x, flowline.y, 'default', nan);
		hnew = plot(X(ice_mask<0), data_interp(ice_mask<0), '-', 'Color', cmap(find(times==plot_times(i)),:));
      terminus = min(find(ice_mask<0));
      plot([X(terminus) X(terminus)], [bed(terminus) data_interp(terminus)], '-','Color',cmap(find(times==plot_times(i)),:));
	end
	
   % Plot bed
   h=[h;hnew];
	hnew=plot(X,bed,'k-');
   
   colormap(ax1,'jet');
   c = colorbar(ax1);
   cticks = linspace(md.timestepping.start_time,md.timestepping.final_time,6);
   cticklabels = {};
   for i = 1:length(cticks)
      cticklabels{i} = sprintf('%4.0f',cticks(i));
   end
   c.Ticks = linmap(cticks,[0,1]);
   c.TickLabels = cticklabels;

   % Plot 1985 and WorldView surfaces (observed)
   % surface_1985 = interpAEROdem(flowline.x,flowline.y);
   % surface_1985(surface_1985 == -9999) = nan;
   % surface_1985(surface_1985 == 0.0) = nan;
   %plot(-X(ice_mask<0), surface_1985(ice_mask<0), 'k--', 'linewidth', 3)

   % surface_worldview = interpWVdem(glacier,flowline.x,flowline.y);
   % if ~isempty(surface_worldview)
   %    surface_worldview(surface_worldview == -9999) = nan;
   % end
   % if ~isempty(surface_worldview), plot(X(ice_mask<0), surface_worldview(ice_mask<0), 'k--', 'linewidth', 3), end

	ylabel('z (m)','FontSize',15);
	%xlabel('Distance to front (km)','FontSize',15);
   box on; grid on;
   set(gca,'xticklabel',[])
   %set(gca,'FontSize',14,'XDir','reverse');
   set(gca,'FontSize',14);

   % Plot dynamic dh
   % SMB anomaly
   [smb_anomaly_rate, smb_anomaly] = SMBanomaly(md, 'transient'); % [m/yr, m] ... SMB anomaly (rate and cumulative) at each forcing time from md.results.TransientSolution
   smb_anomaly_flowline = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, smb_anomaly(1:end-1,:), flowline.x, flowline.y); % interpolated along flowline
   [plot_times_mesh, X_mesh] = meshgrid(plot_times, X);
   smb_times = smb_anomaly_rate(end,:);
   smb_anomaly_interp = interp2(smb_times, X, smb_anomaly_flowline, plot_times_mesh, X_mesh); % m ... SMB anomaly interpolated to the plot times
   
   % SMB absolute
   smb_times = [md.results.TransientSolution(:).time];
   smb_absolute = [md.results.TransientSolution(:).SmbMassBalance];
   smb_cumulative = cumtrapz(smb_times, smb_absolute, 2);
   smb_cumulative_flowline = InterpFromMeshToMesh2d(md.mesh.elements, md.mesh.x, md.mesh.y, smb_cumulative, flowline.x, flowline.y);
   smb_cumulative_interp = interp2(smb_times, X, smb_cumulative_flowline, plot_times_mesh, X_mesh);

   ice_mask = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,md.results.TransientSolution(end).MaskIceLevelset,flowline.x,flowline.y);
   s_init = surface_matrix(:,1);
   times = [md.results.TransientSolution(:).time];
   plot_dh     = nan * ones(size(s_init,1),length(plot_times));
   plot_dh_dyn = nan * ones(size(s_init,1),length(plot_times));
   for i = 1:numel(plot_times)
      s = surface_matrix(:,find(times==plot_times(i)));
      plot_dh(ice_mask<0,i) = s(ice_mask<0) - s_init(ice_mask<0);
      plot_dh_dyn(:,i) = plot_dh(:,i) - smb_anomaly_interp(:,i);
   end

   % Plot dynamic thinning
   hlegend = [];
   ax2 = subplot(3,1,2); hold on
   
   % Observed
   %firstThickening = find(dh>0, 1, 'first');
   %thinningidx = ice_mask<0;
   %thinningidx(firstThickening:end) = false;
   %colormap(ax2,cm_plasma)
   %scatter(X(thinningidx), dh(thinningidx), [], unitvolthinpercent(thinningidx))
   %plot(X(firstThickening:end), dh(firstThickening:end), 'k.-', 'markersize', 15);
   h = plot(dh_Pe_obs.d/1000, dh_Pe_obs.dyn_dh, 'k-', 'LineWidth', 2);
   hlegend = [hlegend h'];
   
   % Plot average thinning between 2012-2014 (approximate dates of ArcticDEM)
   idx1 = find(times >= 2011, 1, 'first');
   idx2 = find(times <= 2014, 1, 'last');
   [~, icolor] = min(abs(2011 - times));
   surf_mean = mean(surface_matrix(:,idx1:idx2),2);
   smb_anomaly_mean = mean(smb_anomaly_flowline(:,idx1:idx2),2);
   plot_dh_mean = surf_mean - s_init;
   plot_dh_dyn_mean = plot_dh_mean - smb_anomaly_mean;
   plot_dh_dyn_mean(ice_mask>0) = nan;
   h = plot(X, plot_dh_dyn_mean, '-', 'LineWidth', 2, 'Color', cmap(icolor,:));
   hlegend = [hlegend h];

   for i = 1:numel(plot_times)
      [~, icolor] = min(abs(plot_times(i) - times));
      h = plot(X, plot_dh_dyn(:,i), '-', 'LineWidth', 2, 'Color', cmap(icolor,:));
      hlegend = [hlegend h'];
   end
	%plot(dh_Pe_obs.d/1000, dh_Pe_obs.dh, 'k.-')
   %plot(X, plot_dh, '.-')

   % Plot unit volume loss {{{
   % cm_plasma=plasma(100);
   unitvolthinpercent = cumulative_volume_thinning_percent(plot_dh_dyn);
   for i = 1:numel(plot_times)
      [~, idx] = min(abs(unitvolthinpercent(:,i) - 90));
      y1 = plot_dh_dyn(idx) - 100;
      y2 = plot_dh_dyn(idx) + 100;
      plot([X(idx) X(idx)], [y1 y2], '-', 'Color', cmap(find(times==plot_times(i)),:), 'linewidth', 2);
   end
      unitvolthinpercent_obs = cumulative_volume_thinning_percent(dh_Pe_obs.dyn_dh');
      [~, idx] = min(abs(unitvolthinpercent_obs - 90));
      y1 = dh_Pe_obs.dyn_dh(idx) - 100;
      y2 = dh_Pe_obs.dyn_dh(idx) + 100;
      plot([dh_Pe_obs.d(idx) dh_Pe_obs.d(idx)]/1000, [y1 y2], 'k-', 'linewidth', 2);
   %%}}}

   % Legend, labels, axes{{{
   legend_cell = {'observed'; 'model mean (2011-2014)'};
   for i = 1:numel(plot_times)
      legend_cell{i+2} = sprintf('model (%4d)', plot_times(i));
   end
   lgd = legend(hlegend,legend_cell,'Location','southeast');
   lgd.FontSize = 12;

   ylabel('dyn. dh (m)','FontSize',15);
	%xlabel('Distance to front (km)','FontSize',15);
   box on; grid on;
   set(gca,'xticklabel',[])
   %set(gca,'FontSize',14,'XDir','reverse');
   set(gca,'FontSize',14);
   %%}}}

   % Plot observed Pe
   ax3 = subplot(3,1,3);
   plot(dh_Pe_obs.d/1000, dh_Pe_obs.Pe, 'k.-')

   ylabel('Pe','FontSize',15);
	xlabel('along-flow distance (km)','FontSize',15);
   box on; grid on;
   %set(gca,'FontSize',14,'XDir','reverse');
   set(gca,'FontSize',14);

   linkaxes([ax1, ax2, ax3], 'x')
   lastvalid = sum(~isnan(plot_dh(:,1))); xlim([-1 X(lastvalid)+1]);
   
   set(ax1,'position',[0.13, 0.6,  0.775, 0.3])
   set(ax2,'position',[0.13, 0.35, 0.775, 0.2])
   set(ax3,'position',[0.13, 0.1,  0.775, 0.2])

   %%{{{ xlim, ylim
   xlim(ax1, [0 X(lastvalid)]);
   switch glacier
      case 'KAK'
         xlim(ax2, [-5 205]);
         ylim(ax2, [-125 25]);
      case 'KLG'
         ylim(ax2, [-250 50]);
      case 'HEL'
         %ylim(ax2, [-100 50]);
      case 'UMI'
         %xl = xlim(ax1); yl = ylim(ax1);
         %xlim(ax1,[0 50])
         %ylim(ax1,[-600 2500])
         xl = xlim(ax2); yl = ylim(ax2);
         dh_minima = nanmin(nanmin(plot_dh));
         dh_maxima = nanmax(nanmax(plot_dh));
         ylim(ax2,[dh_minima-5 dh_maxima+5])
      case 'ING'
         xl = xlim(ax1); yl = ylim(ax1);
         xlim(ax1,[0 50])
      case 'SIL'
         xl = xlim(ax1); yl = ylim(ax1);
         xlim(ax1,[0 35])
         xl = xlim(ax2); yl = ylim(ax2);
         dh_minima = [nanmin(dh_10yr) nanmin(dh_end) nanmin(dh_Pe_obs.dyn_dh)];
         dh_maxima = [nanmax(dh_10yr) nanmax(dh_end) nanmax(dh_Pe_obs.dyn_dh)];
         ylim(ax2,[nanmin(dh_minima)-5 nanmax(dh_maxima)+5])
      case 'KAN'
         xl = xlim(ax1); yl = ylim(ax1);
         xlim(ax1,[0 38])
         xl = xlim(ax2); yl = ylim(ax2);
         ylim(ax2,[nanmin(dh_end)-5 nanmax(dh_end)+5])
   end
   %}}}

   thinning_diff = compare_flowline_thinning(dh_Pe_obs.d, dh_Pe_obs.dyn_dh, X*1000, plot_dh_dyn_mean);
   fprintf('\n\tThinning mod - obs along flowline (RMSE) = %5.2f m\n\n', rms(thinning_diff(~isnan(thinning_diff))));

   % Save plot to pdf
	filename = ['Graphics/' glacier '_transient_flowline.pdf'];
	export_fig(filename) %,'-nocrop');
   return
	filename = ['Graphics/' glacier '_transient_flowline.png'];
	export_fig(filename,'-nocrop');

   %%{{{ Plot -- map view
   if false
      figure(2)
	   clf; set(gcf,'color','w')
      smb_stable = interpRACMO(md.mesh.x, md.mesh.y, 'smb_downscaled', {'stable'});
      pos1 = find( smb_stable~=0 & md.mask.ice_levelset<0 );
      pos2 = find( smb_stable==0 & md.mask.ice_levelset<0 );
      smb_stable(pos2) = griddata(md.mesh.x(pos1),md.mesh.y(pos1),smb_stable(pos1),md.mesh.x(pos2),md.mesh.y(pos2),'nearest');
      dt = mean(diff([md.results.TransientSolution(:).time]));
      time_span = md.results.TransientSolution(end).time - md.results.TransientSolution(1).time;
      smb_sum = sum([md.results.TransientSolution(:).SmbMassBalance].*dt,2);
      smb_sum_anomaly = smb_sum - smb_stable.*time_span;
      plotmodel(md, 'data', md.results.TransientSolution(end).Thickness - md.results.TransientSolution(1).Thickness - smb_sum_anomaly, ...
         'mask', md.results.TransientSolution(end).MaskIceLevelset<0, ...
         'caxis', [-50 50], ...
         'colorbarYlabel', 'meters', ...
         'title', [glacier ': dynamic dh anomaly'], ...
         'fontsize', 16, ...
         'unit', 'km', ...
         'contourlevels', {0}, 'contourcolor', 'k', 'linewidth', 3, ...
         'figure', 2, 'figposition', 'fullscreen')
      expdisp(['Exp/' glacier '_flowline.exp'], 'linestyle', 'k.-', 'multiplier', 10^-3)
      axis([min(md.mesh.x) max(md.mesh.x) min(md.mesh.y) max(md.mesh.y)]/1000)
      
      % Save plot to pdf
	   filename = ['Graphics/' glacier '_transient_DynDhAnomaly.pdf'];
	   export_fig(filename,'-nocrop');
	   filename = ['Graphics/' glacier '_transient_DynDhAnomaly.png'];
	   export_fig(filename,'-nocrop');
   end
   %%}}}

end%}}}

if perform(org,'TransientCalving'),% {{{ STEP 8

	md = loadmodel(org,'Relaxation');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep %%{{{
	md.inversion.iscontrol=0;
	md.timestepping.time_step  = .01;
   switch glacier
      case 'KLG'
         %md.timestepping.time_step  = .005;
         md.timestepping.time_step  = .02;
      case 'HEL'
         md.timestepping.time_step  = .005;
   end
   md.timestepping.start_time = glacier_epoch;
	
   % % 5 time steps and output each one
   % md.timestepping.final_time = md.timestepping.start_time + 5 * md.timestepping.time_step; md.settings.output_frequency = 1;
   % To 2000 (1 year)
   md.timestepping.final_time = 2000; md.settings.output_frequency = 1; %(1/md.timestepping.time_step)/2; % forward run to 2000
   % % To 2002
   % md.timestepping.final_time = 2002; md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2002
   % % To 2015
   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2015
   % % To 2100
   % md.timestepping.final_time = 2100;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   %%}}}

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);
   md.transient.ismovingfront = 1;
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);

   % Melt rate
   filename = ['OceanForcing/' glacier '_runoff.mat'];
   if ~exist(filename, 'file') %%{{{
      disp(['Runoff mat file ' filename ' does not exist. Exiting.'])
      return
   else
      disp('Assigning Qsg');
      runoff_matfile = load(filename);
      Qsg = [];
      for y = md.timestepping.start_time:md.timestepping.final_time-1
         for m = 1:12
            % Months from RACMO epoch (1958-01-15)
            numdays = datenum(y,m,15) - datenum(1958,1,15);
            numdaysvec = datevec(numdays);
            nummonths = numdaysvec(1) * 12 + numdaysvec(2);
            Qsg(end+1) = eval(['runoff_matfile.runoff.basin' basinNum '(nummonths + 1)']);
         end
      end
   end
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   %%}}}

   % disp('Assigning TF');
   % filename = ['Exp/' glacier '_meltrate.exp'];
   % if exist(filename, 'file')
   %    disp(['Applying melt rate using: ' filename]);
   %    pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,1));
   %    md.calving.meltingrate(pos) = 1095; % 3 m/day -> m/yr
   % end
   
   %TF = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   %TF(end,:) =  timestamps;

   %A = 3e-4;
	%B = 0.15;
	%Alpha = 0.39;
	%Beta  = 1.18;
	%BED = -repmat(min(0,md.geometry.bed),[1 numel(timestamps)]);
   %md.calving.meltingrate=zeros(md.mesh.numberofvertices+1,numel(timestamps));
   %md.calving.meltingrate(1:end-1,:) = (A*BED.*Qsg.^Alpha + B).*TF.^Beta;
   %md.calving.meltingrate(end,:) = timestamps;
   
   % % From Mathieu's AGU2016 runme.m:
   % melt = 1;
	% depth = 300;
	% meltingrate = max(0,min(melt,-md.geometry.bed/depth*melt));
	% delta_forcing = 1/12;
	% ntimes = ceil(md.timestepping.final_time/delta_forcing);
	% md.calving.meltingrate = NaN(md.mesh.numberofvertices+1,ntimes);
	% for i=1:ntimes
	% 	time = i*delta_forcing;
	% 	md.calving.meltingrate(:,i) = [meltingrate*(1+sin(2*pi*time))/2;time];
	% end
	% pos = find(md.geometry.bed>0);
	% md.calving.meltingrate(pos,:)=0;
	% pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/NoMelting.exp',2));
	% md.calving.meltingrate(pos,:)=0;

   melt_rate_max = 12; % m/day
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   meltingrate_monthly = melt_rate_max * [0 0 0 0 1/3 2/3 3/3 2/3 1/3 0 0 0];
   md.calving.meltingrate = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   idx = 1;
   for y = md.timestepping.start_time:md.timestepping.final_time-1
      for m = 1:12
         md.calving.meltingrate(1:end-1,idx) = meltingrate_monthly(m);
         idx = idx + 1;
      end
   end
	pos = find(md.geometry.bed>0);
	md.calving.meltingrate(pos,:) = 0;
   md.calving.meltingrate(1:end-1,:) = md.calving.meltingrate(1:end-1,:) * 365.25; % m/day to m/yr
   md.calving.meltingrate(end,:) = timestamps;

   % % Does the front advance?
   % md.calving.stress_threshold_floatingice = 1e50;
   % md.calving.stress_threshold_groundedice = 1e50;
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);
   
	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','CalvingCalvingrate'};
	%md.stressbalance.requested_outputs={'default'};
   %md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}
   %md.stressbalance.requested_outputs={'default','StrainRateparallel','StrainRateperpendicular','StrainRateeffective'};
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 0.5*3600;
      fprintf('Check cluster params:\n')
      cluster
      s = input('-> Continue (y/n)?','s');
      if ~strcmpi(s(1),'y')
         return
      end
   end
	md.verbose.solution=1;
	md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   %md=solve(md,'transient','batch',batch);
   md=solve(md,'transient');

   %if waitonlock == 0 || isnan(waitonlock)
   %   loadandsaveresultsfromcluster;
   %else
      savemodel(org,md);
   %end
end%}}}

