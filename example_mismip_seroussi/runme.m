addpath '/Users/dfelikso/Research/Projects/GrIS_Calibrated_SLR/issm';

% TODO
% AIS
% 1. instead of calving, crank up the melt just like it's actually done in the MISMIP+ experiment
% 2. recent dan goldberg paper showing that where shelf melt is applied is important ... review for EGU talk (possibly)

% GrIS
% 1. re-run the Transient_Steadystate_noFloating_extra to ACTUALLY get it to steady-state and then re-run the rest of the steps

steps = [8];

modelnum = 2;
switch(modelnum),
	case 1,  modelname = 'AIS';
	case 2,  modelname = 'GrIS';
	%case 3,  modelname = '500m_viscous';
	%case 4,  modelname = '250m_viscous';
	%case 5,  modelname = '2km_coulomb';
	%case 6,  modelname = '1km_coulomb';
	%case 7,  modelname = '500m_coulomb';
	%case 8,  modelname = '250m_coulomb';
	%case 9,  modelname = '125m_viscous';
	%case 10, modelname = '125m_coulomb';
end

%Hard coded parameters
loadonly = 1;
clusterName = 'oibserve';
interactive = 0;
waitonlock = nan;
switch clusterName %%{{{
   case {'','gs15serac'}
      cluster = generic('name', oshostname(), 'np', 2);
      waitonlock = 10;

   case 'oibserve'
      cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 6, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
      cluster.interactive = interactive;

   case 'discover'
      cluster=discover;
      cluster.name='discover.nccs.nasa.gov';
      cluster.login='dfelikso';
      cluster.numnodes=1;
      cluster.cpuspernode=16;
      cluster.time=1*60;
      cluster.interactive=interactive;
      cluster.processor='sand';
      cluster.queue='allnccs';
      cluster.codepath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/bin';
      cluster.executionpath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/execution';
      cluster.email='denis.felikson@nasa.gov';

end
%%}}}

%Run Steps
org=organizer('repository',['./Models_' modelname ],'prefix',['MISMIP_' modelname '_'],'steps',steps);

%Initialization
if perform(org,'Mesh'),% {{{1                                    STEP 1

   if modelnum == 1
      md=bamg(model,'domain',['Exp/Domain_' modelname '.exp'],'hmax',2000,'splitcorners',1);
   elseif modelnum == 2
      md=bamg(model,'domain',['Exp/Domain_' modelname '.exp'],'hmax',200,'splitcorners',1);
   end
	%if modelnum==1 | modelnum==5,
	%	md=bamg(model,'domain','Domain.exp','hmax',2000,'splitcorners',1);
	%elseif modelnum==2 | modelnum==6,
	%	md=bamg(model,'domain','Domain.exp','hmax',1000,'splitcorners',1);
	%elseif modelnum==3 | modelnum==7,
	%	md=bamg(model,'domain','Domain.exp','hmax',500,'splitcorners',1);
	%elseif modelnum==4 | modelnum==8,
	%	md=bamg(model,'domain','Domain.exp','hmax',250,'splitcorners',1);
	%elseif modelnum==9 | modelnum==10,
	%	md=bamg(model,'domain','Domain.exp','hmax',1000,'splitcorners',1,'maxnbv',4*10^6);
	%	hvertices=NaN*ones(md.mesh.numberofvertices,1);
	%	pos=find(md.mesh.x>200000 & md.mesh.x<800000);
	%	hvertices(pos)=500;
	%	pos=find(md.mesh.x>300000 & md.mesh.x<600000);
	%	hvertices(pos)=125;
	%	md=bamg(md,'hVertices',hvertices,'hmax',1000);
	%elseif modelnum==11,
	%	md=bamg(model,'domain','Domain.exp','hmax',500,'splitcorners',1,'maxnbv',4*10^6);
	%	hvertices=NaN*ones(md.mesh.numberofvertices,1);
	%	pos=find(md.mesh.x>300000 & md.mesh.x<600000);
	%	hvertices(pos)=125;
	%	md=bamg(md,'hVertices',hvertices,'hmax',500);
	%else
	%	error('model not supported yet');
	%end
	md.miscellaneous.name=['MISMIP_' modelname ''];

	savemodel(org,md);
end % }}}
if perform(org,'Parameterization'),% {{{1                        STEP 2

	md=loadmodel(org,'Mesh');

	md=setmask(md,'','');
	md=parameterize(md, ['MISMIP_' modelname '.par']);

	savemodel(org,md);
end% }}}
if perform(org,'Transient_Steadystate'),% {{{1                   STEP 3

	md=loadmodel(org,'Parameterization');

	md=setflowequation(md,'SSA','all');

	if(modelnum==5 | modelnum==6 | modelnum==7 | modelnum==8 | modelnum==10), %%{{{
		md.friction=frictioncoulomb();
		md.friction.coefficient=sqrt(3.160*10^6)*ones(md.mesh.numberofvertices,1); %q=1.
		md.friction.coefficientcoulomb=sqrt(0.5)*ones(md.mesh.numberofvertices,1); %q=1.
		md.friction.p=3*ones(md.mesh.numberofelements,1);
		md.friction.q=zeros(md.mesh.numberofelements,1);
	end
	if modelnum==4 | modelnum==8,
		if modelnum==4,
			md500=loadmodel('Models_500m_viscous/MISMIP_500m_viscous_Transient_steadystate8.mat');
		elseif modelnum==8,
			md500=loadmodel('Models_500m_coulomb/MISMIP_500m_coulomb_Transient_steadystate8.mat');
		end
		md.geometry.base=InterpFromMeshToMesh2d(md500.mesh.elements,md500.mesh.x,md500.mesh.y,md500.results.TransientSolution(end-1).Base,md.mesh.x,md.mesh.y);
		md.geometry.surface=InterpFromMeshToMesh2d(md500.mesh.elements,md500.mesh.x,md500.mesh.y,md500.results.TransientSolution(end-1).Surface,md.mesh.x,md.mesh.y);
		md.mask.groundedice_levelset=InterpFromMeshToMesh2d(md500.mesh.elements,md500.mesh.x,md500.mesh.y,md500.results.TransientSolution(end-1).MaskGroundediceLevelset,md.mesh.x,md.mesh.y);
		pos=find(md.mask.groundedice_levelset>0);
		md.geometry.base(pos)=md.geometry.bed(pos);
		md.geometry.thickness=md.geometry.surface-md.geometry.base;
	elseif modelnum==9 | modelnum==10,
		if modelnum==9,
			md250=loadmodel('Models_250m_viscous/MISMIP_250m_viscous_Transient_steadystate3.mat');
		elseif modelnum==10,
			md250=loadmodel('Models_250m_coulomb/MISMIP_250m_coulomb_Transient_steadystate3.mat');
		end
		md.geometry.base=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).Base,md.mesh.x,md.mesh.y);
		md.geometry.surface=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).Surface,md.mesh.x,md.mesh.y);
		md.mask.groundedice_levelset=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).MaskGroundediceLevelset,md.mesh.x,md.mesh.y);
		pos=find(md.mask.groundedice_levelset>0);
		md.geometry.base(pos)=md.geometry.bed(pos);
		md.geometry.thickness=md.geometry.surface-md.geometry.base;
	elseif modelnum==11,
		if modelnum==11,
			md250=loadmodel('Models_250m_viscous/MISMIP_250m_viscous_Transient_steadystate3.mat');
		elseif modelnum==12,
			md250=loadmodel('Models_250m_coulomb/MISMIP_250m_coulomb_Transient_steadystate3.mat');
		end
		md.geometry.base=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).Base,md.mesh.x,md.mesh.y);
		md.geometry.surface=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).Surface,md.mesh.x,md.mesh.y);
		md.mask.groundedice_levelset=InterpFromMeshToMesh2d(md250.mesh.elements,md250.mesh.x,md250.mesh.y,md250.results.TransientSolution(end-1).MaskGroundediceLevelset,md.mesh.x,md.mesh.y);
		pos=find(md.mask.groundedice_levelset>0);
		md.geometry.base(pos)=md.geometry.bed(pos);
		md.geometry.thickness=md.geometry.surface-md.geometry.base;
	end
	%%}}}
   
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

	md.timestepping.time_step=1;
	md.timestepping.final_time=50000; %200000;
	md.settings.output_frequency=500;
	
   md.timestepping.time_step=0.05;
	md.timestepping.final_time=50;
	md.settings.output_frequency=10;
   
   % Adaptive timestepping
   %md.timestepping = timesteppingadaptive(md.timestepping);
	%md.timestepping.time_step_min = 0.05;

   md.stressbalance.maxiter=10;
	md.stressbalance.abstol=NaN;
	md.stressbalance.restol=1;
	md.verbose=verbose('convergence',false,'solution',true);
	md.cluster=cluster;
   md.settings.waitonlock=waitonlock;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org, 'Transient_Steadystate_noFloating'),% {{{1       STEP 4

   md = loadmodel(org, 'Transient_Steadystate');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 100;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 50;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Start spclevelset with the original ice_levelset
   md.levelset.spclevelset          = [md.mask.ice_levelset; md.timestepping.start_time];

   % Remove all floating ice over the first year
   levelset = md.mask.groundedice_levelset;
   levelset(md.mask.groundedice_levelset<0) = +1;
   levelset(md.mask.groundedice_levelset>0) = -1;
   md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time+1];
   md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

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
         levelset(pos) = -levelset(pos);
         md.levelset.spclevelset(1:end-1,i) = levelset;
      end
   end

   % Set the requested outputs
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');

   % Save
   savemodel(org,md);

end % }}}
if perform(org, 'Transient_Steadystate_noFloating_extra'),% {{{1 STEP 5

   md = loadmodel(org, 'Transient_Steadystate_noFloating');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 100;
   md.timestepping.time_step = 0.01;
	md.settings.output_frequency=50;

   % Go solve
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');

   % Combine with Transient_Steadystate_noFloating results
   tmpResults = [md.results.TransientSolution3 md.results.TransientSolution];
   md.results.TransientSolution = tmpResults;
   md.results = rmfield(md.results, 'TransientSolution3');
   clear tmpResults

   % Save - cheat and save it back to Transient_Steadystate_noFloating
   save([org.repository '/' md.miscellaneous.name '_Transient_Steadystate_noFloating.mat'],'md','-v7.3');
   %savemodel(org,md);

end % }}}

%Experiments - specified terminus positions
experiment_style = 1;
switch experiment_style
   case 1
      retreat_3yrAvg_positions  = [26000 25000 24000];
      retreat_3yrSummer_offsets = [  500   500   500];
      retreat_3yrWinter_offsets = [  500   500   500];
   case 2
      retreat_3yrAvg_positions  = [24000 22000 20000];
      retreat_3yrSummer_offsets = [  500   500   500];
      retreat_3yrWinter_offsets = [  500   500   500];
end

if perform(org, ['Retreat_specified1_experiment' num2str(experiment_style)]),% {{{1

   md = loadmodel(org, 'Transient_Steadystate_noFloating');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 25;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   levelset = md.mask.ice_levelset;
   filename = ['Exp/' modelname '_Retreat3yrAvg_x' num2str(retreat_3yrAvg_positions(3)) '.exp'];
   pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
   levelset(pos) = +1;

   md.levelset.spclevelset          = [md.mask.ice_levelset; md.timestepping.start_time];
   md.levelset.spclevelset(:,end+1) = [            levelset; md.timestepping.start_time+3];
   
   % Remove all floating ice over the first year
   md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Convert to signed distance fields
   if size(md.levelset.spclevelset,2)>1,
      disp('Converting levelsets to signed distance fields');
      for i=1:size(md.levelset.spclevelset,2)
         levelset = md.levelset.spclevelset(1:end-1,i);
         pos      = find(levelset<0);

         if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
         expcontourlevelzero(md,levelset,0,'./TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp'));
         delete('./TEMP.exp');
         levelset(pos) = -levelset(pos);
         md.levelset.spclevelset(1:end-1,i) = levelset;
      end
   end

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Retreat_specified2_experiment' num2str(experiment_style)]),% {{{1

   md = loadmodel(org, 'Transient_Steadystate_noFloating');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 25;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Initialize
   md.levelset.spclevelset = [md.mask.ice_levelset; md.timestepping.start_time];

   % NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
   % This will setup the levelset such that the terminus continues to oscillate
   % (summer/winter) using the last positions.
   % NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE

   % % Set up seasonal levelsets
   % ds = dir(['Exp/' modelname '_Retreat3yrSummer_x' num2str(retreat_3yrAvg_position) '.exp']);
   % for i = 1:numel(ds)
   %    % Find last summer
   %    years(i) = str2num(ds(i).name(13));
   % end
   % [last_year, last_year_idx] = max(years);
   
   % Remove ice (year 3 summer)
   filename = ['Exp/' modelname '_Retreat3yrSummer_x' num2str(retreat_3yrAvg_positions(3)) '-' num2str(retreat_3yrSummer_offsets(3)) '.exp'];
   levelset = md.mask.ice_levelset;
   pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
   levelset(pos) = +1;
   md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time + 3.50];
      
   % Find corresponding winter
   filename_winter = ['Exp/' modelname '_Retreat3yrWinter_x' num2str(retreat_3yrAvg_positions(3)) '+' num2str(retreat_3yrSummer_offsets(3)) '.exp'];
   if ~exist(filename_winter,'file')
      disp('ERROR')
   else
      % Remove ice (winter)
      levelset = md.mask.ice_levelset;
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename_winter,2));
      levelset(pos) = +1;
      md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time + 4.00];
   end

   % Oscillate around the last two (summer/winter) positions until final_time
   last_spc_time = md.levelset.spclevelset(end,end);
   while last_spc_time < md.timestepping.final_time
      levelsets = md.levelset.spclevelset(1:end-1,end-1:end);
      new_times = [floor(last_spc_time) + 1.50 floor(last_spc_time) + 2.00];
      md.levelset.spclevelset = [md.levelset.spclevelset [levelsets; new_times]];
      last_spc_time = md.levelset.spclevelset(end,end);
   end

   % No calving / no melting
   md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Convert to signed distance fields
   if size(md.levelset.spclevelset,2)>1,
      disp('Converting levelsets to signed distance fields');
      for i=1:size(md.levelset.spclevelset,2)
         levelset = md.levelset.spclevelset(1:end-1,i);
         pos      = find(levelset<0);

         if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
         expcontourlevelzero(md,levelset,0,'./TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp'));
         delete('./TEMP.exp');
         levelset(pos) = -levelset(pos);
         md.levelset.spclevelset(1:end-1,i) = levelset;
      end
   end

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Retreat_specified4_experiment' num2str(experiment_style)]),% {{{1

   md = loadmodel(org, 'Transient_Steadystate_noFloating');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 25;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Initialize
   md.levelset.spclevelset = [md.mask.ice_levelset; md.timestepping.start_time];

   % NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE
   % This will setup the levelset such that the terminus continues to oscillate
   % (summer/winter) using the last positions.
   % NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE

   % Set up seasonal levelsets
   for year = 1:3
      % Remove ice (summer)
      filename = ['Exp/' modelname '_Retreat' num2str(year) 'yrSummer_x' num2str(retreat_3yrAvg_positions(year)) '-' num2str(retreat_3yrSummer_offsets(year)) '.exp'];
      levelset = md.mask.ice_levelset;
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename,2));
      levelset(pos) = +1;
      md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time + year + 0.50];

      % Find corresponding winter
      filename_winter = ['Exp/' modelname '_Retreat' num2str(year) 'yrWinter_x' num2str(retreat_3yrAvg_positions(year)) '+' num2str(retreat_3yrSummer_offsets(year)) '.exp'];
      if ~exist(filename_winter,'file')
         disp('ERROR')
      else
         % Remove ice (winter)
         levelset = md.mask.ice_levelset;
         pos = find(ContourToNodes(md.mesh.x,md.mesh.y,filename_winter,2));
         levelset(pos) = +1;
         md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time + year + 1.00];
      end
   end

   % Oscillate around the last two (summer/winter) positions until final_time
   last_spc_time = md.levelset.spclevelset(end,end);
   while last_spc_time < md.timestepping.final_time
      levelsets = md.levelset.spclevelset(1:end-1,end-1:end);
      new_times = [floor(last_spc_time) + 1.50 floor(last_spc_time) + 2.00];
      md.levelset.spclevelset = [md.levelset.spclevelset [levelsets; new_times]];
      last_spc_time = md.levelset.spclevelset(end,end);
   end

   % No calving / no melting
   md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Convert to signed distance fields
   if size(md.levelset.spclevelset,2)>1,
      disp('Converting levelsets to signed distance fields');
      for i=1:size(md.levelset.spclevelset,2)
         levelset = md.levelset.spclevelset(1:end-1,i);
         pos      = find(levelset<0);

         if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
         expcontourlevelzero(md,levelset,0,'./TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp'));
         delete('./TEMP.exp');
         levelset(pos) = -levelset(pos);
         md.levelset.spclevelset(1:end-1,i) = levelset;
      end
   end

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');

   % Save
   savemodel(org,md);

end % }}}

%Experiments - melting/calving laws
melt_rate_max = 12;
meltingrate_monthly = melt_rate_max * [0 0 0 0 1/3 2/3 3/3 2/3 1/3 0 0 0];
meltingrate_annual  = mean(meltingrate_monthly);
if perform(org,'Calving_Steadystate'),% {{{

	md = loadmodel(org,'Transient_Steadystate_noFloating');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 1;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Melt rate
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   md.calving.meltingrate = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   md.calving.meltingrate(1:end-1,:) = md.calving.meltingrate(1:end-1,:) * 365.25; % m/day to m/yr
   md.calving.meltingrate(end,:) = timestamps;

   % % Does the front advance?
   % md.calving.stress_threshold_floatingice = 1e50;
   % md.calving.stress_threshold_groundedice = 1e50;
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);
   
   % % Does the front retreat with a prescribed calving rate?
   % md.calving = calving();
   % md.calving.calvingrate = 10000*ones(md.mesh.numberofvertices,1);
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate','CalvingMeltingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}
if perform(org,['Calving_monthly_maxMeltRate' num2str(melt_rate_max)]),% {{{

	md = loadmodel(org,'Calving_Steadystate');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 100;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 10;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Melt rate
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
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
   
   % % Does the front retreat with a prescribed calving rate?
   % md.calving = calving();
   % md.calving.calvingrate = 10000*ones(md.mesh.numberofvertices,1);
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate','CalvingMeltingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}
if perform(org,['Calving_annual_maxMeltRate' num2str(melt_rate_max)]),% {{{

	md = loadmodel(org,'Calving_Steadystate');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 10;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);
   %FIXME
   pos = find(md.geometry.bed>0);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   %FIXME

   % Melt rate
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   md.calving.meltingrate = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   idx = 1;
   for y = md.timestepping.start_time:md.timestepping.final_time-1
      for m = 1:12
         md.calving.meltingrate(1:end-1,idx) = meltingrate_annual;
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
   
   % % Does the front retreat with a prescribed calving rate?
   % md.calving = calving();
   % md.calving.calvingrate = 10000*ones(md.mesh.numberofvertices,1);
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate','CalvingMeltingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}

melt_rate_pulse = 50;
melt_rate_pulse_start_time = 2.4; % melt pulse will be "inserted" at this time ...
melt_rate_pulse_duration  = 0.02; % ... and maintained for this duration
if perform(org,['Calving_monthly_maxMeltRate' num2str(melt_rate_max) '_meltPulse_time' num2str(melt_rate_pulse_start_time) '_duration' num2str(melt_rate_pulse_duration)]),% {{{

	md = loadmodel(org,'Calving_Steadystate');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 10;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   %FIXME
   pos = find(md.geometry.bed>0);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   %FIXME
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Melt rate
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   md.calving.meltingrate = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   idx = 1;
   for y = md.timestepping.start_time:md.timestepping.final_time-1
      for m = 1:12
         md.calving.meltingrate(1:end-1,idx) = meltingrate_monthly(m);
         idx = idx + 1;
      end
   end
   md.calving.meltingrate(end,:) = timestamps;

   % Add the melt pulse
   last_idx_before_pulse = find(md.calving.meltingrate(end,:)<=md.timestepping.start_time+melt_rate_pulse_start_time,1,'last');
   first_idx_after_pulse = find(md.calving.meltingrate(end,:)>=md.timestepping.start_time+melt_rate_pulse_start_time+melt_rate_pulse_duration,1,'first');
   md.calving.meltingrate(:,last_idx_before_pulse+1:first_idx_after_pulse-1) = [];
   melt_pulse = [melt_rate_pulse * ones(md.mesh.numberofvertices,2); [md.timestepping.start_time+melt_rate_pulse_start_time md.timestepping.start_time+melt_rate_pulse_start_time+melt_rate_pulse_duration]];
   md.calving.meltingrate = [md.calving.meltingrate(:,1:last_idx_before_pulse) melt_pulse md.calving.meltingrate(:,first_idx_after_pulse:end)];
	
   pos = find(md.geometry.bed>0);
	md.calving.meltingrate(pos,:) = 0;
   md.calving.meltingrate(1:end-1,:) = md.calving.meltingrate(1:end-1,:) * 365.25; % m/day to m/yr

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate','CalvingMeltingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}
if perform(org,['Calving_annual_maxMeltRate' num2str(melt_rate_max) '_meltPulse_time' num2str(melt_rate_pulse_start_time) '_duration' num2str(melt_rate_pulse_duration)]),% {{{

	md = loadmodel(org,'Calving_Steadystate');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 10;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Melt rate
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   md.calving.meltingrate = zeros(md.mesh.numberofvertices+1,numel(timestamps));
   idx = 1;
   for y = md.timestepping.start_time:md.timestepping.final_time-1
      for m = 1:12
         md.calving.meltingrate(1:end-1,idx) = meltingrate_annual;
         idx = idx + 1;
      end
   end
   md.calving.meltingrate(end,:) = timestamps;
   
   % Add the melt pulse
   last_idx_before_pulse = find(md.calving.meltingrate(end,:)<=md.timestepping.start_time+melt_rate_pulse_start_time,1,'last');
   first_idx_after_pulse = find(md.calving.meltingrate(end,:)>=md.timestepping.start_time+melt_rate_pulse_start_time+melt_rate_pulse_duration,1,'first');
   md.calving.meltingrate(:,last_idx_before_pulse+1:first_idx_after_pulse-1) = [];
   melt_pulse = [melt_rate_pulse * ones(md.mesh.numberofvertices,2); [md.timestepping.start_time+melt_rate_pulse_start_time md.timestepping.start_time+melt_rate_pulse_start_time+melt_rate_pulse_duration]];
   md.calving.meltingrate = [md.calving.meltingrate(:,1:last_idx_before_pulse) melt_pulse md.calving.meltingrate(:,first_idx_after_pulse:end)];
	
	pos = find(md.geometry.bed>0);
	md.calving.meltingrate(pos,:) = 0;
   md.calving.meltingrate(1:end-1,:) = md.calving.meltingrate(1:end-1,:) * 365.25; % m/day to m/yr

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate','CalvingMeltingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}

return
if perform(org,'TransientCalving'),% {{{

	md = loadmodel(org,'Transient_Steadystate_noFloating');

	% Transfer resuts to model fields
	md = transientrestart(md);

   % Set start/end times and timestep
   md.timestepping.final_time = md.timestepping.start_time + 1;
   md.timestepping.time_step = 0.01;
   md.settings.output_frequency = 1;

   % Calving modes (code from test540 - PigTranCalvingDevSSA2d)
   disp('Setting calving model')
   md.transient.isgroundingline = 1;
   md.transient.ismovingfront = 1;
   %md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
   md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
   pos = find(md.mesh.vertexonboundary);
   md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);
   md.calving = calvingvonmises();
   md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

   %% Melt rate %%{{{
   %filename = ['OceanForcing/' glacier '_runoff.mat'];
   %if ~exist(filename, 'file')
   %   disp(['Runoff mat file ' filename ' does not exist. Exiting.'])
   %   return
   %else
   %   disp('Assigning Qsg');
   %   runoff_matfile = load(filename);
   %   Qsg = [];
   %   for y = md.timestepping.start_time:md.timestepping.final_time-1
   %      for m = 1:12
   %         % Months from RACMO epoch (1958-01-15)
   %         numdays = datenum(y,m,15) - datenum(1958,1,15);
   %         numdaysvec = datevec(numdays);
   %         nummonths = numdaysvec(1) * 12 + numdaysvec(2);
   %         Qsg(end+1) = eval(['runoff_matfile.runoff.basin' basinNum '(nummonths + 1)']);
   %      end
   %   end
   %end
   %timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];

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
   %%}}}

   melt_rate_max = 12; % m/day
   %md.levelset.calving_max = 10000; % 10 km/yr
   timestamps = [md.timestepping.start_time : 1/12 : md.timestepping.final_time-1/12];
   meltingrate_monthly = melt_rate_max * [0 0 0 0 1/3 2/3 3/3 2/3 1/3 0 0 0];
   meltingrate_monthly = mean(meltingrate_monthly) * ones(1,12);
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
   
   % % Does the front retreat with a prescribed calving rate?
   % md.calving = calving();
   % md.calving.calvingrate = 10000*ones(md.mesh.numberofvertices,1);
   % md.calving.meltingrate = zeros(md.mesh.numberofvertices,1);

	% Set the requested outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','CalvingCalvingrate'};
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
   md=solve(md,'transient');

   savemodel(org,md);
end%}}}

% TODO TRY THIS
md.timestepping = timesteppingadaptive(md.timestepping);

% Other stuff I was trying
if perform(org, 'Inversion_specifySurface'),% {{{1

   md = loadmodel(org, 'Transient_Steadystate_noFloating');

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.geometry.surface = md.geometry.surface + 50;
   md.geometry.thickness = md.geometry.thickness + 50;
   md.timestepping.final_time = md.timestepping.start_time+1;

	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);

   % Obs
   md.inversion.vx_obs = md.initialization.vx;
   md.inversion.vy_obs = md.initialization.vy;
   md.inversion.vel_obs = md.initialization.vel;

	%Cost functions
	md.inversion.cost_functions=[101 103 501]; %Abs, Log, reg
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=40;
   md.inversion.cost_functions_coefficients(:,3)=1.6e-06;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=50;
	md.inversion.maxiter =50;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;
   md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy'}

   % Go solve
	md.verbose.solution=1;
	md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'Stressbalance');

   % Save
   savemodel(org,md);

end % }}}
if perform(org,'Transient_Inversion_specifySurface'),% {{{1

	md=loadmodel(org,'Inversion_specifySurface');
   
   %Put results of inversion back into the model for forward runs
   md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;

   % Set parameters
   md.inversion.iscontrol=0;
   md.timestepping.time_step  = .01;
   md.timestepping.final_time = md.timestepping.start_time + 10;
   md.settings.output_frequency = (1/md.timestepping.time_step)/5; % 5/yr

	% Set the requested outputs
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Solve
	md.cluster=cluster;
   md.settings.waitonlock=waitonlock;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}

%Experiments SEP1
if perform(org,'Experiment0'),% {{{1

	md=loadmodel(org,'Transient_Steadystate');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(9).Vx;
	md.initialization.vy=md.results.TransientSolution(9).Vy;
	md.initialization.vel=md.results.TransientSolution(9).Vel;
	md.geometry.thickness=md.results.TransientSolution(9).Thickness;
	md.geometry.base=md.results.TransientSolution(9).Base;
	md.geometry.surface=md.results.TransientSolution(9).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(9).MaskGroundediceLevelset;

	md.timestepping.time_step=0.25;
	md.timestepping.final_time=200;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=10;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.cluster=cluster;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP0'];
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1r'),% {{{1

	md=loadmodel(org,'Experiment0');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(101).Vx;
	md.initialization.vy=md.results.TransientSolution(101).Vy;
	md.initialization.vel=md.results.TransientSolution(101).Vel;
	md.geometry.thickness=md.results.TransientSolution(101).Thickness;
	md.geometry.base=md.results.TransientSolution(101).Base;
	md.geometry.surface=md.results.TransientSolution(101).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(101).MaskGroundediceLevelset;

	md.basalforcings.meltrate_factor=0.2;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=100;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};

	md.miscellaneous.name=['MISMIP_' modelname '_MP1r'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	md.results=[];
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1ra_short'),% {{{1

	md=loadmodel(org,'Experiment1r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings.meltrate_factor=0.0;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=10;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP1ra'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1ra'),% {{{1

	md=loadmodel(org,'Experiment1r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings.meltrate_factor=0.0;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=900;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=10;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP1ra'];
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1rr'),% {{{1

	md=loadmodel(org,'Experiment1r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings.meltrate_factor=0.2;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=900;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP1rr'];
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'Experiment1rr_short'),% {{{1

	md=loadmodel(org,'Experiment1r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings.meltrate_factor=0.2;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP1rr'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
if perform(org,'Experiment2r'),% {{{1

	md=loadmodel(org,'Experiment0');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end-1).Vx;
	md.initialization.vy=md.results.TransientSolution(end-1).Vy;
	md.initialization.vel=md.results.TransientSolution(end-1).Vel;
	md.geometry.thickness=md.results.TransientSolution(end-1).Thickness;
	md.geometry.base=md.results.TransientSolution(end-1).Base;
	md.geometry.surface=md.results.TransientSolution(end-1).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end-1).MaskGroundediceLevelset;

	md.basalforcings=basalforcings();
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	pos=find(md.mesh.x>480000);
	md.basalforcings.floatingice_melting_rate(pos)=100;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};

	md.miscellaneous.name=['MISMIP_' modelname '_MP2r'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
if perform(org,'Experiment2ra_short'),% {{{1

	md=loadmodel(org,'Experiment2r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings=basalforcings();
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP2ra'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
if perform(org,'Experiment2ra'),% {{{1

	md=loadmodel(org,'Experiment2r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings=basalforcings();
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=900;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP2ra'];
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'Experiment2rr'),% {{{1

	md=loadmodel(org,'Experiment2r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings=basalforcings();
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	pos=find(md.mesh.x>480000);
	md.basalforcings.floatingice_melting_rate(pos)=100;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=900;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP2rr'];
	md.cluster=cluster;
	md=solve(md,'tr');

	savemodel(org,md);
end% }}}
if perform(org,'Experiment2rr_short'),% {{{1

	md=loadmodel(org,'Experiment2r');

	md=setflowequation(md,'SSA','all');
	md.initialization.vx=md.results.TransientSolution(end).Vx;
	md.initialization.vy=md.results.TransientSolution(end).Vy;
	md.initialization.vel=md.results.TransientSolution(end).Vel;
	md.geometry.thickness=md.results.TransientSolution(end).Thickness;
	md.geometry.base=md.results.TransientSolution(end).Base;
	md.geometry.surface=md.results.TransientSolution(end).Surface;
	md.mask.groundedice_levelset=md.results.TransientSolution(end).MaskGroundediceLevelset;

	md.basalforcings=basalforcings();
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	pos=find(md.mesh.x>480000);
	md.basalforcings.floatingice_melting_rate(pos)=100;
	md.timestepping.time_step=0.25;
	md.timestepping.final_time=100;
	md.settings.output_frequency=4;
	md.stressbalance.maxiter=20;
	md.stressbalance.restol=1;
	md.stressbalance.reltol=0.001;
	md.stressbalance.abstol=NaN;
	md.verbose.convergence=true;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','BasalforcingsFloatingiceMeltingRate','GroundedArea'};
	md.miscellaneous.name=['MISMIP_' modelname '_MP2rr'];
	md.settings.waitonlock = 0;
	md.cluster = cluster;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	savemodel(org,md);
end% }}}
