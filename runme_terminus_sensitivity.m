% FIXME
% [x] 1. change the time of the spclevelset in relaxation to start_year -- not year 0
% [x] 2. re-run to end_year=2015
% [x] 3. need to make to "end_year"s -- one that specifies the end year of the terminus positions (i.e., the hindcast?) and one that specifies the end year of the simulation (i.e., the forecast?)
% [ ] 4. there is ice in the BedMachine mask that lies beyond the Front0 exp file -- remove it somehow
% [ ] 5. need separate MinExtent exps ... one for yearly termini ... one for monthly

steps = [12];
clusterName = ''; % empty for localhost
clusterName = 'oibserve';
%clusterName = 'discover';

%% Setup
region = 'SAtoES';
experiment = 'NyearlyTest';
switch region
   case 'WestGrIS' %%{{{
      start_year = 1985;
      end_year   = 2010;
      relaxation_years = 1;
      termini_dir = '/Users/dfelikso/Research/Data/GlacierTermini/West_Greenland';
      clear glaciers
      glaciers{ 1} = 'Ing';
      glaciers{ 2} = 'Umiamiko';
      glaciers{ 3} = 'RINK';
      glaciers{ 4} = 'KS';
      glaciers{ 5} = 'Kss';
      glaciers{ 6} = 'Perd';
      glaciers{ 7} = 'Silvia';
      glaciers{ 8} = 'Kang';
      glaciers{ 9} = 'Sermilik';
      glaciers{10} = 'Lille';
      glaciers{11} = 'Store';
      glaciers{12} = 'SA';
      glaciers{13} = 'SK';
      glaciers{14} = 'Kan';
      glaciers{15} = 'Equip';
      %%}}}
   case 'SAtoES' %%{{{
      start_year = 1985;
      end_year   = 2100;
      relaxation_years = 1;
      termini_dir = '/Users/dfelikso/Research/Data/GlacierTermini/West_Greenland';
      end_year_termini = 2015;
      clear glaciers
      glaciers{ 1} = 'SA';
      glaciers{ 2} = 'SK';
      glaciers{ 3} = 'Kan';
      glaciers{ 4} = 'Equip';
%%}}}
end

% To launch a job remotely and not wait, set the following:
%  md.cluster.interactive = 0;
%  md.settings.waitonlock = nan; or = 0;
switch clusterName %%{{{
   case {'','gs15serac'}
      cluster = generic('name', oshostname(), 'np', 2);
      waitonlock = 10;

   case 'oibserve'
      cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 8, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution', ...
         'valgrind', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/install/bin/valgrind', ...
         'valgrindlib', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/install/lib/valgrind/libmpiwrap-amd64-linux.so', ...
         'valgrindsup', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/issm.supp');
      cluster.interactive = 0; %1;
      waitonlock = nan; %90;

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
model_repository_dir = ['./Models/' region '_' experiment];
if isdir(model_repository_dir) %%{{{
   d = dir(model_repository_dir);
   if length(d) > 2
      fprintf(['\n\033[' '103;30' 'mWARNING: Repository directory ' model_repository_dir ' is not empty. \033[0m \n']);
      s = input('Overwrite files (y/n)?','s');
      if ~strcmpi(s(1),'y')
         disp('Exiting ... rerun runme')
         return
      end
   end
else
   mkdir(model_repository_dir);
end
%%}}}
org = organizer('repository', model_repository_dir, 'prefix', [region '_'], 'steps', steps);

%% Processing
if ~exist(['Exp/' region '.exp'],'file') %%{{{
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
%%}}}

if perform(org,'ExpSetup'),% {{{ STEP 1

   clear year
   fprintf('Setting up terminus exps\n');
   for iglacier = 1:numel(glaciers)
      % Names
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      fprintf(' -> %s\n', glacier_long_name);

      % Remove "Front0" and "TerminusYearly" exp files for the glaciers we're simulating
      filename = sprintf('Exp/%s_Front0.exp',            glacier_long_name); if length(dir(filename)) > 0; delete(filename); end;
      filename = sprintf('Exp/%s_TerminusMonthly_*.exp', glacier_long_name); if length(dir(filename)) > 0; delete(filename); end;
      filename = sprintf('Exp/%s_TerminusYearly_*.exp',  glacier_long_name); if length(dir(filename)) > 0; delete(filename); end;
      filename = sprintf('Exp/%s_TerminusMinExtent.exp', glacier_long_name); if length(dir(filename)) > 0; delete(filename); end;
      
      termini = [];

      % Names
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');
      
      % Find all shapefiles
      parent_dir = ['/Users/dfelikso/Research/Data/GlacierTermini/West_Greenland/' glacier_name];
      fjord = [parent_dir '/fjord/' glacier_abbr '_fjord.shp'];
      shps = dir([parent_dir '/*/*shp']);
      for i = 1:numel(shps)
         parts = strsplit(shps(i).name,'_');
         if numel(parts) == 5
            y = str2int( parts{3} );
            m = str2int( parts{4} );
            d = str2int( parts{5}(1:2) );
            if isempty(termini)
               termini.filename = [shps(i).folder '/' shps(i).name];
               termini.dt = datetime(y,m,d);
            else
               termini(end+1).filename = [shps(i).folder '/' shps(i).name];
               termini(end  ).dt = datetime(y,m,d);
            end
         end
      end

      % Catalog all termini (TerminusMonthly) %%{{{
      % Find termini closest to the 15th of every month (of every year)
      dts = [termini(:).dt];
      idxs_select = [];
      for y = start_year:end_year_termini
         year_dir = [parent_dir '/' sprintf('%d', y)];
         for m = 1:12
            [~, idx] = min(abs( datetime(y,m,15) - dts ));
            if year(dts(idx)) ~= y | month(dts(idx)) ~= m
               fprintf('    Warning! Terminus not found for %s for year %d, month %d.\n', glacier_long_name, y, m);
            else
               %fprintf('Terminus found: %s.\n', datestr(dts(idx)));
               if isempty(idxs_select) idxs_select = idx; else idxs_select(end+1) = idx; end
            end
         end
      end

      % Create polygon to remove ice ... polygon extends from end of fjord picks to terminus
      % NOTE: Must set bufferdistance to 0, otherwise it errors
      for idx = idxs_select
         terminus = termini(idx).filename; %[parent_dir '/' datestr(dt(idx), 'yyyy') '/' glacier_abbr '_gl_' datestr(dt(idx), 'yyyy_mm_dd') '.shp'];
         try
            polygon = closedPolygonCreator(terminus, '', fjord, 'bufferDistance', 0, 'utm2stere', true);
            a.x = polygon(:,1);
            a.y = polygon(:,2);
            expwrite(a, ['Exp/' glacier_long_name '_TerminusMonthly_' datestr(dts(idx), 'yyyymmm') '.exp']);
         catch
            fprintf(' Some kind of error encountered in closedPolygonCreator while making %s.\n', ['Exp/' glacier_long_name '_TerminusYearly_' datestr(dts(idx), 'yyyymmmdd') '.exp']);
         end
      end %%}}}

      % Catalog all termini (TerminusYearly) %%{{{
      shps = dir([parent_dir '/*/*shp']);
      for i = 1:numel(shps)
         parts = strsplit(shps(i).name,'_');
         if numel(parts) == 5
            y = str2int( parts{3} );
            m = str2int( parts{4} );
            d = str2int( parts{5}(1:2) );
            if isempty(termini)
               termini.filename = [shps(i).folder '/' shps(i).name];
               termini.dt = datetime(y,m,d);
            else
               termini(end+1).filename = [shps(i).folder '/' shps(i).name];
               termini(end  ).dt = datetime(y,m,d);
            end
         end
      end

      % Find termini closest to August 15 of every year
      dts = [termini(:).dt];
      idxs_select = [];
      for y = start_year:end_year_termini
         year_dir = [parent_dir '/' sprintf('%d', y)];
         [~, idx] = min(abs( datetime(y,8,15) - dts ));
         if year(dts(idx)) ~= y
            fprintf('    Warning! Terminus not found for %s for year %d.\n', glacier_long_name, y);
         else
            %fprintf('Terminus found: %s.\n', datestr(dts(idx)));
            if isempty(idxs_select) idxs_select = idx; else idxs_select(end+1) = idx; end
         end
      end

      % Create polygon to remove ice ... polygon extends from end of fjord picks to terminus
      % NOTE: Must set bufferdistance to 0, otherwise it errors
      for idx = idxs_select
         terminus = termini(idx).filename; %[parent_dir '/' datestr(dt(idx), 'yyyy') '/' glacier_abbr '_gl_' datestr(dt(idx), 'yyyy_mm_dd') '.shp'];
         try
            polygon = closedPolygonCreator(terminus, '', fjord, 'bufferDistance', 0, 'utm2stere', true);
            a.x = polygon(:,1);
            a.y = polygon(:,2);
            expwrite(a, ['Exp/' glacier_long_name '_TerminusYearly_' datestr(dts(idx), 'yyyy') '.exp']);
         catch
            fprintf(' Some kind of error encountered in closedPolygonCreator while making %s.\n', ['Exp/' glacier_long_name '_TerminusYearly_' datestr(dts(idx), 'yyyymmmdd') '.exp']);
         end
      end %%}}}

      % Create exps to refine mesh within furthest extent of ice (Front0) %%{{{
      % Terminus pick closest to August 15 [start_year]
      [~, idx] = min(abs( datetime(start_year,08,15) - dts ));
      terminus1 = termini(idx).filename;

      fprintf('    %30s Front0: %s\n', glacier_long_name, datestr(dts(idx),'yyyymmmdd'));
      % Create a polygon and extend the ice mask
      polygon = closedPolygonCreator(terminus1, '', fjord, 'bufferDistance', 0, 'utm2stere', true);
      a.x = polygon(:,1);
      a.y = polygon(:,2);
      expwrite(a, sprintf('Exp/%s_Front0.exp', glacier_long_name)); %%}}}

      % Find minimum ice extent used to remove ice from original mask to initialize the levelset (TerminusMinExtent) %%{{{
      area_min_ice_extent = nan;
      year_min_ice_extent = nan;
      for y = start_year:end_year
         exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if numel(exps) == 1
            exp_struct = expread([exps(1).folder '/' exps(1).name]);
            exp_area = polyarea(exp_struct.x, exp_struct.y);
            if isnan(area_min_ice_extent) || exp_area < area_min_ice_extent
               year_min_ice_extent = exps(1).name(numel([glacier_long_name '_TerminusYearly_'])+1:numel(exps(1).name)-4);
               area_min_ice_extent = exp_area;
            end
         elseif numel(exps) > 1
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end

      % Write new exp with the area of ice to be removed for min extent
      [~, idx] = min(abs( datetime(str2num(year_min_ice_extent),08,15) - dts ));
      terminus = dir(sprintf('%s/%s/%s/%s_gl_%s*.shp',termini_dir,glacier_name,datestr(dts(idx),'yyyy'),glacier_abbr,datestr(dts(idx),'yyyy_mm_dd')));
      if numel(terminus) == 1
         terminus = [terminus(1).folder '/' terminus(1).name];
      else
         fprintf(' Found none or multiple shapefiles for glacier %s, datestring %4d\n', glacier_long_name, datestr(dts(idx),'yyyy_mm_dd'));
      end

      fjord = sprintf('%s/%s/fjord/%s_fjord.shp', termini_dir, glacier_name, glacier_abbr);
      try
         polygon = closedPolygonCreator('', terminus, fjord, 'bufferDistance', 0, 'utm2stere', true);
         a.x = polygon(:,1);
         a.y = polygon(:,2);
         expwrite(a, ['Exp/' glacier_long_name '_TerminusMinExtent.exp']);
      catch
         fprintf(' Some kind of error encountered in closedPolygonCreator while making %s.\n', ['Exp/' glacier_long_name '_TerminusMinExtent.exp']);
      end %%}}}
   end

end %}}}
if perform(org,'Mesh'),% {{{ STEP 2
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
   disp(['  -- Refining mesh using velocities and ... '])
   hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
   for iglacier = 1:numel(glaciers)
      % Names
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      d = dir(['Exp/' glacier_long_name '_Front0.exp']);
      if numel(d) == 1
         filename = [d.folder '/' d.name];
         if ~exist(filename,'file'),
         end
         disp(['    -> ' filename])
         elem_in = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,filename,'element',1));
         elem_adj = reshape(md.mesh.elementconnectivity(elem_in,:), 3*numel(elem_in), 1);
         elem_adj = reshape(md.mesh.elementconnectivity(elem_adj,:), 3*numel(elem_adj), 1);
         %node_in = reshape(md.mesh.elements(elem_in,:), 3*numel(elem_in), 1);
         node_adj = reshape(md.mesh.elements(elem_adj,:), 3*numel(elem_adj), 1);
         hmaxVertices(node_adj)=hmin;
      else
         error(['    PROBLEM for glacier ' glacier_long_name])
      end
   end
   
   md = bamg(md,'hmin',hmin,'hmax',hmax,'field',vel,'err',2,'hmaxVertices',hmaxVertices);
   %md = bamg(md,'hmin',hmin,'hmax',hmax,'field',vel,'err',2);
   
   [md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg = 3413;

	savemodel(org,md);

end %}}}
if perform(org,'Param'),% {{{ STEP 3

	md=loadmodel(org,'Mesh');

   % Parameterize
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
if perform(org,'Inversion'),% {{{ STEP 4

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
   md.inversion.cost_functions_coefficients(:,3)=1.6e-06;
   switch region
      case 'KAK'
         md.inversion.cost_functions_coefficients(:,3)=1.6e-06;
      case 'KLG'
         md.inversion.cost_functions_coefficients(:,3)=8.0e-05;
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
      md.friction.coefficient(pos1,:)  = 10;
      %md.friction.coefficient(pos1,:)  = 40;
      %md.friction.coefficient(pos1,:)  = 50;
      %md.friction.coefficient(pos1,:)  = 100;
      %md.friction.coefficient(pos1,:)  = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2,:),md.mesh.x(pos1),md.mesh.y(pos1));

      md.inversion.max_parameters(pos1) = md.friction.coefficient(pos1,:);
      
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
   md = solve(md,'Stressbalance');

   % Save
   savemodel(org,md);
end%}}}
if perform(org,['Relaxation' num2str(relaxation_years) 'yr']),% {{{ STEP 5

	md=loadmodel(org,'Inversion');
   
   %Put results of inversion back into the model for forward runs
   md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

   % Special post-processing of inverted friction coefficient
   filename = ['Exp/' region '_coeffront_after_inversion.exp'];
   if exist(filename, 'file')
      pos = find(ContourToNodes(md.mesh.x, md.mesh.y, filename, 1));
      md.friction.coefficient(pos) = 10;
   end

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

   % Save
   savemodel(org,md);
end%}}}

if perform(org,'Control'),% {{{ STEP 6

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % No moving front but need to turn on grounding line migration
	md.transient.isgroundingline=1;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
if perform(org,'TerminusYearly'),% {{{ STEP 7

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   if false
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   % 2.) Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   for y = start_year:end_year_termini
      levelset = levelset_min_extent;
      %fprintf(' Finding termini for year %4d\n', year);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusYearly retreat exps found for %s\n', glacier_long_name);
         elseif numel(exps) == 1
	      	filename = exps(1).name;
	      	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end
      if isnan(md.levelset.spclevelset) 
         md.levelset.spclevelset = [levelset; y];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y];
      end
   end %%}}}
   end
   % Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   levelset = md.mask.ice_levelset;
   for y = start_year:end_year_termini
      %fprintf(' Finding termini for year %4d\n', y);
      %fprintf(' Finding termini for month %2d\n', m);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusMonthly retreat exps found for %s, year %4d, month %s\n', glacier_long_name, y, month_str);
         elseif length(exps) == 1
            % 1.) Remove ice down to the minimum extent for this glacier
            exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
            levelset(pos) = +1;

            % 2.) Add ice using the extent of this month's exp file
            filename = exps(1).name;
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('PROBLEM\n');
            return
         end
      end

      % Populate spclevelset
      [doy,frac] = date2doy(datenum(y,8,15));
      if isnan(md.levelset.spclevelset)
         md.levelset.spclevelset = [levelset; y+frac];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y+frac];
      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
if perform(org,['Terminus' num2str(start_year) '-' num2str(end_year_termini)]),% {{{ STEP 8

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   % 2.) Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   for y = [start_year end_year_termini]
      levelset = levelset_min_extent;
      %fprintf(' Finding termini for year %4d\n', year);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusYearly retreat exps found for %s\n', glacier_long_name);
         elseif numel(exps) == 1
	      	filename = exps(1).name;
	      	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end
      if isnan(md.levelset.spclevelset) 
         md.levelset.spclevelset = [levelset; y];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y];
      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
if perform(org,['TerminusMonthly']),% {{{ STEP 9

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   if false
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   end
   % Use monthly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   levelset = md.mask.ice_levelset;
   for y = start_year:end_year_termini
      %fprintf(' Finding termini for year %4d\n', y);
      for m = 1:12
         levelset_before = levelset; filenames = {};%DEBUG
         %fprintf(' Finding termini for month %2d\n', m);
         month_str = month(datetime(2000,m,1), 'shortName');
         month_str = month_str{1};
         for iglacier = 1:numel(glaciers)
            glacier_name = glaciers{iglacier};
            glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
            glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	         exps = dir(sprintf('Exp/%s_TerminusMonthly_%d%s.exp', glacier_long_name, y, month_str));
            if isempty(exps)
               %fprintf('No TerminusMonthly retreat exps found for %s, year %4d, month %s\n', glacier_long_name, y, month_str);
            elseif length(exps) == 1
               % 1.) Remove ice down to the minimum extent for this glacier
               exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
               pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
               levelset(pos) = +1;

               % 2.) Add ice using the extent of this month's exp file
               filename = exps(1).name;
               pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
               levelset(pos) = -1;
               filenames{end+1} = ['Exp/' filename]; %DEBUG
            else
               fprintf('PROBLEM\n');
               return
            end
         end

         % Populate spclevelset
         [doy,frac] = date2doy(datenum(y,m,15));
         if isnan(md.levelset.spclevelset)
            md.levelset.spclevelset = [levelset; y+frac];
         else
            md.levelset.spclevelset(:,end+1) = [levelset; y+frac];
         end
         
         %DEBUG
         %if ~isempty(filenames)
         %   levelset_after_exp = levelset;
         %   plotmodel(md, 'data', levelset_before, 'figure', 1, 'title', 'before', 'figposition', 'fullscreen')
         %   for i = 1:numel(filenames) expdisp(filenames{i}); end
         %   plotmodel(md, 'data', levelset_after_exp, 'figure', 2, 'title', sprintf('%f',y+frac), 'figposition', 'fullscreen')
         %   for i = 1:numel(filenames) expdisp(filenames{i}); end
         %   filenames
         %   pause
         %end
         %DEBUG

      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}

if perform(org,'Terminus3yearly'),% {{{ STEP 10

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   if false
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   % 2.) Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   for y = start_year:end_year_termini
      levelset = levelset_min_extent;
      %fprintf(' Finding termini for year %4d\n', year);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusYearly retreat exps found for %s\n', glacier_long_name);
         elseif numel(exps) == 1
	      	filename = exps(1).name;
	      	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end
      if isnan(md.levelset.spclevelset) 
         md.levelset.spclevelset = [levelset; y];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y];
      end
   end %%}}}
   end
   % Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   levelset = md.mask.ice_levelset;
   for y = start_year:3:end_year_termini
      %fprintf(' Finding termini for year %4d\n', y);
      %fprintf(' Finding termini for month %2d\n', m);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusMonthly retreat exps found for %s, year %4d, month %s\n', glacier_long_name, y, month_str);
         elseif length(exps) == 1
            % 1.) Remove ice down to the minimum extent for this glacier
            exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
            levelset(pos) = +1;

            % 2.) Add ice using the extent of this month's exp file
            filename = exps(1).name;
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('PROBLEM\n');
            return
         end
      end

      % Populate spclevelset
      [doy,frac] = date2doy(datenum(y,8,15));
      if isnan(md.levelset.spclevelset)
         md.levelset.spclevelset = [levelset; y+frac];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y+frac];
      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
if perform(org,'Terminus5yearly'),% {{{ STEP 11

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   if false
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   % 2.) Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   for y = start_year:end_year_termini
      levelset = levelset_min_extent;
      %fprintf(' Finding termini for year %4d\n', year);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusYearly retreat exps found for %s\n', glacier_long_name);
         elseif numel(exps) == 1
	      	filename = exps(1).name;
	      	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end
      if isnan(md.levelset.spclevelset) 
         md.levelset.spclevelset = [levelset; y];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y];
      end
   end %%}}}
   end
   % Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   levelset = md.mask.ice_levelset;
   for y = start_year:5:end_year_termini
      %fprintf(' Finding termini for year %4d\n', y);
      %fprintf(' Finding termini for month %2d\n', m);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusMonthly retreat exps found for %s, year %4d, month %s\n', glacier_long_name, y, month_str);
         elseif length(exps) == 1
            % 1.) Remove ice down to the minimum extent for this glacier
            exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
            levelset(pos) = +1;

            % 2.) Add ice using the extent of this month's exp file
            filename = exps(1).name;
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('PROBLEM\n');
            return
         end
      end

      % Populate spclevelset
      [doy,frac] = date2doy(datenum(y,8,15));
      if isnan(md.levelset.spclevelset)
         md.levelset.spclevelset = [levelset; y+frac];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y+frac];
      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
if perform(org,'Terminus10yearly'),% {{{ STEP 12

   md=loadmodel(org,['Relaxation' num2str(relaxation_years) 'yr']);

   %Transfer resuts to model fields
   md=transientrestart(md);

   % Set parameters
   md.inversion.iscontrol = 0;
   md.timestepping.time_step = .02;

   % md.timestepping.final_time = 2015; md.settings.output_frequency = (1/md.timestepping.time_step)/5; % forward run to 2015
   % md.timestepping.final_time = 2100 - start_year;  md.settings.output_frequency = (1/md.timestepping.time_step)/2; % forward run to 2100
   md.timestepping.final_time = end_year; md.settings.output_frequency = (1/md.timestepping.time_step)/5;

	% Go through all terminus exp files for the region
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=1;

	disp('Setting up levelset time series');
   if false
   % 1.) Find minimum ice extent and remove ice from original mask to initialize the levelset %%{{{
   levelset_min_extent = md.mask.ice_levelset;
   for iglacier = 1:numel(glaciers)
      glacier_name = glaciers{iglacier};
      glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
      glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

      exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
      pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
      levelset_min_extent(pos) = +1;
   end
   %%}}}
   % 2.) Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   for y = start_year:end_year_termini
      levelset = levelset_min_extent;
      %fprintf(' Finding termini for year %4d\n', year);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%4d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusYearly retreat exps found for %s\n', glacier_long_name);
         elseif numel(exps) == 1
	      	filename = exps(1).name;
	      	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('Multiple TerminusYearly exp files found for %s\n', glacier_long_name);
         end
      end
      if isnan(md.levelset.spclevelset) 
         md.levelset.spclevelset = [levelset; y];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y];
      end
   end %%}}}
   end
   % Use yearly terminus exps to extend ice from the minimum extent to the extent of each year %%{{{
   levelset = md.mask.ice_levelset;
   for y = start_year:10:end_year_termini
      %fprintf(' Finding termini for year %4d\n', y);
      %fprintf(' Finding termini for month %2d\n', m);
      for iglacier = 1:numel(glaciers)
         glacier_name = glaciers{iglacier};
         glacier_long_name = glacier_name_translator(glacier_name, 'terminus_name', 'centerline_name');
         glacier_abbr = glacier_name_translator(glacier_name, 'terminus_name', 'terminus_abbr');

	      exps = dir(sprintf('Exp/%s_TerminusYearly_%d.exp', glacier_long_name, y));
         if isempty(exps)
            %fprintf('No TerminusMonthly retreat exps found for %s, year %4d, month %s\n', glacier_long_name, y, month_str);
         elseif length(exps) == 1
            % 1.) Remove ice down to the minimum extent for this glacier
            exp_min_extent = ['Exp/' glacier_long_name '_TerminusMinExtent.exp'];
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,exp_min_extent,2));
            levelset(pos) = +1;

            % 2.) Add ice using the extent of this month's exp file
            filename = exps(1).name;
            pos = find(ContourToNodes(md.mesh.x,md.mesh.y,['Exp/' filename],2));
            levelset(pos) = -1;
         else
            fprintf('PROBLEM\n');
            return
         end
      end

      % Populate spclevelset
      [doy,frac] = date2doy(datenum(y,8,15));
      if isnan(md.levelset.spclevelset)
         md.levelset.spclevelset = [levelset; y+frac];
      else
         md.levelset.spclevelset(:,end+1) = [levelset; y+frac];
      end
   end
   %%}}}

   % % Check for stranded ice nodes
   % [stranded_nodes clean_mask_or_levelset] = remove_stranded_ice_spclevelset(md,'spclevelset');
   % if sum(sum(stranded_nodes)) > 0
   %    fprintf(['\n\033[' '103;30' 'm   WARNING: Stranded ice nodes exist for some spclevelsets. \033[0m \n', ...
   %               '\033[' '103;30' 'm            Check stranded_nodes and/or clean_mask_or_levelset \033[0m \n\n']);
   %    s = input('Continue by removing ice from stranded nodes (y/n)?','s');
   %    if strcmpi(s(1),'y')
   %       md.levelset.spclevelset = clean_mask_or_levelset;
   %    else
   %       return
   %    end
   % end
      
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
   md.stressbalance.requested_outputs={'default'};

   % Go solve
   if contains(cluster.name, 'discover')
      cluster.time = 12*3600;
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

end %}}}
return
% Extra transient
if perform(org,'Transient200'),% {{{

end
% }}}

