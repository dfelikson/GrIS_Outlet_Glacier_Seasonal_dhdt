
steps = [1];
cluster = ''; % empty for localhost
%cluster = 'oibserve';
%cluster = 'discover';

%% Setup
region = 'WestGrIS';

switch cluster %%{{{
   case ''
      md.cluster = generic('name', oshostname(), 'np', 2);
   case 'oibserve'
      md.cluster = generic('name', 'oibserve', 'np', 4, ...
         'login', 'dfelikso', ...
         'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
         'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
         'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
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

%%}}}

