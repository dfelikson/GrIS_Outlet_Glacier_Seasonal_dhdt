steps=[88];
modelnum =1;

switch(modelnum),
	case 1,  modelname = '2km_viscous';
	case 2,  modelname = '1km_viscous';
	case 3,  modelname = '500m_viscous';
	case 4,  modelname = '250m_viscous';
	case 5,  modelname = '2km_coulomb';
	case 6,  modelname = '1km_coulomb';
	case 7,  modelname = '500m_coulomb';
	case 8,  modelname = '250m_coulomb';
	case 9,  modelname = '125m_viscous';
	case 10, modelname = '125m_coulomb';
end
%Hard coded parameters
loadonly = 1;
cluster=generic('name','murdo','np',20);

%Run Steps
org=organizer('repository',['./Models_' modelname ],'prefix',['MISMIP_' modelname '_'],'steps',steps);

%Initialization
if perform(org,'Mesh_generation'),% {{{1

	if modelnum==1 | modelnum==5,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',2000,'splitcorners',1);
	elseif modelnum==2 | modelnum==6,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',1000,'splitcorners',1);
	elseif modelnum==3 | modelnum==7,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',500,'splitcorners',1);
	elseif modelnum==4 | modelnum==8,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',250,'splitcorners',1);
	elseif modelnum==9 | modelnum==10,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',1000,'splitcorners',1,'maxnbv',4*10^6);
		hvertices=NaN*ones(md.mesh.numberofvertices,1);
		pos=find(md.mesh.x>200000 & md.mesh.x<800000);
		hvertices(pos)=500;
		pos=find(md.mesh.x>300000 & md.mesh.x<600000);
		hvertices(pos)=125;
		md=bamg(md,'hVertices',hvertices,'hmax',1000);
	elseif modelnum==11,
		md=bamg(model,'domain','./Exp_Par/Domain.exp','hmax',500,'splitcorners',1,'maxnbv',4*10^6);
		hvertices=NaN*ones(md.mesh.numberofvertices,1);
		pos=find(md.mesh.x>300000 & md.mesh.x<600000);
		hvertices(pos)=125;
		md=bamg(md,'hVertices',hvertices,'hmax',500);
	else
		error('model not supported yet');
	end
	md.miscellaneous.name=['MISMIP_' modelname ''];

	savemodel(org,md);
end % }}}
if perform(org,'Parameterization'),% {{{1

	md=loadmodel(org,'Mesh_generation');

	md=setmask(md,'','');
	md=parameterize(md,'./Exp_Par/Mismip.par');

	savemodel(org,md);
end% }}}
if perform(org,'Transient_Steadystate'),% {{{1

	md=loadmodel(org,'Parameterization');

	md=setflowequation(md,'SSA','all');

	if(modelnum==5 | modelnum==6 | modelnum==7 | modelnum==8 | modelnum==10),
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
	
	md.timestepping.time_step=1;
	md.timestepping.final_time=200000;
	md.settings.output_frequency=500;
	md.stressbalance.maxiter=10;
	md.stressbalance.abstol=NaN;
	md.stressbalance.restol=1;
	md.verbose=verbose('convergence',false,'solution',true);
	md.cluster=cluster;
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
