function C = basal_friction_initial_guess(md)
   
   disp('   Initialize basal friction using driving stress');
   disp('      -- Compute surface slopes and use 10 L2 projections');
   [sx,sy,s]=slope(md); sslope=averaging(md,s,10);
   disp('      -- Process surface velocity data');
   vel = md.inversion.vel_obs;
   flags=(vel==0); pos1=find(flags); pos2=find(~flags);
   vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
   %velmax = max(vel);
   %vel(vel==0 & md.mask.ice_levelset<0) = velmax;
   disp('      -- Filling in missing ice velocity with MEaSUREs mosaic');
   [velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
   vel = sqrt( velx.^2 + vely.^2 );
   idx = md.mask.ice_levelset < 0 & isnan(vel);
   vel(idx) = sqrt( velx(idx).^2 + vely(idx).^2 );
   vel=max(vel,0.1);
   disp('      -- Calculate effective pressure');
   Neff = md.constants.g .* (md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base);
   Neff(find(Neff<=0))=1;
   % -- NOTE --
   Neff=max(Neff,5e4);
   pos1 = find(Neff==5e4 & vel <100 & md.mask.ice_levelset<0);
   pos2 = find(Neff >5e4 & vel>=100 & md.mask.ice_levelset<0);
   vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1),'nearest');
   % -- NOTE --
   %Neff = effectivepressure(md);
   %vel = basal_vel_from_surface_vel_obs(md);

   disp('      -- Deduce friction coefficient');
   q = averaging(md, md.friction.q, 1);
   p = averaging(md, md.friction.p, 1);
   r = q./p;
   s = 1./p;
   C=sqrt(md.materials.rho_ice*md.constants.g*md.geometry.thickness.*(sslope)./(Neff.^r.*(vel/md.constants.yts).^s));
   C = min(C, 200);

   %Cmax = sqrt(1./vel.^(s-1)) .* 200;
   %C=min(C,Cmax);
   
   %switch mean(p)
   %   case 1
   %      C=min(C,200);
   %   case 2
   %      C=min(C,1e5);
   %   case 5
   %      C=min(C,1e13);
   %   case 20
   %      C=min(C,1e55);
   %   otherwise
   %      disp('sorry!')
   %      C = nan;
   %      return
   %end

   disp('      -- Extrapolate on ice free and floating ice regions');
   flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
   C(pos1) = 1;
   pos=find(isnan(C));
   C(pos)  = 1;

end % main function

