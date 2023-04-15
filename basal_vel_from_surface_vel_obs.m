function vb = basal_vel_from_surface_vel_obs(md)

   [sx,sy,s]=slope(md); sslope=averaging(md,s,10);
   vel = md.inversion.vel_obs;
   flags=(vel==0); pos1=find(flags); pos2=find(~flags);
   vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
   %velmax = max(vel);
   %vel(vel==0 & md.mask.ice_levelset<0) = velmax;
   [velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
   vel = sqrt( velx.^2 + vely.^2 );
   idx = md.mask.ice_levelset < 0 & isnan(vel);
   vel(idx) = sqrt( velx(idx).^2 + vely(idx).^2 );
   vel=max(vel,0.1);

   Neff = effectivepressure(md);

   pos1 = find(Neff==5e4 & vel <100 & md.mask.ice_levelset<0);
   pos2 = find(Neff >5e4 & vel>=100 & md.mask.ice_levelset<0);
   vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1),'nearest');

   vb = vel;

end % main function

