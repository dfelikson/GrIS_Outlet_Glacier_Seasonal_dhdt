function N = effectivepressure(md)

   N = nan * ones(md.mesh.numberofvertices,1);
   N = md.materials.rho_ice .* md.constants.g .* md.geometry.thickness + md.materials.rho_water .* md.constants.g .* md.geometry.base;

   N(find(N<=0))=1;
   N=max(N,5e4);

   % Where there is no ice, set effective pressure to nan
   %N(md.mask.ice_levelset>0) = nan;

   % md.friction.coupling
   %  == 0: uniform sheet (negative pressure ok, default)
   %  == 1: ice pressure only
   %  == 2: water pressure assuming uniform sheet (no negative pressure)
   %  == 3: use provided effective_pressure
   %  == 4: use coupled model (not implemented yet)

   % the code below is for md.friction.coupling == 1
   %N(md.geometry.bed< 0) = md.materials.rho_ice * md.constants.g * md.geometry.thickness(md.geometry.bed<0 ) + md.materials.rho_water * md.constants.g * md.geometry.bed(md.geometry.bed<0);
   %N(md.geometry.bed>=0) = md.materials.rho_ice * md.constants.g * md.geometry.thickness(md.geometry.bed>=0);

   % Make sure effective pressure is positive
   %N(N<0) = 0;

end % main function

