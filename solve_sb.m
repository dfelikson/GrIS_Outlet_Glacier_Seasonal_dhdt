function md = solve_sb(md)
   md.cluster = generic;
   md.cluster.interactive = 1;
   md.cluster.np = 4;

   %md.cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 28, ...
   %   'login', 'dfelikso', ...
   %   'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
   %   'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
   %   'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution', ...
   %   'valgrind', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/install/bin/valgrind', ...
   %   'valgrindlib', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/install/lib/valgrind/libmpiwrap-amd64-linux.so', ...
   %   'valgrindsup', '/home/dfelikso/Software/ISSM/trunk-jpl/externalpackages/valgrind/issm.supp');
   %md.cluster.interactive = 0; %1;
   %md.settings.waitonlock = nan; %90;

   %md.verbose.convergence = true;
   %md.verbose.control = false;
   %md.verbose.qmu = false;

   md.settings.waitonlock = Inf;

   md = solve(md, 'sb');
end

