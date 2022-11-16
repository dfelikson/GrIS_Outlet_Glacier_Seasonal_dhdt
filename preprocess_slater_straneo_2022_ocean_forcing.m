load /Users/dfelikso/Downloads/donaldaslater-slater_2022_submelt-97d9702/data/twglaciers.mat

morlighem_ids = [13, 52, 53, 90];

for morlighem_id = morlighem_ids
	eval(sprintf('load Cheat_matfiles/glacier%04d.mat;', morlighem_id));
	
	idx = find([twglaciers(:).morlighem_number] == morlighem_id);

	%twglaciers(idx).runoff.RACMO.t
	%twglaciers(idx).runoff.RACMO.Q

	eval(sprintf('glacier%04d.EN4.t  = twglaciers(idx).ocean.EN4.t;', morlighem_id));
	eval(sprintf('glacier%04d.EN4.TF = twglaciers(idx).ocean.EN4.TF_GL;', morlighem_id));

	eval(sprintf('save Cheat_matfiles/glacier%04d.mat glacier%04d;', morlighem_id, morlighem_id));
end

