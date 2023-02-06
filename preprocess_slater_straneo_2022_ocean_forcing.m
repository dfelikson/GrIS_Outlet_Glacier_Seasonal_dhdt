load /Users/dfelikso/Research/Data/donaldaslater-slater_2022_submelt-97d9702/data/twglaciers.mat

morlighem_ids = [13, 52, 53, 90];

for morlighem_id = morlighem_ids
	eval(sprintf('load Cheat_matfiles/glacier%04d.mat;', morlighem_id));
	
	idx = find([twglaciers(:).morlighem_number] == morlighem_id);
	eval(sprintf('glacier = twglaciers(idx);', morlighem_id));

	%twglaciers(idx).runoff.RACMO.t
	%twglaciers(idx).runoff.RACMO.Q

	% Plot the various ocean thermal forcings
	h = figure; hold on;
	plot(glacier.ocean.ORAS5.t, glacier.ocean.ORAS5.TF_GL, '.-');
	plot(glacier.ocean.EN4.t, glacier.ocean.EN4.TF_GL, '.-');
	plot(glacier.ocean.ASTE.t, glacier.ocean.ASTE.TF_GL, '.-');
	plot(glacier.ocean.CHORE.t, glacier.ocean.CHORE.TF_GL, '.-');

	title(glacier.name)
	xlabel('year')
	ylabel('ocean TF (deg. C)')
	legend('ORAS5', 'EN4', 'ASTE', 'CHORE')
	set(h,'color','w');
	continue

	eval(sprintf('glacier%04d.EN4.t  = twglaciers(idx).ocean.EN4.t;', morlighem_id));
	eval(sprintf('glacier%04d.EN4.TF = twglaciers(idx).ocean.EN4.TF_GL;', morlighem_id));

	eval(sprintf('save Cheat_matfiles/glacier%04d.mat glacier%04d;', morlighem_id, morlighem_id));
end

