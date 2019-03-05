function plot_bed_and_termini(md, varargin)

   %%{{{ Parse input args
   p = inputParser;
   addRequired(p, 'md');
   addOptional(p, 'index', '');
   addOptional(p, 'step', 1);
   addOptional(p, 'fig', '');
   addOptional(p, 'legend', false);
   addOptional(p, 'title', false);
   parse(p, md, varargin{:});
   md = p.Results.md;
   index = p.Results.index;
   step = p.Results.step;
   fig = p.Results.fig;
   legendFlag = p.Results.legend;
   titleFlag = p.Results.title;
   %%}}}
   
   if isempty(fig)
      f = figure;
   else
      f = fig;
   end

   bedColormap;
   plotmodel(md, 'data', md.geometry.bed, 'caxis', [-400 400], 'colormap', cmap, 'figure', f)
   
   if isempty(index)
      index = 1:step:size(md.levelset.spclevelset,2);
   end

   nTermini = length(index);
   if nTermini < 4
      cmapTermini = spring;
   else
      cmapTermini = cbrewer('qual', 'Paired', nTermini);
   end
   cmapTermini = interp1(1:size(cmapTermini,1), cmapTermini, linspace(1, size(cmapTermini,1), nTermini));
   
   if legendFlag
      h = zeros(length(index),1);
      legend_str = {''};
   end
   
   k = 1;
   for i = index
      levelset = md.levelset.spclevelset(1:end-1,i);
      expcontourlevelzero(md,levelset,0,'TEMP.exp');
      expdisp('TEMP.exp', 'linestyle', cmapTermini(k,:), 'figure', f)
      time = md.levelset.spclevelset(end,i);
      dtime = time - md.timestepping.start_time;
      if legendFlag
         h(k) = plot(NaN,NaN,'-','color',cmapTermini(k,:));
         legend_str{k} = sprintf('%5.1f years', dtime);
      end
      k = k + 1;
   end
   
   if legendFlag
      legend(h, legend_str, 'Location', 'East');
   end

   if ~isempty(index) && titleFlag
      title(sprintf('%5.1f years', dtime))
   end

   delete('TEMP.exp'); 

end % main function

