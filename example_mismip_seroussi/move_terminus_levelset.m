function signeddistance_out = move_terminus_levelset(md, icelevelset_in, magnitude, direction)

   % NOTE
   %  direction == -1: less ice (make signed distance values more positive)
   %            == +1: more ice (make signed distance values more negative)

   signeddistance_out = nan * ones(size(icelevelset_in));
   signeddistance_out(end,:) = icelevelset_in(end,:);

   % Create signed distance field
   if size(icelevelset_in,2)>1,
      disp('Converting levelsets to signed distance fields');
      for i=1:size(icelevelset_in,2)
         levelset = icelevelset_in(1:end-1,i);
         pos      = find(levelset<0);

         if exist('TEMP.exp','file'), delete('TEMP.exp'); end
         expcontourlevelzero(md, levelset, 0, 'TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
         %delete('TEMP.exp');
         levelset(pos) = -levelset(pos);
         signeddistance_out(1:end-1,i) = levelset;
      end
   end

   % (possibly) Check that there is only one zero-level contour

   % Alter the signed distance field based on magnitude and direction
   switch sign(direction)
      case -1
         disp('Removing ice from the levelset(s)');
      case +1
         disp('Adding ice to the levelset(s)');
   end
   signeddistance_out(1:end-1,:) = signeddistance_out(1:end-1,:) - sign(direction) * magnitude;

end % main function

