% Setup
region = 'WestGrIS';
start_year = 1985;

% Mesh sizing
triangleresolution = 1000;

% Mesh
md = model()
md = triangle(md,['./Exp/' region '.exp'],triangleresolution);

savevars(['./Models/' region '_mesh'], 'md', md)


% 
