function makeprior_BayMEM(baymem_file,m40_file)
filebase=regexprep(baymem_file,'\.BayMEM','');

if nargin<2
    m40_file=[filebase '.m40'];
    disp(['m40-filename not given. Will look for: ' m40_file])
end

outputfile=[filebase '.prior'];
priorm81_file=[filebase '_prior.m81'];

%read BayMEM file
fid = fopen(baymem_file);
if fid==-1; error('BayMEM-file not found'); return; end


x=fread(fid,'*char')';
fclose(fid);

m40=readm40new(m40_file,5);
atomstr='';
for i=1:length(m40.name);
    name=m40.name{i};
    occ=m40.occupancy(i);
    type=['n' regexprep(name,'[0-9]+.*','')];
    position=m40.position(i,:);
    adp=m40.adp{i}(1:6);
    atomstr=[atomstr sprintf('%s %s %9.6f%9.6f%9.6f%9.6f%9.6f%9.6f%9.6f%9.6f%9.6f%9.6f\n',name,type,occ,position,adp)];
end
prior_str=['prioratoms\n' atomstr '\nendprioratoms\noutputprior ' priorm81_file ];
x=regexprep(x, 'endsymmetry',['endsymmetry\n' prior_str]);

fid=fopen(outputfile,'w');
fprintf(fid,'%s',x);
fclose(fid);




end