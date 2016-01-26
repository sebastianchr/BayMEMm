function A=getcell_BayMEM(inputfile)

%read BayMEM file
fid = fopen(inputfile);
if fid==-1; error('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid);

aa=regexp(x, '(?-s)(?m)cell[\t ]+([0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+)','tokens');
cellstr=aa{1}{1};
cell=str2num(cellstr);

%from giacovazzo s 75
A=cell2A(cell);
%first row in A is the first lattice vector and so forth.
%Be carefull som may use the transpose definition
end

function A=cell2A(cell)
%return vector representation of unit cell. 
a=cell(1); b=cell(2); c=cell(3); alpha=cell(4); beta=cell(5); gamma=cell(6);
cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
V=a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5;
cstar=a*b*sind(gamma)/V;
%her er noget galt med A(3,3). Har sammenlignet med giacovazzo og det ser
%ud til at være ok
A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
%first row in A is the first lattice vector and so forth.
%Be carefull som may use the transpose definition
end