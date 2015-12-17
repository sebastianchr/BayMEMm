function sigmareject_BayMEM(inputfile,limit,outputfile)
%Keep only reflections with significance (F_norm /sigma) higher than the
%limit

if nargin == 2
    outputfile=inputfile;
    disp('WARNING: Output will be written to inputfile')
elseif nargin==3

else
    error('incorrect number of parameters')
end

[hkl A B sigma]=readFobs_BayMEM(inputfile);
F_norm=sqrt(A.^2+B.^2);

ind_keep=F_norm > limit*sigma;

writeFobs_BayMEM(inputfile,outputfile,hkl(ind_keep,:),A(ind_keep),B(ind_keep),sigma(ind_keep));

end