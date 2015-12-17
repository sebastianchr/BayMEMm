function cutoff_BayMEM(inputfile,sintheta_lambda,outputfile)
%function cutoff_BayMEM(inputfile,sintheta_lambda,outputfile)
%if sintheta_lambda contains only one value it is interpreted at the upper
%value. If it contains 2 is is a lower and upper value
%Dependencies: readFobs_BayMEM, writeFobs_BayMEM, getcell_BayMEM
if nargin == 2
    outputfile=inputfile;
    disp('WARNING: Output will be written to inputfile')
elseif nargin==3

else
    error('incorrect number of parameters')
end

q_limit=4*pi*sintheta_lambda;

[hkl A B sigma]=readFobs_BayMEM(inputfile);
M=getcell_BayMEM(inputfile);
Mrl=2*pi*(M^-1)';
qnorm_aa=sqrt(sum((hkl*Mrl).^2,2));

if length(q_limit)==1
    lower_limit=0;
    upper_limit=q_limit;
else
   q_limit=sort(q_limit);
   lower_limit=q_limit(1);
   upper_limit=q_limit(2);
end

within_limits=qnorm_aa==0 | (qnorm_aa > lower_limit & qnorm_aa <upper_limit);

writeFobs_BayMEM(inputfile,outputfile,hkl(within_limits,:),A(within_limits),B(within_limits),sigma(within_limits));

end