function clean_BayMEM(inputfile,sintheta_lambda_observation_limit,outputfile)
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



[hkl A B sigma]=readFobs_BayMEM(inputfile);
Fnorm=sqrt(A.^2+B.^2);

Acell=getcell_BayMEM(inputfile);
M=2*pi*(Acell^-1);
qnorm_aa=sqrt(sum((hkl*M').^2,2));
if length(sintheta_lambda_observation_limit)==1
    q_limit=4*pi*sintheta_lambda_observation_limit;
    lower_limit=0;
    upper_limit=q_limit;
elseif length(sintheta_lambda_observation_limit)>=2
    
    q_limit=4*pi*sintheta_lambda_observation_limit(1:2);
    q_limit=sort(q_limit);
    lower_limit=q_limit(1);
    upper_limit=q_limit(2);
    if length(sintheta_lambda_observation_limit)==3
        observation_limit=sintheta_lambda_observation_limit(3);
        
        Int=Fnorm.^2;
        sigma_Int=2*Fnorm.*sigma;
        
        observed=Int>sigma_Int*observation_limit;
    else
        observed=ones(size(Fnorm));
    end
    
end

within_limits=qnorm_aa==0 | (qnorm_aa > lower_limit & qnorm_aa <upper_limit);

keep=within_limits & observed;
writeFobs_BayMEM(inputfile,outputfile,hkl(keep,:),A(keep),B(keep),sigma(keep));

end