function CalcGint(inputfile,outputfile)
%function CalcGint(inputfile,outputfile)
%29-06-12 Sebastian Christensen - sebastian@chem.au.dk
%Known bugs/missing elements: The multiplicity of the reflections has not
%been tested for lower symmetry than orthorhombic.

%Make sure that not too few or too many reflections are removed from
%Fconstraints. Reflections only be either in Fconstraint or in
%Gconstraints.

%Please check your results and report errors and/or suggestions to changes
%and improvements

%Setting up G constraints for BayMEM
%Give name of the BayMEM file used as input and a name for the outputfile
% run like: CalcGint("input.BayMEM","output.BayMEM")
%Inputfile: Should contain all F-constraints. Reflections should be given like:
% gbegin
% ggroup
%    0   0   2
%    2   0   0
% ggroup
%    0   2   2
%    2   2   0
% ggroup
%    1   1   3
%    3   1   1
% ggroup
%    1   3   1
%....
% endg

% %Program will change this to:
% gbegin
% ggroup 143.2020001     0.0665475
%    0   0   2
%    2   0   0
% ggroup 52.6714917     0.0977953
%    0   2   2
%    2   2   0
% ggroup 51.4950546     0.0719168
%    1   1   3
%    3   1   1
% ggroup 34.5584944     0.1403766
%    1   3   1
%....
%endg
%F-constraints

%Latest changes
%27012015 - Added the d-spacing in commandprompt output for missing
%reflections.
%16052015 - search for hkl equivalent to the ones written in G-constraint.
%Not only that specific hkl
%16052015 - a Gconstraint is kept even is some hkls are missing. If all hkls are missing the g-constraint is removed%
do_adjust_filenames=1;

if nargin == 1
    outputfile=inputfile;
    disp('Overwriting inputfile')
end

%read BayMEM file
fid = fopen(inputfile);
if fid==-1; error('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid)
%
outfile_base=regexprep(outputfile,'\.[a-zA-Z]*(?!.*\.)','');

%read symmetry operators of the laue class to determine symmetry equivalent reflections
laueR=makelaue(inputfile);
numlaue=size(laueR,3);
%read the unitcell dimensions
aa=regexp(x, '(?-s)(?m)cell[\t ]+([0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+)','tokens');
cellstr=aa{1}{1};
cell=str2num(cellstr);
a=cell(1); b=cell(2); c=cell(3); alpha=cell(4); beta=cell(5); gamma=cell(6);
cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
V=a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5;
cstar=a*b*sind(gamma)/V;
%her er noget galt med A(3,3)
A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
Gstar=(A^2)^-1;
%get positions of keywords nesting the F- and G-constraints sections
gbegin=regexp(x, '(?-s)(?m)gbegin', 'start');
endg=regexp(x, '(?-s)(?m)endg', 'end');
fbegin=regexp(x, '(?-s)(?m)fbegin', 'start');
endf=regexp(x, '(?-s)(?m)endf', 'end');

if isempty(gbegin); error('Keyword gbegin not found') ; end
if isempty(endg); error('Keyword endg not found') ; end
if isempty(fbegin); error('Keyword fbegin not found') ; end
if isempty(endf); error('Keyword endf not found'); end

G=x(gbegin:endg);
Gnew=G; %for rewriting and adding intensities
F=x(fbegin:endf);
Fnew=F;
%get get index separating the different ggroups
T = regexp(G, '(?-s)(?m)ggroup', 'start');
numG=length(T); %count number of g-groups
if numG==0; error('no g groups found'); end
T=[T length(G)];
%read hkl values
ghkl=[];
bad_g=[];
for i=1:numG  %number of g-groups
    str=G(T(i):T(i+1));
    missing_reflections=0;
    t=regexp(str,'([0-9]+[\t ]+[0-9]+[\t ]+[0-9]+)', 'tokens');
    numhkl=length(t); HKL=zeros(numhkl,3); I=zeros(numhkl,3); m=zeros(numhkl,1);
    gd=[];
    ind_found=zeros(numhkl,1);
    for k=1:numhkl;
        hkl=str2num(t{k}{1});
        
        %convert hkl string to the correct format: "   0   8   1"
        hklstr_org=sprintf('%d%4d%4d', hkl);
        if ismember(hkl,ghkl,'rows')
            error(['Reflection - ' hklstr_org ' - is a dublicate'])
        end
        ghkl=[ghkl; hkl];
        %calculate d for each hkl to make sure they are close to each other
        d=(hkl*Gstar*hkl')^-0.5;
        gd=[gd; d];
        
        %calculate multiplicity of the reflection
        hkleq=zeros(numlaue,3);
        for j=1:numlaue
            hkleq(j,:)=hkl*laueR(:,:,j);
        end
        hkleq=unique(hkleq,'rows');
        m(k)=size(hkleq,1);
        
        %read structure factor and sigmas of the reflection
        reflection_found=0;
        for ii=1:m(k)
            hklstr=sprintf(' %d%4d%4d ', hkleq(ii,:));
            intsig=regexp(F,['(?-s)(?m)^[ ]*' hklstr '[\t ]+([0-9\.\-]+)[\t ]+([0-9\.\-]+)[\t ]+([0-9\.\-]+)'], 'tokens');
            if isempty(intsig) ==0
                reflection_found=1;
                break
            end
            
        end
        
        if reflection_found==1;
            I(k,:)=str2num(char(intsig{1}))';
            %remove the relevant F-constraint
            Fnew=regexprep(Fnew,['(?m)^[ ]*' hklstr '[\d \t\-\.]*[\r\n]+'],'');
            HKL(k,:)=hkl;
            ind_found(k)=1;
        else
            disp(['Reflection - ' hklstr_org ' - not found. (' num2str(d) ' Å)']);
            Gnew=regexprep(Gnew,['(?m)^[ ]*' hklstr_org '[ \r]*\n'],'');
            
        end
    end
    
    if max(gd)/min(gd) > 1.1; error('Large difference in dspacings. Are you sure these should be grouped?'); end
    
    ind_found=(ind_found==1);
    if all(ind_found==0)
        bad_g=[bad_g; i];
        Gint=0;
        Gsig=0;
    else
    Gint=sqrt(sum(m(ind_found).*(I(ind_found,1).^2+I(ind_found,2).^2))./sum(m(ind_found)));
    Gsig=1./Gint.*sqrt(sum((sqrt(I(ind_found,1).^2+I(ind_found,2).^2).*m(ind_found).*I(ind_found,3)).^2))./sum(m(ind_found));
    end
    if isnan(Gsig) || isinf(Gsig)
        Gsig=max(I(:,3));
        
    end
    %prepare string for insertion in baymem file
    Gstr=num2str([Gint Gsig],'%14.7f');
    Gnew=regexprep(Gnew,'ggroup[\w\d\. \-]*',['ggroup ' Gstr],i);

end

%remove empty G-constraints
%get beginning of group.
group_start=regexp(Gnew, '(?-s)(?m)ggroup', 'start');
group_end=[group_start(2:end)-1 regexp(Gnew, '(?-s)(?m)endg', 'start')-1];

bad_group_start=group_start(bad_g);
bad_group_end=group_end(bad_g);

number_of_bad_groups=length(bad_group_start);

rmstr=zeros(size(Gnew));
for j=1:number_of_bad_groups
    rmstr(bad_group_start(j):bad_group_end(j))=1;
end
Gnew(find(rmstr))='';

xnew=regexprep(x,G,Gnew);       %exchange the G-constraint section
xnew=regexprep(xnew,F,Fnew);    %exchange the F-constraint section

%Replace the names to conform with the filename
%title
if do_adjust_filenames
    xnew=regexprep(xnew,'(title\s+)\S+',['$1' outfile_base]);
    xnew=regexprep(xnew,'(outputfile\s+)\S+(\.m81)',['$1' outfile_base '$2']);
    xnew=regexprep(xnew,'(outputprior\s+)\S+(_prior\.m81)',['$1' outfile_base '$2']);
end
%outputprior pbte_105K_anh_prior.m81

fidout=fopen(outputfile,'w'); %make file for output
fprintf(fidout,'%s',xnew); %write updated BayMEM file
fclose(fidout); %close file
end

function [Rall Tall]=readsym(inputfile)
%read symmetry operators from baymemfile
%read BayMEM file
fid = fopen(inputfile);
x=fread(fid,'*char')';
fclose(fid);
symops=regexp(x, '( [\-]*x[1-6][\+\d\/]*[\t ]+[\-]*x[1-6][\+\d\/]*[\t ]+[\-]*x[1-6][\+\d\/]*)', 'tokens');
Rall=zeros(3,3,length(symops)); Tall=zeros(3,1,length(symops));
for k=1:length(symops);
    opstr=symops{k}{1};
    transtr=regexprep(opstr,'[\-]*x[\d]','0');
    T=str2num(transtr)';
    rotstr=regexprep(opstr,'[\-\+][\/\d]+',' ');
    rotstr=regexprep(rotstr,'[^-]x1[\+\-\/\d]*','1 0 0 ;');
    rotstr=regexprep(rotstr,'[\-]x1[\+\-\/\d]*','-1 0 0 ;');
    rotstr=regexprep(rotstr,'[^-]x2[\+\-\/\d]*','0 1 0 ;');
    rotstr=regexprep(rotstr,'[\-]x2[\+\-\/\d]*','0 -1 0 ;');
    rotstr=regexprep(rotstr,'[^-]x3[\+\-\/\d]*','0 0 1 ;');
    rotstr=regexprep(rotstr,'[\-]x3[\+\-\/\d]*','0 0 -1 ;');
    R=str2num(rotstr)';
    
    
    Rall(:,:,k)=R; Tall(:,:,k)=T;
end

end

function laueR=makelaue(inputfile)
%returns the rotation matrices of the laue class
R=readsym(inputfile); %read symmetry
numsym=size(R,3);
%compare matrices and add centrosymmetry
resR=[reshape(R,9,numsym)' ; -1*reshape(R,9,numsym)'];
redR=unique(resR,'rows');
%reform 3x3xnumsym matrix
laueR=reshape(redR',3,3,[]);
end