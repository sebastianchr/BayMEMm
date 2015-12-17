%read fcalc from bmout-file
function [hkl0 A0 B0]=readBMout(inputfile)
% NOt finished (not really started)
% 
% inputfile='pbte_oc3_proffit_aff_po_aim0p06.BMout';

fid = fopen(inputfile);
if fid==-1; error('Inputfile not found'); return; end
x=fread(fid,'*char')';


fbegin=regexp(x, '(?-s)(?m)F-Constraints:', 'start');
endf=regexp(x, '(?m)F-Constraints:.*?\n\n', 'end');
fstr=x(fbegin:endf);

gbegin=regexp(x, '(?-s)(?m)G-Constraints:', 'start');
endg=regexp(x, '(?m)G-Constraints:.*?\n\n\n', 'end');
gstr=x(gbegin:endg);


%read ffor F-constraints
% trim fstr
newlines=regexp(fstr,'\n','start');
fstr=fstr(newlines(2):end);
%  #     h   k   l  A/Gobs B A/Gmem B DeltaF Sigma DeltaF/Sigma sinth/l  (Gw)

dat=sscanf(fstr,'%i:%i%i%i%f%f%f%f%f%f%f%f\n');
FdatF=reshape(dat,12,[])';

% read from G-constraints
groups=regexp(gstr,'(?-s)(?m)group nr', 'start');

FdatG=[];
Gdat=[];
for i=1:length(groups)-1
    groupstr=gstr(groups(i):groups(i+1));
    groupnewline=regexp(groupstr,'\n','start');
%        #   group Gobs GMEM deltaG  sigma  deltaG
    gcon_str=groupstr(1:groupnewline);
    gdat=sscanf(gcon_str,'group nr. %i%f%f%f%f%f\n');
    
    Gdat=[Gdat; gdat'];
    groupfstr=groupstr(groupnewline:end);
%      #  h k l AMEM BMEM sigma 
    fdat=sscanf(groupfstr,'%i:%i%i%i%f%f%f\n');
    fdat=reshape(fdat,7,[])';
    
    FdatG=[FdatG; fdat];
end

    groupstr=gstr(groups(end):end);
    groupnewline=regexp(groupstr,'\n','start');
    gcon_str=groupstr(1:groupnewline);
    gdat=sscanf(gcon_str,'group nr. %i%f%f%f%f%f\n');
    
    Gdat=[Gdat; gdat'];
    groupfstr=groupstr(groupnewline:end);
    fdat=sscanf(groupfstr,'%i:%i%i%i%f%f%f\n');
    fdat=reshape(fdat,7,[])';
    
    FdatG=[FdatG; fdat];
noFF=FdatF(:,1);
noFG=FdatG(:,1);

hkl0=[FdatF(:,2:4); FdatG(:,2:4)];
A0=[FdatF(:,5); FdatG(:,5)];
B0=[FdatF(:,6); FdatG(:,6)];

end

    
    
    



