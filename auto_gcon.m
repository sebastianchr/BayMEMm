function auto_gcon(inputfile, outputfile, selection_criterion)
%function auto_gcon(inputfile, outputfile, selection_criterion)
% Description: Generate a list of G-constraints based on the position and
% width of adjacent peaks. The information is taken from the jana prf-file.
% Inputparameters:
% inputfile: Name of the prf file. Include fileextension
% outputfile: Name of 
% selection_criterion: 2 peaks are bundled if tth1-tth2 <
% (fwhm1+fwhm2)/2*selection_criterion

[~,~,prfdat]=read_prf(inputfile);

exp_tth_min=min(prfdat.twotheta_pattern);
exp_tth_max=max(prfdat.twotheta_pattern);
exp_range= prfdat.twotheta_hkl > exp_tth_min & prfdat.twotheta_hkl< exp_tth_max;



hkl=prfdat.hkl(exp_range,:);
tth=prfdat.twotheta_hkl(exp_range);
width=prfdat.peakwidth_hkl(exp_range);


tth_diff=tth(2:end)-tth(1:end-1);
width_sum=(width(2:end)+width(1:end-1))/2;


inrange=tth_diff<=width_sum*selection_criterion;

gcon_all={};
gcon=[];
for i=1:length(inrange)
   hkli=hkl(i,:);
   gcon=[gcon; hkli];
   if inrange(i)==0
       if size(gcon,1)>1
          gcon_all{end+1}=gcon;
          gcon=[];
       else
           gcon=[];
       end
   end
end


% write suggested gconstraints to file.
fid=fopen(outputfile,'w');
fprintf(fid,'%s\n','gbegin');
for g=gcon_all
   gi=g{1};
   fprintf(fid,'%s\n','ggroup ');
   fprintf(fid,'%4d%4d%4d\n',gi');
    
end
fprintf(fid,'%s','endg');
fclose(fid);
disp('')
disp(['Suggested G-contraints have been written to: ' outputfile])

end
