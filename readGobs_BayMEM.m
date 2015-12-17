function [gobs sigma hkl]=readGobs_BayMEM(inputfile)

fid = fopen(inputfile);
if fid==-1; disp('inputfile not found'); return; end
x=fread(fid,'*char')';

gbegin=regexp(x, '(?-s)(?m)gbegin', 'start');
endg=regexp(x, '(?-s)(?m)endg', 'end');

g=x(gbegin:endg);
intsig=regexp(g,['(?-s)ggroup[\t ]+([0-9\.\-]+)[\t ]+([0-9\.\-]+)[\n\r ]+([[0-9]+[ ]+[0-9]+[ ]+[0-9]+[ \n\r]*]+)'], 'tokens');
n=length(intsig);

gobs=zeros(n,1); sigma=zeros(n,1);

hkl={};
for i =1:n
   gobs(i)=str2double(intsig{i}{1}); sigma(i)=str2double(intsig{i}{2});
   hkl{i}=str2num(intsig{i}{3});
end

end