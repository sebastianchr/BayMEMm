function [F_out G_out]=compare_BayMEM(filename1, filename2)
% function compare_BayMEM(filename1, filename2)
% Compare F_obs and G_obs from 2 BayMEM files by plotting the ratio
% F_obs1/F_obs2, and G_obs1/G_obs2 including propagated errors.
%
% Input: filename1, filename2: The full filenames of the 2 BayMEM files you
%        wish to compare.
%
% Dependencies: getcell_BayMEM.m, readFobs_BayMEM.m, readGobs_BayMEM.m

%%
%Compare F_obs

[hkl1_0 A1_0 B1_0 sigF1_0]=readFobs_BayMEM(filename1);
[hkl2_0 A2_0 B2_0 sigF2_0]=readFobs_BayMEM(filename2);

[is_common,ind] = ismember(hkl1_0,hkl2_0,'rows');
%iscommon: entry i is 1 if hkl1(i,:) is a row in hkl2
%ind(i) is location in hkl2 of row number i in hkl1

hkl1=hkl1_0(is_common,:);
A1=A1_0(is_common); B1=B1_0(is_common); sigF1=sigF1_0(is_common);
hkl2=hkl2_0(ind(is_common),:);
A2=A2_0(ind(is_common)); B2=B2_0(ind(is_common)); sigF2=sigF2_0(ind(is_common));

F1=sqrt(A1.^2+B1.^2);
F2=sqrt(A2.^2+B2.^2);

%%
%list of hkl values not in common
%%

A=getcell_BayMEM(filename1);
Gstar=(A^2)^-1;
d_inv=dot((hkl1*Gstar)',hkl1')'.^(1/2);
    adata0=[];
    adata0.hkl=hkl1;
    adata0.dinv=d_inv;
figure('OuterPosition',[100,100,1000,600])
subplot(1,2,1)
y=F1./F2;
sig_y=sqrt(sigF1.^2./F2.^2+(F1./F2.^2).^2.*sigF2.^2);
errorbarplus(d_inv,y,sig_y,adata0,'.')
xlim([min(d_inv) max(d_inv)])
line([min(d_inv) max(d_inv)],[1 1],'color','r','linestyle','--')
xlabel('d^{-1} (Å^{-1})')
ylabel('F_{obs,1} / F_{obs,2}')
% title(outputname,'Interpreter','none')
F_out=[hkl1,d_inv,y,sig_y];

%%
%compare G constraints

[G1 sigG1 hklG1]=readGobs_BayMEM(filename1);
[G2 sigG2 hklG2]=readGobs_BayMEM(filename2);

if isempty(G1) | isempty(G2)
disp('No Gconstraints were found in one of the files')
    G_out=[];
else
G1c=[]; sigG1c=[];
G2c=[]; sigG2c=[];
ind_commonG=[];
for i=1:length(G1)
    is_match=0;
    for j=1:length(G2)
        try
            if all(hklG1{i}==hklG2{j})
                is_match=1;
                break
            end
        catch Error
            %disp('caught error');
        end
    end
    if is_match
        G1c(end+1,:)=G1(i);
        sigG1c(end+1,:)=sigG1(i);
        G2c(end+1,:)=G2(j);
        sigG2c(end+1,:)=sigG2(j);
        ind_commonG(end+1,:)=[i j];
    end
end

subplot(1,2,2)
y=G1c./G2c;
sig_y=sqrt(sigG1c.^2./G2c.^2+(G1c./G2c.^2).^2.*sigG2c.^2);
adata1=[];
errorbarplus(ind_commonG(:,1),y,sig_y,adata1,'.')
xlim([min(ind_commonG(:,1)) max(ind_commonG(:,1))])
line([min(ind_commonG(:,1)) max(ind_commonG(:,1))],[1 1],'color','r','linestyle','--')
xlabel('Index of Gcon in file no. 1')
ylabel('G_{obs,1} / G_{obs,2}')
% title(outputname,'Interpreter','none')
G_out=[ind_commonG(:,1),y,sig_y];
end
end