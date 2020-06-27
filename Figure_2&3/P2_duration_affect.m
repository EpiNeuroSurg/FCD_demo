clear;clc;
addpath /home/mo/Desktop/surfstat
D = readtable('/media/mo/SAMSUNG1/2020/demographic.txt'); 
name = table2array(D(1:64,1));
Group = table2array(D(1:64,2)); group = term(Group);
Sex = table2array(D(1:64,3)); sex = term(Sex); 
Age = table2array(D(1:64,4)); age = term(Age);
Side = table2array(D(1:64,7)); lside = find(Side==1); rside = find(Side==0); hcside = find(Side==2);
Duration = table2array(D(1:64,5)); duration = term(Duration);
Duration2 = Duration.*Duration; duration2 = term(Duration2);
Duration3 = Duration2.*Duration; duration3 = term(Duration3);
Frequency = table2array(D(1:64,6));
Severity = Duration.*Frequency; severity = term(Severity);
Severity2 = Severity.*Severity; 
Severity3 = Severity2.*Severity;

%thickness w-g.pct complexity gradient gm_FLAIR_0.5 myeline_2 gm_PET_0.5 
feature = 'gradient';
ipsi = cell(64,1); 
for k = 1:length(lside)
    j = lside(k,1);  
    ipsi{j,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/lh.',char(feature),'_z_on_lh.sm.mgh'];
end
for k = 1:length(rside)
    j = rside(k,1);
    ipsi{j,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/rh.',char(feature),'_z_on_lh.sm.mgh'];
end

Y = SurfStatReadData(ipsi);
[u,k,subj] = unique([floor(Age) group.pt],'rows');
Subj = term(var2fac(subj));
M = 1+severity+age+sex+random(Subj)+I;
filesleft = SurfStatListDir('/usr/local/freesurfer_de/subjects/fsaverage_sym/surf/lh.pial');
avsurf = SurfStatAvSurf(filesleft);
slm = SurfStatLinMod(Y,M,avsurf);
slm = SurfStatT(slm, Severity);

[a,l_label,b] = read_annotation('/usr/local/freesurfer_de/subjects/fsaverage_sym/label/lh.aparc.annot');
mask(1,163842) = 0; 
a1 = find(l_label == 2647065); a2 = find(l_label == 660700); a3 = find(l_label == 9231540); a4 = find(l_label == 7874740); 
a5 = find(l_label == 3302560); a6 = find(l_label == 3988500); a7 = find(l_label == 14474380); a8 = find(l_label == 11146310); 
a9 = find(l_label == 13145750); A12 = union(a1,a2);A123 = union(A12,a3);A1234 = union(A123,a4);
A12345 = union(A1234,a5);A123456 = union(A12345,a6);A1234567 = union(A123456,a7);
A12345678 = union(A1234567,a8);A123456789 = union(A12345678,a9);
mask(1,A123456789) = 1;
mask = logical(mask);

p = 1-tcdf(slm.t,slm.df);
pp = p<0.05;
slm.d = (slm.t*2)/(sqrt(slm.df+1));
SurfStatView(slm.d.*mask.*pp,avsurf,'d value')
% figure
% [pval, peak, clus] = SurfStatP(slm,mask);
% pval.PP = pval.P<0.05;
% aaaa = find(pval.PP==1);
% SurfStatView(slm.t.*mask.*pval.PP, avsurf, 'T value' );
% SurfStatColormap( 'jet' );
figure
all_d = find(slm.d == slm.d.*pp);
meand = Y(:,all_d); 
meand = mean(meand,2);
SurfStatPlot(duration,meand,1,age,'LineWidth',2, 'MarkerSize',12 )
