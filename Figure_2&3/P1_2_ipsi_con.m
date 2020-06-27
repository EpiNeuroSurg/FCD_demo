clear;clc;
addpath /home/mo/Desktop/surfstat

D = readtable('/media/mo/SAMSUNG1/2020/demographic.txt'); 
Group = table2array(D(1:128,2)); group = term(Group);
Side = table2array(D(1:64,7)); lside = find(Side==1); rside = find(Side==0);
name = table2array(D(1:64,1));

%thickness w-g.pct complexity gradient gm_FLAIR_0.5 myeline_2 gm_PET_0.5 
feature = 'gm_PET_0.5';
ipsi = cell(128,1); 
for k = 1:length(lside)
    j = lside(k,1);  
    ipsi{j,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/lh.',char(feature),'_z_on_lh.sm.mgh'];
end
for k = 1:length(rside)
    j = rside(k,1);
    ipsi{j,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/rh.',char(feature),'_z_on_lh.sm.mgh'];
end
for k = 1:length(lside)
    j = lside(k,1);  
    ipsi{j+64,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/rh.',char(feature),'_z_on_lh.sm.mgh'];
end
for k = 1:length(rside)
    j = rside(k,1);
    ipsi{j+64,1} = ['/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/',char(name(j,1)),'/xhemi/surf/lh.',char(feature),'_z_on_lh.sm.mgh'];
end

Y = SurfStatReadData(ipsi);
M = 1+group;
filesleft = SurfStatListDir('/usr/local/freesurfer_de/subjects/fsaverage_sym/surf/lh.pial');
avsurf = SurfStatAvSurf(filesleft);
slm = SurfStatLinMod(Y,M,avsurf);
slm = SurfStatT(slm, group.pt-group.hc);

[a,l_label,b] = read_annotation('/usr/local/freesurfer_de/subjects/fsaverage_sym/label/lh.aparc.annot');
mask(1,163842) = 0; single_mask = mask;
a1 = find(l_label == 2647065); a2 = find(l_label == 660700); a3 = find(l_label == 9231540); a4 = find(l_label == 7874740); 
a5 = find(l_label == 3302560); a6 = find(l_label == 3988500); a7 = find(l_label == 14474380); a8 = find(l_label == 11146310); 
a9 = find(l_label == 13145750); A12 = union(a1,a2);A123 = union(A12,a3);A1234 = union(A123,a4);
A12345 = union(A1234,a5);A123456 = union(A12345,a6);A1234567 = union(A123456,a7);
A12345678 = union(A1234567,a8);A123456789 = union(A12345678,a9);
mask(1,A123456789) = 1; single_mask(1,a8) = 1;
mask = logical(mask); single_mask = logical(single_mask);

% SurfStatView(single_mask, avsurf, 'T value' );
% cohen's d value calculation
slm.d = (slm.t*2)/(sqrt(slm.df+1));
[pval, peak, clus] = SurfStatP(slm,mask);
pval.PP = pval.P<0.05;
all_dvalue = slm.d.*single_mask.*pval.PP;mean_dvalue = mean(all_dvalue(all_dvalue>0));
SurfStatView(slm.d.*mask.*pval.PP, avsurf, 'T value' );
SurfStatColormap( 'jet' );
