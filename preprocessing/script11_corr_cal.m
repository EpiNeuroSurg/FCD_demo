clear;clc;
addpath /home/mo/Desktop/surfstat
path = '/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a';

filesleft = SurfStatListDir(['/usr/local/freesurfer_de/subjects/fsaverage_sym/surf/lh.pial']);
lhavsurf = SurfStatAvSurf(filesleft);
[~,l_label,~] = read_annotation(['/usr/local/freesurfer_de/subjects/fsaverage_sym/label/lh.aparc.annot']);
lhmask(1,163842) = 0; 
a1 = find(l_label == 2647065); a2 = find(l_label == 660700); a3 = find(l_label == 9231540); a4 = find(l_label == 7874740); 
a5 = find(l_label == 3302560); a6 = find(l_label == 3988500); a7 = find(l_label == 14474380); a8 = find(l_label == 11146310); 
a9 = find(l_label == 13145750); A12 = union(a1,a2); A123 = union(A12,a3); A1234 = union(A123,a4); A12345 = union(A1234,a5);
A123456 = union(A12345,a6); A1234567 = union(A123456,a7); A12345678 = union(A1234567,a8);A123456789 = union(A12345678,a9);
lhmask(1,A123456789) = 1; lhmask = logical(lhmask);

filesright = SurfStatListDir(['/usr/local/freesurfer_de/subjects/fsaverage_sym/surf/rh.pial']);
rhavsurf = SurfStatAvSurf(filesright);
[~,r_label,~] = read_annotation(['/usr/local/freesurfer_de/subjects/fsaverage_sym/label/rh.aparc.annot']);
rhmask(1,163842) = 0;
b1 = find(r_label == 660700); b2 = find(r_label == 2647065); b3 = find(r_label == 3302560); b4 = find(r_label == 3988500); 
b5 = find(r_label == 7874740); b6 = find(r_label == 9231540); b7 = find(r_label == 11146310); b8 = find(r_label == 13145750); 
b9 = find(r_label == 14474380); B12 = union(b1,b2); B123 = union(B12,b3); B1234 = union(B123,b4); B12345 = union(B1234,b5); 
B123456 = union(B12345,b6); B1234567 = union(B123456,b7); B12345678 = union(B1234567,b8); B123456789 = union(B12345678,b9);
rhmask(1,B123456789) = 1;rhmask = logical(rhmask);

%% all patients 
load('/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/SEEG_patient.mat')
all1 = []; all2 = []; all3 = [];
entire1 = []; entire2 = []; entire3 = [];
for i = 1:14
    t1 = SurfStatReadData({[path,'/',name{i,1},'/xhemi/surf/',name{i,2},'.gradient_z_on_',name{i,2},'.sm.mgh']});
    t2 = SurfStatReadData({[path,'/',name{i,1},'/xhemi/surf/',name{i,2},'.gm_PET_0.5_z_on_',name{i,2},'.sm.mgh']});
    t3 = SurfStatReadData({[path,'/',name{i,1},'/mapping_TT.mgh']}); mask_all = eval([name{i,2},'mask']);
    tt3 = t3.*mask_all;position = find(tt3>0);stru = t1(position);meta = t2(position);elec = t3(position);z_elec = normalize(elec);
    all1 = [all1,stru]; all2 = [all2,meta]; all3 = [all3,z_elec];
    entire1 = [entire1,t1]; entire2 = [entire2,t2]; entire3 = [entire3,t3]; 
end

entire1 = normalize(entire1); entire2 = normalize(entire2);

figure
scatter(all1,all2,20,all3,'fill');colormap(jet);
colorbar;caxis([0,6])
title('Structure-Metabolism-Epileptogenicity Correlation','FontSize',12,'FontWeight','Bold')
xlabel('FLAIR Gradient','FontSize',11,'FontWeight','Bold');
ylabel('PET Value','FontSize',11,'FontWeight','Bold');
c = colorbar;
c.Label.String = 'Ictal High-frequency Power';c.Label.FontSize = 11;c.Label.FontWeight = 'Bold';

figure
scatter(entire1,entire2,20,entire3,'fill');colormap(jet);
colorbar;caxis([0,6])
title('Structure-Metabolism-Epileptogenicity Correlation','FontSize',12,'FontWeight','Bold')
xlabel('FLAIR Gradient','FontSize',11,'FontWeight','Bold');
ylabel('PET Value','FontSize',11,'FontWeight','Bold');
c = colorbar;
c.Label.String = 'Ictal High-frequency Power';c.Label.FontSize = 11;c.Label.FontWeight = 'Bold';

%% individual stru-meat-elec-surg 
i = 9;
t1 = SurfStatReadData({[path,'/',name{i,1},'/xhemi/surf/',name{i,2},'.gradient_z_on_',name{i,2},'.sm.mgh']});
t2 = SurfStatReadData({[path,'/',name{i,1},'/xhemi/surf/',name{i,2},'.gm_PET_0.5_z_on_',name{i,2},'.sm.mgh']});
t3 = SurfStatReadData({[path,'/',name{i,1},'/mapping_TT.mgh']}); 
mask = eval([name{i,2},'mask']);
tt3 = t3.*mask; position = find(tt3>0); 
stru = zeros(1,163842); stru(position) = t1(position); 
meta = zeros(1,163842); meta(position) = t2(position); 
elec = zeros(1,163842); elec(position) = t3(position); z_elec = normalize(elec);
surg = SurfStatReadData({[path,'/',name{i,1},'/mapping_surg.mgh']});

figure
scatter(stru,meta,20,elec,'fill');colormap(jet);colorbar
title('Structural-Metabolic-Epileptogenic Correlation','FontSize',13,'FontWeight','Bold')
xlabel('FLAIR Gradient','FontSize',11,'FontWeight','Bold');
ylabel('PET Value','FontSize',11,'FontWeight','Bold');
c = colorbar;
c.Label.String = 'Ictal High-frequency Power';c.Label.FontSize = 11;c.Label.FontWeight = 'Bold';

figure
t1 = normalize(t1); t2 = normalize(t2);
scatter(t1,t2,20,mask,'fill');colormap(jet);colorbar
title('Structural-Metabolic-Epileptogenic Correlation','FontSize',13,'FontWeight','Bold')
xlabel('FLAIR Gradient','FontSize',11,'FontWeight','Bold');
ylabel('PET Value','FontSize',11,'FontWeight','Bold');
c = colorbar;
c.Label.String = 'Ictal High-frequency Power';c.Label.FontSize = 11;c.Label.FontWeight = 'Bold';

figure
SurfStatView(t1.*mask, eval([name{i,2},'avsurf']), 'Cort Thick (mm), FreeSurfer data');
figure
SurfStatView(t2.*mask, eval([name{i,2},'avsurf']), 'Cort Thick (mm), FreeSurfer data');
figure
SurfStatView(z_elec.*mask, eval([name{i,2},'avsurf']), 'Cort Thick (mm), FreeSurfer data');

figure
SurfStatView(surg.*mask, eval([name{i,2},'avsurf']), 'Cort Thick (mm), FreeSurfer data');

surg_t1 = surg.*t1; surg_t2 = surg.*t2; surg_t3 = surg.*elec; 
tt1 = t1.*mask; tt2 = t2.*mask; tt3 = elec.*mask;
tt1 (tt1 == 0) = 100; tt2 (tt2 == 0) = 100; tt3(tt3 == 0) = 100;
surg_t1 (surg_t1 == 0) = 100; surg_t2 (surg_t2 == 0) = 100; surg_t3 (surg_t3 == 0) = 100; surg_t3 (surg_t3 < 0) = 100;

map = brewermap(3,'Set1'); 
figure
histogram(tt2,-3:0.01:3,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
hold on
histogram(surg_t2,-3:0.01:3,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
set(gca,'YLim',[0 200])
box off
t1_02 = find (tt1 > 0 & tt1 < 2); t2_02 = find (tt2 > -2 & tt2 < 20); overlay1 = intersect (t1_02, t2_02); 
surg_res = find (surg == 1);overlay2 = intersect (overlay1, surg_res); rate1 = (size(overlay2,2))./(size(overlay1,2));
overlay3 = intersect (position,surg_res); rate2 = (size(overlay3,2))./(size(position,2));
% [R,P] = corrcoef(stru,meta); 
% text(2,2,['R = ',num2str(R(1,2)),char(13,10)','P = ',num2str(P(1,2))])
% stru_position = find(0<t1&t1<2); stru1 = zeros(1,163842); stru1(stru_position) = t1(stru_position);
% meta_position = find(-2<t2&t2<0); meta1 = zeros(1,163842); meta1(meta_position) = t2(meta_position);
% elec_position = find(t3>0); elec1 = zeros(1,163842); elec1(elec_position) = t3(elec_position);
s = SurfStatReadSurf({...
    '/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/21wuzhisen/surf/lh.pial',...
    '/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/21wuzhisen/surf/rh.pial'});
t = SurfStatReadData({...
    '/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/21wuzhisen/surf/lh.complexity_z.sm.mgh',...
    '/media/mo/SAMSUNG1/2020/temporal_cases/recon/FCD3a/21wuzhisen/surf/rh.complexity_z.sm.mgh'});
figure
SurfStatView( t, s, 'Cort Thick (mm), FreeSurfer data' );