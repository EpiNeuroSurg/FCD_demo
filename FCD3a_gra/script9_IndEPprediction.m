clear;clc
Neuropath = '/media/mo/SAMSUNG1/FCD3a_gra/recon/fcd3a';
SEEGpath = '/media/mo/SAMSUNG1/FCD3a_gra/SEEG';
Surgpath = '/media/mo/SAMSUNG1/FCD3a_gra/images';

sub = '18xiaozhiguo';side = 'lh';
s = SurfStatReadSurf({['/usr/local/freesurfer/subjects/fsaverage_sym/surf/',char(side),'.pial']});
Measures = {'.thickness_z_on_lh.sm.mgh';'.gradient_z_on_lh.sm.mgh';...
    '.area_z_on_lh.sm.mgh';'.complexity_z_on_lh.sm.mgh';'.gm_FLAIR_0.5_z_on_lh.sm.mgh';...
    '.gm_PET_0.5_z_on_lh.sm.mgh';'.w-g.pct_z_on_lh.sm.mgh';'.FLAIRgradient_z_on_lh.sm.mgh';'.PETgradient_z_on_lh.sm.mgh';};
 
Y_thick = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{1,1})]});
Y_gradient = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{2,1})]});
Y_area = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{3,1})]});
Y_comp = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{4,1})]});
Y_gmflair = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{5,1})]});
Y_gmpet = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{6,1})]});
Y_blur = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{7,1})]});
Y_flairgra = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{8,1})]});
Y_petgra = SurfStatReadData({[char(Neuropath),'/',char(sub),'/xhemi/surf/',char(side),char(Measures{9,1})]});
Y_em = SurfStatReadData({[char(SEEGpath),'/',char(sub),'/SEEG/SPM_EI_SZ1_80_250_5_0/mapping_EI_on_lh.mgh']});

% epi prediction
non_0 = find(Y_em~=0); em = Y_em(non_0); EM = em';
% EM_ori = em'; EM = normalize(EM_ori);
ima = [Y_thick',Y_gradient',Y_area',Y_comp',Y_gmflair',Y_gmpet',Y_flairgra',Y_petgra']; 
IMA = ima(non_0,:);
Mdl = fitrsvm(IMA,EM,'KernelFunction','gaussian');
% conv = Mdl.ConvergenceInfo.Converged;
% lStd = resubLoss(Mdl);
% iter = Mdl.NumIterations;
fit = predict(Mdl,ima);
surg = SurfStatReadData({[char(Surgpath),'/fcd3a/',char(sub),'/surg_mapping.mgh']});

figure;SurfStatView( fit, s, 'Predict');
figure;Y_em (Y_em < 0) = 0; SurfStatView(Y_em,s,'EI')
surg = -surg;surg = surg+1;
figure;SurfStatView(surg,s,'surg')

FIT = fit';
[R,P] = corrcoef (Y_em(non_0), FIT(non_0));
