clear;clc;
addpath /home/mo/Desktop/surfstat
addpath /media/mo/SAMSUNG1/2020/temporal_cases/FCD3a_codes
age=[19;24;22;25;34;54;12;34];Age=term(age);
duration=[3;5;6;1;3;6;5;6];Duration=term(duration);
group=[{'female'},{'female'},{'female'},{'female'},{'male'},{'male'},{'male'},{'male'}];Group=term(group);
M=1+Group+Age+Duration;

feature='thickness';
fcd3a=({...
    '/home/mo/Desktop/recon/FCD3a/01wuxuehui/xhemi/surf/rh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/02zhouhezhen/xhemi/surf/lh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/03zhanghanyu/xhemi/surf/lh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/04caoyan/xhemi/surf/lh.thickness_z_on_lh.sm.mgh';...
    });
control=({...
    '/home/mo/Desktop/recon/FCD3a/01wuxuehui/xhemi/surf/lh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/02zhouhezhen/xhemi/surf/rh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/03zhanghanyu/xhemi/surf/rh.thickness_z_on_lh.sm.mgh';...
    '/home/mo/Desktop/recon/FCD3a/04caoyan/xhemi/surf/rh.thickness_z_on_lh.sm.mgh';...
    });
contrast = Group.female - Group.male;
Y = SurfStatReadData([fcd3a;control]);
filesleft  =  SurfStatListDir('/home/mo/Desktop/recon/FCD3a/fsaverage_sym/surf/lh.pial' );
avsurf = SurfStatAvSurf( filesleft );
slm = SurfStatLinMod( Y, M, avsurf );
slm = SurfStatT( slm, contrast );
 
[a,l_label,b] = read_annotation('/home/mo/Desktop/recon/FCD3a/fsaverage_sym/label/lh.aparc.annot');
mask(1,:) = 0; 
a1 = find(l_label == 2647065); a2 = find(l_label == 660700); a3 = find(l_label == 9231540); a4 = find(l_label == 7874740); 
a5 = find(l_label == 3302560); a6 = find(l_label == 3988500); a7 = find(l_label == 14474380); a8 = find(l_label == 11146310); 
a9 = find(l_label == 13145750); A12 = union(a1,a2);A123 = union(A12,a3);A1234 = union(A123,a4);
A12345 = union(A1234,a5);A123456 = union(A12345,a6);A1234567 = union(A123456,a7);
A12345678 = union(A1234567,a8);A123456789 = union(A12345678,a9);
mask(1,A123456789) = 1;
col=find(mask==1);
%% permutation
[epval,t_orig,crit_t,est_alpha,seed_state]=mult_comp_perm_t2(Y(1:4,col),Y(5:8,col),1000);
ppval=zeros(1,length(mask));
ppval(1,col)=epval;
SurfStatView(ppval,avsurf, 'P value');
% SurfStatColLim([0 1])
%% cohen's d
cohen=zeros(1,length(mask));num=0;
for i=1:length(col)
    num=col(1,i);
    cohen(1,num)=computeCohen_d(Y(1:4,num),Y(5:8,num));
end
% cohen(1,col) = computeCohen_d(Y(1:4,:),Y(5:8,:)); 
SurfStatView(cohen,avsurf,'Cohens d');
%% reproducibility
rpval=zeros(1000,163842);prob=rpval(1,:);
for i = 1:1000
    permutation=randperm(8);
    [h,p]=ttest2(Y(permutation(1:4),col),Y(permutation(5:8),col));
    rpval(i,col)=h;
end
for j = 1:163842
    temp=rpval(:,j);
    num=find(temp>0 & temp<0.05);
    len=length(num);
    prob(1,j)=len./1000;
end
SurfStatView(prob,avsurf,'reproducibility');
%% partial correlation
pcrho=zeros(1,163842);pcval=pcrho;
for i = 1:163842
    [rho,pval] = partialcorr(Y(:,i),seve,age);
    pcrho(1,i)=rho;
    pcval(1,i)=pval;
end
SurfStatView(pcrho,avsurf,'Rho');
SurfStatView(pcval,avsurf,'p');
