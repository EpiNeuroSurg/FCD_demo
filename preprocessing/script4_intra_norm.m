clear;clc;
SUBJECTS_DIR = '/media/mo/SAMSUNG1/2020/temporal_cases/recon/control';
cd(SUBJECTS_DIR)
setenv SUBJECTS_DIR .
addpath /usr/local/freesurfer/matlab
Subs=dir;
subs=cell(length(Subs),1);
for s = 1:length(Subs)
    subs{s}=Subs(s).name;
end
% Measures={'.gradient.sm5.mgh';'.gm_FLAIR_0.5.sm5.mgh';'.myeline_2.sm5.mgh';...
%         '.thickness.sm5.mgh';'.w-g.pct.sm5.mgh';'.complexity.sm5.mgh';};

Measures={'.myeline_1.sm5.mgh'};
NumberOfMeasures=length(Measures);
for s=3:length(subs)
    sub=subs(s);
    sub=cell2mat(sub);   
    Cortexlh=read_label(['',sub,''],['lh.cortex']);
    Cortexlh=Cortexlh(:,1)+1;
    Cortexrh=read_label(['',sub,''],['rh.cortex']);
    Cortexrh=Cortexrh(:,1)+1;
    for L=1:NumberOfMeasures
        M_lh=MRIread(['',sub,'/surf/lh',Measures{L},'']);
        M_rh=MRIread(['',sub,'/surf/rh',Measures{L},'']); 
        % Normalise by mean and std
        M=[M_lh.vol(Cortexlh),M_rh.vol(Cortexrh)];
        meanM=mean(M);
        stdM=std(M);
        M_lh.vol(1,:)=(M_lh.vol(1,:)-meanM)./stdM;
        M_rh.vol(1,:)=(M_rh.vol(1,:)-meanM)./stdM;
        %create string for saving output after removing last 9 characters.
        Meas_output = Measures{L}(1:end-8);
        MRIwrite(M_lh,['',sub,'/surf/lh',Meas_output,'_z.sm.mgh']);
        MRIwrite(M_rh,['',sub,'/surf/rh',Meas_output,'_z.sm.mgh']);
    end
end   