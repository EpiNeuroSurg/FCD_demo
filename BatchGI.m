clear;clc;
Root='C:\Users\jiaji\Desktop\fcd3aSEEG';
%Ici on rentre les infos des patient
I=1;
Patient{I}.Name='forfig1';
Patient{I}.File{1}='SZ1';
Patient{I}.Baseline{1}='Baseline_SZ1';
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[80 250];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=5;
Patient{I}.sMRI=fullfile(Root,Patient{I}.Name,'imaging','wBNI_2017_Buchuntao.nii');
Patient{I}.Pre='';
Patient{I}.ElectrodeName='BNI_2017_Buchuntao_MNI_Name.txt';
Patient{I}.ElectrodePos='BNI_2017_Buchuntao_MNI_Pos.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ImaGIN
%En-dessous, on ne touche  plus
ThDelay=0.05;
for i0=1%2:length(Patient)
        %Convert TRC to SPM
        %[files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name,'Crises'),'.TRC');
    for i1=1:length(Patient{i0}.File)
        %Convert EDF to SPM
        clear S
        S.dataset=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.edf']);
        S.Atlas='Human';
        S.channel=[];
        S.coarse=1;
        S.FileOut=[Patient{i0}.File{i1}];
        S.SEEG='No';
        S.Bipolar='No';
%       S.NeventType=0;
        D = ImaGIN_spm_eeg_converteeg2mat(S);
        %Rename labels
        D=spm_eeg_load(fullfile(Root,Patient{i0}.Name,'SEEG',S.FileOut));
        Sensors=sensors(D,'EEG');
        tmp=Sensors.label;
        for i2=1:length(tmp)
            if length(tmp{i2})>4
                if strcmp(tmp{i2}(end-2:end),'ref')
                    tmp{i2}=tmp{i2}(1:end-3);
                end
            end
        end
        Sensors.label=tmp;
        D=sensors(D,'EEG',Sensors);
        D=chanlabels(D,1:length(Sensors.label),Sensors.label);
        save(D) 
        %specify channels type
        Cnames=chanlabels(D);
        Bad=[];
        EMG=[];
        ECG=[];
        Other=[];
        for i2=1:length(Cnames)
             a=lower(Cnames{i2});
             switch a
                 case {'emg','emg1','emg2','emg11','emg21'}
                     EMG=[EMG i2];
                  case {'dc10'}
                     Other=[Other i2];
                case {'ecg','ecg1','ecg2'}
                     ECG=[ECG i2];
                 otherwise
                     v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
                     if (isempty(v_num))
                         Bad=[Bad i2];
                     end
             end
             if strcmp(a(1),'$')
                     Other=[Other i2];
             end    
         end
         D=chantype(D,[Bad Other],'Other');
         D=chantype(D,[EMG],'EMG');
         D=chantype(D,[ECG],'ECG');
        save(D) 
        %Add electrode position
        clear S
        S.Fname=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']);
        S.filenamePos=fullfile(Root,Patient{i0}.Name,'implantation',Patient{i0}.ElectrodePos);
        S.filenameName=fullfile(Root,Patient{i0}.Name,'implantation',Patient{i0}.ElectrodeName);
        S.FileOut=S.Fname;
        D = ImaGIN_Electrode(S);
        %Longitudinal bipolar montage
        clear S
        S.Filename=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']);
        S.SaveFile=Patient{i0}.File{i1};
        S.FileOut=S.Filename;
        D = ImaGIN_BipolarMontage(S);
    end    
end

%ImaGIN add here marker of seizure onset, Baseline_
for i0=1  %2:length(Patient)
    %Set bad channels
    for i1=1:length(Patient{i0}.File)
        D=spm_eeg_load(fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']));
        D=badchannels(D,Patient{i0}.BadChannel,1);
        save(D);
    end
    for i1=1:length(Patient{i0}.Baseline)
        D=spm_eeg_load(fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.Baseline{i1} '.mat']));
        D=badchannels(D,Patient{i0}.BadChannel,1);
        save(D);
    end
end
%Compute wavelet
for i0=1:length(Patient)
    cd(fullfile(Root,Patient{i0}.Name,'SEEG'))
    TimeResolutionTF=0.05;
    TimeResolution=0.1;
    for i1=1:length(Patient{i0}.File)
        clear SS
        SS.D=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']);
        D=spm_eeg_load(SS.D);
        SS.Synchro='No';
        SS.Pre='';
        SS.Method='Morlet wavelet';
        SS.frequencies=10:3:230;
        SS.FactMod=0;
        SS.Mfactor=20;
        SS.Width=0;
        SS.TimeWindow=-10:.2:10;
        SS.TimeWindowWidth=SS.TimeWindow(2)-SS.TimeWindow(1);
        SS.Coarse=0;
        SS.channels=1:D.nchannels;
        SS.TimeResolution=TimeResolutionTF;
        ImaGIN_spm_eeg_tf(SS);
        clear SS2
        SS2.D=fullfile(D.path,['w1_' SS.Pre '_' D.fname]);
        SS2.B=[-10 -1];
        ImaGIN_NormaliseTF(SS2);
    end
        
    [files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name),'^nw1.*\.mat$');
    if size(files,1)>1
            cd(fullfile(Root,Patient{i0}.Name,'Crises'))
        clear S
    %     for i1=1:size(files,1)
    %         S.D{i1}=deblank(files(i1,:));
    %     end
        S.D=files;
        S.Method='Mean';
        S.NewName='Mean';
        D = ImaGIN_AverageTF(S);
    end
end
%Epileptogenicity
for i0=1:length(Patient)
    clear S
    for i1=1:length(Patient{i0}.File)
        if i1==1
            S.D=fullfile(Root,[Patient{i0}.Name],'SEEG',[Patient{i0}.File{i1} '.mat']);
            S.B=fullfile(Root,[Patient{i0}.Name],'SEEG',[Patient{i0}.Baseline{i1} '.mat']);
        else
            S.D=char(S.D,fullfile(Root,[Patient{i0}.Name],'SEEG',[Patient{i0}.File{i1} '.mat']));
            S.B=char(S.B,fullfile(Root,[Patient{i0}.Name],'SEEG',[Patient{i0}.Baseline{i1} '.mat']));
        end
    end
    S.TimeWindow=[0:0.01:Patient{i0}.TimeConstant+1+max(Patient{i0}.Latency)];
    S.FreqBand=Patient{i0}.FreqBand;
    S.HorizonT=Patient{i0}.TimeConstant;
    S.BadChannel=Patient{i0}.BadChannel;
    try
        S.FileName=Patient{i0}.Pre;
    catch
        S.FileName='';
    end
    S.Latency=Patient{i0}.Latency;
    S.TimeResolution=0.1;
    S.ThDelay=ThDelay;
    S.Atlas='Human';
    S.AR=0;
    S.Latency=0;
    S.sMRI=Patient{i0}.sMRI;
    S.CorticalMesh=1;
    ImaGIN_Epileptogenicity(S);
end

% % Compute multi-taper
% for i0=3%5:length(Patient)
%     cd(fullfile(Root,Patient{i0}.Name,'SEEG'))
%     maxTime=10;
% %     maxTime=2;
%     minTime=-10;
%     for i1=1:length(Patient{i0}.File)
%         clear SS
%         SS.D=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']);
%         D=spm_eeg_load(SS.D);
%         minTime=max([minTime min(time(D))+0.5]);
%         maxTime=min([maxTime max(time(D))-0.5]);
%     end
%     for i1=1:length(Patient{i0}.File)
%         clear SS
%         SS.D=fullfile(Root,Patient{i0}.Name,'SEEG',[Patient{i0}.File{i1} '.mat']);
%         D=spm_eeg_load(SS.D);
%         SS.Pre='';
%         SS.Method='Multitaper';
%         SS.Taper='hanning';
%         SS.TimeResolution=0.1;
%         SS.frequencies=10:3:230;
%         SS.FactMod=10;
%         SS.TimeWindow=[minTime maxTime];
%         SS.TimeWindowWidth=1;
%         SS.channels=1:D.nchannels;
%         SS.NSegments=1;
%         ImaGIN_spm_eeg_tf(SS);
%         clear SS2
%         SS2.D=fullfile(D.path,['m1_' SS.Pre '_' D.fname]);
%         SS2.B=[-6-1];
%         ImaGIN_NormaliseTF(SS2);
%     end
%     
%     
%     [files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name),'^nm1.*\.mat$');
%     if size(files,1)>1
%             cd(fullfile(Root,Patient{i0}.Name))
%         clear S
%     %     for i1=1:size(files,1)
%     %         S.D{i1}=deblank(files(i1,:));
%     %     end
%         S.D=files;
%         S.Method='Mean';
%         S.NewName='Mean';
%         D = ImaGIN_AverageTF(S);
%     end
% end
