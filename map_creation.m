function map_creation(inputs,vendor,RFflag,B1flag)
% caller function for choosing the correct processing pipeline depending on
% available input and vendor specific defaults/B1 mapping
switch vendor
    case 'p'
        if contains(B1flag,'noB1')
            loc_conf = 'hmri_local_defaults_NISCI_THS_PH_UNICORT_all';% set UNICORT defaults
            job_fn = sprintf('u_%s_job.m',RFflag.name);
        else
            prep_type = get_metadata_val(inputs{5,1}{1,1},'PrepulseType');
            switch deblank(prep_type)
                case 'NO' % Yarnikh B1 mapping
                    loc_conf = 'hmri_local_defaults_NISCI_THS_PH_thrA10e8'; % set B1 mapping defaults
                    job_fn = sprintf('p_%s_job.m',RFflag.name);
                case 'SAT' % saturated double angle method -> likely no CSK for MT as well
                    loc_conf = 'hmri_local_defaults_NISCI_THS_PH_thrA10e8_MT30'; % set B1 mapping defaults
                    job_fn = sprintf('p_SDAM_%s_job.m',RFflag.name);
            end
        end
    case 's'
        if contains(B1flag,'noB1')
            loc_conf = 'hmri_local_defaults_NISCI_THS_nocleanup_seg_UNICORT_all';% set UNICORT defaults
            job_fn = sprintf('u_%s_job.m',RFflag.name);
        else
            loc_conf = 'hmri_local_defaults_NISCI_THS_nocleanup_seg'; % set B1 mapping defaults
            job_fn = sprintf('s_%s_job.m',RFflag.name);
        end
end
% collect inputs for SPM jobman
[cp,~,~] = fileparts(mfilename('fullpath'));
spm_inputs{1,1} = {fullfile(cp,'loc_conf',strcat(loc_conf,'.m'))};
spm_inputs{2,1} = inputs{1,1};
switch RFflag.name
    case 'RF_US'
        if contains(B1flag,'noB1')
            spm_inputs{3,1} = inputs{6,1};
            spm_inputs{4,1} = inputs{7,1};
            spm_inputs{5,1} = inputs{8,1};
        else
            spm_inputs{3,1} = inputs{5,1};
            spm_inputs{4,1} = inputs{6,1};
            spm_inputs{5,1} = inputs{7,1};
            spm_inputs{6,1} = inputs{8,1};
        end
    case 'RFonce'
        if contains(B1flag,'noB1')
            spm_inputs{3,1} = inputs{RFflag.RFcont,1};
            spm_inputs{4,1} = inputs{6,1};
            spm_inputs{5,1} = inputs{7,1};
            spm_inputs{6,1} = inputs{8,1};
        else
            spm_inputs{3,1} = inputs{RFflag.RFcont,1};
            spm_inputs{4,1} = inputs{5,1};
            spm_inputs{5,1} = inputs{6,1};
            spm_inputs{6,1} = inputs{7,1};
            spm_inputs{7,1} = inputs{8,1};
        end
    case 'RFsens'
        if contains(B1flag,'noB1')
            spm_inputs{3,1} = inputs{2,1};
            spm_inputs{4,1} = inputs{3,1};
            spm_inputs{5,1} = inputs{4,1};
            spm_inputs{6,1} = inputs{6,1};
            spm_inputs{7,1} = inputs{7,1};
            spm_inputs{8,1} = inputs{8,1};
        else
            spm_inputs{3,1} = inputs{2,1};
            spm_inputs{4,1} = inputs{3,1};
            spm_inputs{5,1} = inputs{4,1};
            spm_inputs{6,1} = inputs{5,1};
            spm_inputs{7,1} = inputs{6,1};
            spm_inputs{8,1} = inputs{7,1};
            spm_inputs{9,1} = inputs{8,1};
        end
end
jobfile = fullfile(fileparts(mfilename('fullpath')),'job_files',job_fn);
jobs = repmat(jobfile, 1, 1);
spm('defaults', 'FMRI');
spm_jobman('run', jobs, spm_inputs{:});
end
