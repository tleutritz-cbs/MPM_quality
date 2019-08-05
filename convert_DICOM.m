function convert_DICOM(dcm_path)
cd(dcm_path);
files = dir(strcat(strcat('**',filesep,'*')));
files = files(~[files.isdir]); % make sure to only have files in the list
files = fullfile({files.folder},{files.name});
files = reshape(files,[numel(files),1]);
spm_inputs{1,1} = files;
spm_inputs{2,1} = {fullfile(fileparts(dcm_path),'NIFTI')};
jobfile = fullfile(fileparts(mfilename('fullpath')),'job_files','dcm_conv_job.m');
jobs = repmat(jobfile, 1, 1);
spm('defaults', 'FMRI');
spm_jobman('run', jobs, spm_inputs{:});
end
