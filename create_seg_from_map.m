function create_seg_from_map(file)
Vsave = spm_vol(file); 
% code from hmri_create_MTPRot:
MTtemp = spm_read_vols(Vsave);
% The 5 outer voxels in all directions are nulled in order to remove
% artefactual effects from the MT map on segmentation:
MTtemp(1:5,:,:)=0; MTtemp(end-5:end,:,:)=0;
MTtemp(:,1:5,:)=0; MTtemp(:,end-5:end,:)=0;
MTtemp(:,:,1:5)=0; MTtemp(:,:,end-5:end)=0;
Vsave.fname = spm_file(Vsave.fname,'suffix','_outer_suppressed');
spm_write_vol(Vsave,MTtemp);

% use unified segmentation with uniform defaults across the toobox:
job_brainmask = hmri_get_defaults('segment');
job_brainmask.channel.vols = {Vsave.fname};
job_brainmask.channel.write = [0 0]; % no need to write BiasField nor BiasCorrected image
job_brainmask.tissue(4).native = [0 0]; % don't need c4/c5 segments
job_brainmask.tissue(5).native = [0 0];
job_brainmask.warp.write = [1 0]; % add inverse deformation field
spm_preproc_run(job_brainmask);
end