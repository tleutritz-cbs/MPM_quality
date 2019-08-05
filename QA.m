function QA(DICOM_path)
[cp,~,~] = fileparts(mfilename('fullpath'));
atlas_ROI_dir = fullfile(cp,'ROIs');
%% DICOM conversion and protocol check
%convert_DICOM(DICOM_path);
NIFTI_path = DICOM_path;%fullfile(fileparts(DICOM_path),'NIFTI');
rename_folders(NIFTI_path) % for a nicer order
[inputs,vendor] = check_protocol_and_create_inputs(NIFTI_path);
%% auto-reorient data for better segmentation
auto_reorient(inputs);
%% map creation depending on input
empty = zeros(1:5); B1flag = 'B1';
for n = 2:5 % check inputs for B1 and RFsens
   empty(n) = isempty(inputs{n,1});
end
if sum(empty(2:4)) == 0
    RFflag.name = 'RFsens';
elseif sum(empty(2:4)) == 3
    RFflag.name = 'RF_US';
elseif sum(empty(2:4)) < 3
    RFflag.name = 'RFonce';
    RFcont = find(empty(2:4));
    RFflag.RFcont = RFcont(1); % limit to first appearing contrast
end
if empty(5)
    B1flag = 'noB1';
end
map_creation(inputs,vendor,RFflag,B1flag)
%% map and multi-echo data quality
% unzip ROIs
gunzip_or_gzip(atlas_ROI_dir,'nii')
% run mpm_quality script depending on B1 mapping method used
if contains(B1flag,'noB1')
    mpm_quality_w_def_atlas_ROI_UNICORT_newTB_ext_v4b(fullfile(fileparts(NIFTI_path),'maps','Results'))
else
    mpm_quality_w_def_atlas_ROI_B1maps_newTB_ext_v4c(fullfile(fileparts(NIFTI_path),'maps','Results'))
end
% run multi-echo quality script
mpm_quality_multiecho_data_atlas_ROI_v4(fullfile(fileparts(NIFTI_path),'maps','Results'))
end
