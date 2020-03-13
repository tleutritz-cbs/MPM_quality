function QA_NIFTI(NIFTI_path)
[cp,~,~] = fileparts(mfilename('fullpath'));
atlas_ROI_dir = fullfile(cp,'ROIs');
%% protocol check
fprintf('NIFTI directory: %s\n',NIFTI_path)
rename_folders(NIFTI_path) % for a nicer order
fprintf('==== check protocol parameters ====\n')
[inputs,vendor] = check_protocol_and_create_inputs(NIFTI_path);
fprintf('==== report written to %s ====\n',fullfile(fileparts(NIFTI_path),'protocol_check.htm'))
fid = fopen(fullfile(NIFTI_path,'inputs.txt'),'w+');
for n=2:8, try fprintf(fid,'%s\n%s\n',inputs{n,1}{1},inputs{n,1}{2}); catch err, end, end
%% auto-reorient data for better segmentation
fprintf('==== reorient to MNI templates ====\n')
inputs = copy_input_files(NIFTI_path,inputs,'AR');
auto_reorient(inputs);
fid = fopen(fullfile(NIFTI_path,'inputs_AR.txt'),'w+');
for n=2:8, try fprintf(fid,'%s\n%s\n',inputs{n,1}{1},inputs{n,1}{2}); catch err, end, end
%% head-masking input data
fprintf('==== head mask data ====\n')
create_head_mask2(inputs{8,1}{1},0.4,10,false,true)
inputs = head_mask_input_data(NIFTI_path,inputs);
fid = fopen(fullfile(NIFTI_path,'inputs_HM.txt'),'w+');
for n=2:8, try fprintf(fid,'%s\n%s\n',inputs{n,1}{1},inputs{n,1}{2}); catch err, end, end
%% map creation depending on input
fprintf('==== create quantitative maps ====\n')
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
%% run mpm_quality script depending on B1 mapping method used
fprintf('==== running ROI based QA on maps ====\n')
if contains(B1flag,'noB1')
    mpm_quality_w_def_atlas_ROI_UNICORT_newTB_ext_v4b(fullfile(fileparts(NIFTI_path),'maps','Results'))
else
    mpm_quality_w_def_atlas_ROI_B1maps_newTB_ext_v4c(fullfile(fileparts(NIFTI_path),'maps','Results'))
end
% run multi-echo quality script
fprintf('==== running ROI based QA on multi-echo data ====\n')
mpm_quality_multiecho_data_atlas_ROI_v4(fullfile(fileparts(NIFTI_path),'maps','Results'))
report_mpm_qlt(fullfile(fileparts(NIFTI_path),'maps','Results'))
fprintf('==== report written to %s ====\n',fullfile(fileparts(NIFTI_path),'maps','Results','Supplementary','quality_report.htm'))
end
