function inputs = head_mask_input_data(NIFTI_path,inputs)
mask = spm_vol(strrep(inputs{8,1}{1},'.nii','_headmask.nii'));
inputs = copy_input_files(NIFTI_path,inputs,'HM');
for inp_ind = 2:8
    for fil_ind = 1:numel(inputs{inp_ind,1})
        [cur_path,cur_file,ext] = fileparts(inputs{inp_ind,1}{fil_ind});
        spm_imcalc([spm_vol(fullfile(cur_path,strcat(cur_file,ext))) mask],fullfile(cur_path,strcat(cur_file,ext)),'i1.*i2');
    end
end

end
