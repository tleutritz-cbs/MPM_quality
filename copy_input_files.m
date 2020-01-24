function inputs = copy_input_files(NIFTI_path,old_inputs,sub_folder)
inputs = old_inputs;
for inp_ind = 2:8
    for fil_ind = 1:numel(old_inputs{inp_ind,1})
        [cur_path,cur_file,ext] = fileparts(old_inputs{inp_ind,1}{fil_ind});
        fs_ind = strfind(cur_path,filesep); % indices of filesep to find folders
        ar = fullfile(cur_path(1:fs_ind(end)),sub_folder);
        if ~exist(ar,'dir')
            mkdir(ar);
        end
        
        % restructure inputs with new paths and copy files
        new_path = fullfile(ar,cur_path(fs_ind(end):end));
        new_file = fullfile(new_path,strcat(cur_file,ext));
        old_file = fullfile(cur_path,strcat(cur_file,ext));
        if contains(new_file,',1')
            new_file = strrep(new_file,',1','');
            old_file = strrep(old_file,',1','');
        end
        inputs{inp_ind,1}{fil_ind} = new_file;
        if ~exist(new_file,'file')
            if ~exist(new_path,'dir')
                mkdir(new_path);
            end
            copyfile(old_file,new_file); % explicitly specify destination folder,
            % otherwise just contents will end up there and not the folder structure itself
        end
        if ~exist(fullfile(new_path,strcat(cur_file,'.json')),'file')
            copyfile(strrep(old_file,'.nii','.json'),strrep(new_file,'.nii','.json'));
        end
    end
end
save(fullfile(NIFTI_path,strcat('inputs_',sub_folder,'.mat')),'inputs')
