function auto_reorient_T1(inputs)
nrun = size(inputs,2);
jobfile = fullfile(fileparts(mfilename('fullpath')),'auto_reorient_temp_job.m');
spm_pth = fileparts(which('spm'));
template = fullfile(spm_pth,'canonical','avg152T1.nii');
%% create job_file
fid = fopen(jobfile, 'w');
str.ref = 'matlabbatch{%i}.spm.tools.hmri.autoreor.reference = {''%s''};\n';
str.tpl = 'matlabbatch{%i}.spm.tools.hmri.autoreor.template = {''%s''};\n';
str.ot1 = 'matlabbatch{%i}.spm.tools.hmri.autoreor.other = {\n';
str.otf = '                                                ''%s''\n';
str.ot2 = '                                                };\n';
str.dir = 'matlabbatch{%i}.spm.tools.hmri.autoreor.output.indir = ''yes'';\n';
str.dep = 'matlabbatch{%i}.spm.tools.hmri.autoreor.dep = ''individual'';\n';
for crun = 1:nrun
    fprintf(fid,sprintf(str.ref,crun,char(inputs{8,crun}(1))));
    fprintf(fid,sprintf(str.tpl,crun,template));
    fprintf(fid,sprintf(str.ot1,crun));
    other = [];
    for row = 2:8
        if ~isempty(inputs{row,crun})
            if isempty(other)
                other = [inputs{row,crun}];
            else
                other = [other;inputs{row,crun}]; %#ok<AGROW>
            end
        end
    end
    for cf = 1:numel(other)
        fprintf(fid,sprintf(str.otf,char(other(cf))));
    end
    fprintf(fid,sprintf(str.ot2));
    fprintf(fid,sprintf(str.dir,crun));
    fprintf(fid,sprintf(str.dep,crun));
end
jobs = repmat({jobfile}, 1, 1);
spm('defaults', 'FMRI');
spm_jobman('run', jobs);
delete(jobfile)
end
