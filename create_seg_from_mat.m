function create_seg_from_mat(seg8_file)
res = load(seg8_file);
tc = false(numel(res.tpm),4);
tc(1:3) = true; % get GM, WM and CSF segmentations
bf = false(1,2);
df = [1 0];
% % path & file name modification, if files were copied to Supplementary path
% supplpath = spm_file(seg8_file,'path');
% respath = strrep(supplpath,'Supplementary','');
% [~,file] = fileparts(res.image(1).fname);
% file = strrep(file,'_outer_suppressed','');
% res.image(1).fname = fullfile(respath,[file '.nii']);
% % overwrite tpm filenames because of reading access restrictions...
% for n=1:6
%     res.tpm(n).fname = hmri_get_defaults('TPM');
% end
mrf = 1;
cleanup = 1;
% from spm_preproc_write8 itself (modified a bit):
bb = NaN(2,3);
vx = NaN;
[bb1,vx1] = spm_get_bbox(res.tpm(1), 'old');
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end
spm_preproc_write8(res,tc,bf,df,mrf,cleanup,bb,vx);
end