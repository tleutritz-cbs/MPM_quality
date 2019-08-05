function create_head_mask2(file,clp,ndier,skip_dilation_erosion,skip_inter)
% according to qap pipeline (see inline hints)
% adapted to Matlab/SPM abilities
% input: file - filename of NIFTI-file to be masked (mandatory)
%        optional: clp - parameter for clip level, default = 0.4
%                        lower values give bigger outlines
%                        bigger values are better for noisy data
%                  ndier - number of dilations/erosions, default = 10
%                  skip_dilation_erosion - default = false
%                        set to true, if clipped outline is already good
%                        otherwise noise can be reduced a lot (espescially
%                        Philips data requires that!)
%                  skip_inter - default = true
%                        set to false, if you want the clipped mask (first 
%                        step of the processing) and the to be written as output
% output: file_clipped.nii (clipped file -> easy mask -> check first)
%         file_headmask.nii - whole head mask
%         file_ghostmask.nii - mask with 4 values: 
%              ghost = 0, background = 1, 2 = ghost + head, 3 = head only
% by Tobias Leutritz, 2018-05-30
% 2018-09-26: refined erosion/dilation steps and output structure
%             raised ndier and modified padding
% 2018-09-27: add option to skip dilation/erosion (default = false); 
%             input arguments; meta-data
% 2018-11-01: add option to skip intermediate results (default = true)
% 2018-11-09: changed defaults to be more robust: clp=0.4; ndier = 10

%% settings
if nargin < 5
    skip_inter = true;
    if nargin < 4
        skip_dilation_erosion = false;
        if nargin < 3
            ndier = 10;
            if nargin < 2
                clp = 0.4;
            end
        end
    end
end

% define kernel for dilation/erosion:
afni_nh(:,:,1) = [0 1 0; 1 1 1; 0 1 0];
afni_nh(:,:,2) = [1 1 1; 1 1 1; 1 1 1];
afni_nh(:,:,3) = [0 1 0; 1 1 1; 0 1 0];

%fill_kernel = afni_nh;

pd = ndier + 2; % offset for padding

% set variables for metadata:
proc_par.clp = clp; proc_par.num_of_dilation_erosion = ndier;
proc_par.skip_dilation_erosion = skip_dilation_erosion;
json = struct('extended',false,'separate',true,'anonym','none',...
    'overwrite',true, 'indent','\t');

%% load data
hdr = spm_vol(file);
tmp_dim = hdr.dim; % save dimensions for reset after padding
input_file = hdr.fname;
dta = spm_read_vols(hdr);
hdr.pinfo(1) = 1; % reset scaling
xdim = size(dta,1);
ydim = size(dta,2);
zdim = size(dta,3);

%% mask data
cll = find_clip_level(dta,clp);
dta(dta >= cll) = 1;
dta(dta ~= 1) = 0;
if ~skip_inter
    hdr.fname = strrep(hdr.fname,'.nii','_clipped.nii');
    hdr.dt(1) = 2; % set data type
    spm_write_vol(hdr,dta);
    
    % Set and write metadata
    Output_hdr = init_mpm_output_metadata(input_file, proc_par);
    Output_hdr.history.output.imtype = 'Mask (clipped)';
    set_metadata(hdr.fname,Output_hdr,json);
end

% h = msgbox(sprintf('Have a look at the clipped file and continue, if noise should be removed:','hMRI toolbox','modal');
% uiwait(h);

if ~skip_dilation_erosion
    %% find contour and get rid of noise voxels
    tmp = imgradient3(dta);
    if ~skip_inter
        hdr.fname = strrep(hdr.fname,'.nii','_imgrad.nii');
        spm_write_vol(hdr,tmp);
    end
    
    cll = ceil(find_clip_level(tmp,0.83)); % this value needs to be fixed (tested at Philips data)
    tmx = max(tmp(:));
    tmp(tmp == 0) = tmx+1;
    tmp(tmp <= cll) = 1;
    dta(tmp == 1) = 0;
    % tmp(tmp < cll) = 0;
    % tmp(tmp >= cll) = 1;
    % dta = tmp;
    
    if ~skip_inter
        hdr.fname = strrep(hdr.fname,'.nii','_clipped.nii');
        spm_write_vol(hdr,dta);
    end
    
    tmp(tmp > tmx) = 0;
    tmp(tmp == 1) = 0;
    dta = dta+tmp; % add up contour
    dta = imbinarize(dta);
    % % dta = tmp;
    % hdr.fname = strrep(hdr.fname,'.nii','_addup.nii');
    % spm_write_vol(hdr,dta);
    
    %% dilate and erode ndier-times in order to get rid of unwanted parts, partially fill
    % makes use of hints from:
    % https://blogs.mathworks.com/steve/2007/03/23/pad-values-in-dilation-and-erosion/
    % https://de.mathworks.com/help/images/operations-that-combine-dilation-and-erosion.html
    %% preparation
    dta_padded = padarray(dta,[pd pd pd], 0);
    dta_padded(:,1:2,:) = 1; % to avoid cutted neck-slice
    hdr.dim = tmp_dim + 2*pd;
    
    %% initial dilation to strengthen contour line 
    for n = 1:1
        dta_padded = imdilate(dta_padded,afni_nh);
        %     hdr.fname = strrep(hdr.fname,'.nii',sprintf('_dilated_%i.nii',n));
        %     spm_write_vol(hdr,dta_padded);
    end

    %% and afterwards get rid of remaining noise
    for n = 1:2
        dta_padded = imerode(dta_padded,afni_nh);
        %     hdr.fname = strrep(hdr.fname,'.nii',sprintf('_eroded_%i.nii',n));
        %     spm_write_vol(hdr,dta_padded);
    end
    
    %% now do the repetaed dilation/erosion to fill up
    for n = 1:ndier
        dta_padded = imdilate(dta_padded,afni_nh);
        %     hdr.fname = strrep(hdr.fname,'.nii',sprintf('_dilated_%i.nii',n));
        %     spm_write_vol(hdr,dta_padded);
    end
    for n = 1:ndier
        dta_padded = imerode(dta_padded,afni_nh);
%             hdr.fname = strrep(hdr.fname,'.nii',sprintf('_eroded_%i.nii',n));
%             spm_write_vol(hdr,dta_padded);
    end
    
%     %% fill holes
%     dta_padded = imclose(dta_padded,afni_nh);
%     if ~skip_inter
%         hdr.fname = strrep(hdr.fname,'.nii','_closed.nii');
%         spm_write_vol(hdr,dta_padded);
%     end
    
    %% undo padding
    hdr.dim = tmp_dim;
    dta = dta_padded((pd+1):(size(dta_padded,1)-pd), ...
        (pd+1):(size(dta_padded,2)-pd), ...
        (pd+1):(size(dta_padded,3)-pd));
    
%     %% additional filling
%     %for slice = 1:xdim, dta(slice,:,:) = imfill(dta(slice,:,:)); end
%     for slice = 1:xdim, dta(slice,:,:) = imclose(dta(slice,:,:),fill_kernel); end
%     hdr.fname = strrep(hdr.fname,'.nii','_x_filled.nii');
%     spm_write_vol(hdr,dta);
% %     for slice = 1:ydim, dta(:,slice,:) = imfill(squeeze(dta(:,slice,:))); end
%     for slice = 1:ydim, dta(:,slice,:) = imclose(dta(:,slice,:),fill_kernel); end
%     hdr.fname = strrep(hdr.fname,'.nii','_y_filled.nii');
%     spm_write_vol(hdr,dta);
% %     for slice = 1:zdim, dta(:,:,slice) = imfill(dta(:,:,slice)); end
%     for slice = 1:zdim, dta(slice,:,:) = imclose(dta(slice,:,:),fill_kernel); end
%     hdr.fname = strrep(hdr.fname,'.nii','_z_filled.nii');
%     spm_write_vol(hdr,dta);
    
%     hdr.fname = strrep(hdr.fname,'_clipped','');
    
end % of dilation erosion process

%% write headmask
hdr.fname = strrep(hdr.fname,'.nii','_headmask.nii');
hdr.dt(1) = 2; % set data type
spm_write_vol(hdr,dta);
% Set and write metadata
Output_hdr = init_mpm_output_metadata(input_file, proc_par);
Output_hdr.history.output.imtype = 'Whole Head Mask (clipped)';
if ~skip_dilation_erosion
    Output_hdr.history.output.imtype = 'Whole Head Mask';
end
set_metadata(hdr.fname,Output_hdr,json);

%% create additionally a ghost mask, asssuming  phase encoding in y direction
% again: adapted from qap pipeline, namely spatial_qc.py
n2_mask_data = zeros(size(dta));
n2 = floor(size(dta,1)/2);
for xcnt = 1:xdim
    for ycnt = 1:ydim
        for zcnt = 1:zdim
            if xcnt <= n2
                n2_mask_data(xcnt,ycnt,zcnt) = 1-dta(xcnt+n2-1,ycnt,zcnt);
%                 if xcnt == 1
%                     n2_mask_data(xcnt,ycnt,zcnt) = 1;
%                 end
            else
                n2_mask_data(xcnt,ycnt,zcnt) = 1-dta(xcnt-n2,ycnt,zcnt);
%                 if xcnt == xdim
%                     n2_mask_data(xcnt,ycnt,zcnt) = 1;
%                 end
            end
        end
    end
end

% create non-ghost background region
n2_mask_data = n2_mask_data + dta;

% add up original mask 
n2_mask_data = n2_mask_data + dta;

%% write ghostmask
hdr.fname = strrep(hdr.fname,'headmask','ghostmask');
spm_write_vol(hdr,n2_mask_data);
% Set and write metadata
Output_hdr = init_mpm_output_metadata(input_file, proc_par);
Output_hdr.history.output.imtype = 'Ghost Mask (clipped)';
if ~skip_dilation_erosion
    Output_hdr.history.output.imtype = 'Ghost Mask';
end
Output_hdr.history.output.units = 'INT 0 = ghost, 1= background, 2 = ghost + head, 3 = head only';
set_metadata(hdr.fname,Output_hdr,json);

end


function metastruc = init_mpm_output_metadata(input_file, proc_params)

proc.descrip = [mfilename '.m - head mask creation'];
proc.params = proc_params;
proc.version = 'v2';
output.imtype = 'mask';
output.units = 'INT';

metastruc = init_output_metadata_structure(input_file, proc, output);

end
