function mpm_quality_w_def_atlas_ROI_B1maps_newTB_ext_v4c(pnames)

% This function calculates some important measures and compares them to 
% values from the NISCI THS study with the same protocol at 3 clinical sites,
% measured at 5 healthy subjects including test and retest measurements.
% Tobias Leutritz, 12.01.2018
% 2018-08-23: use of segmentations from MPMCalc, flexibilize vox. size
% 2018-08-24: fix labelling warning (Ignoring extra legend entries.)
%             + add whole brain histograms
% 2018-08-27: restructure struct to have values with ROI name
% 2018-11-07: fix errounous struct
% 2018-11-15: simplify segmentation input
% 2018-12-21: MT -> MTsat due to NI paper revision
% 2019-05-20: add MTsat segmentation
% 2019-08-07: change ref. values to THS results

% get files to be analyzed (segmented data from subjects, each in a folder)
if nargin < 1
    pnames = cfg_getfile(Inf,'dir','Select folders with segmented MPM images to analyse');
end

% get number of selected path 
[npth,~] = size(pnames);

% Define maps to be read as well as GM, WM ,and CSF segmentations (i.e. c1*,
% c2*, and c3*] and the {inverse} deformation field [{i}y_...] to unwarp the
% ROI masks from normalised space back to the maps space in order to avoid 
% manual ROI definition for each subject.
% *seg8.mat file can be sufficient to reproduce missing files.
map.name_suffix = {'R1' 'PD' 'R2s_OLS' 'MTsat' 'MTsat' 'MTsat' 'MTsat' 'MTsat' 'T1w_OLSfit_TEzero'};
% #                  1    2    3      4       5    6    7    8    9    10 
map.name_prefix = {'' '' '' '' 'c1' 'c2' 'iy_' 'c3' ''};
% #                 1  2  3  4  5  6    7    8     9   10
MPM.name = {'R1' 'PD' 'R2s_OLS' 'MT' 'T1w'};
% #           1   2     4        5    6
% numbering of MPMs according to map.name_suffix
for nm = 1:numel(MPM.name)
    mpname = MPM.name{nm};
    pos.(mpname) = find(contains(map.name_suffix,mpname).*cellfun('isempty',map.name_prefix));
    if numel(pos.(mpname)) > 1
        pos.(mpname) = pos.(mpname)(1); % fetch the doubling of R2s...
    end
end
% numbering of segmentation files according to map.name_prefix
pos.c1 = find(contains(map.name_prefix,'c1'));
pos.c2 = find(contains(map.name_prefix,'c2'));
pos.c3 = find(contains(map.name_prefix,'c3'));
pos.idef = find(contains(map.name_prefix,'iy_'));
map.subfolder = {'' '' '' '' '' '' '' '' 'Supplementary'};
% #              1  2   3              4  5  6  7  8  9  10
mc = numel(map.name_suffix);
[cp,~,~] = fileparts(mfilename('fullpath'));
atlas_ROI_dir = fullfile(cp,'ROIs');

% define ROIs and MPM names for evaluation
ROI.name = {'GM' 'WM' 'GM_CN'  'WM_CC' ...
            'WM_CST' 'GM_S1' 'GM_M1' ...
            'GM_cerebellum' 'WM_cerebellum' 'GM_M1S1' ...
            'GM_Hippocampi' 'wholebrain'};
pos.wb = numel(ROI.name);

% reference values ordered as GM, CN, WM, CC as in Weiskopf et al. 2013
mean_lit.R1 = [0.62388, 0.64858, 0.98738, 1.09377];
mean_lit.R2 = [16.10719, 19.19238, 22.13997, 25.4136];
mean_lit.MT = [2.05309, 2.07694, 4.00709, 4.52533];
mean_lit.T1 = [301.80378, 332.90888, 384.23419, 417.14780];
mean_lit.PD = [80.46574, 88.21052, 68.76486, 68.79135];
cov_lit.R1 = [0.13273, 0.06633, 0.10520, 0.05873];
cov_lit.R2 = [0.38411, 0.17581, 0.18344, 0.17953];
cov_lit.MT = [0.16252, 0.09130, 0.10888, 0.09366];
cov_lit.T1 = [0.10100, 0.04789, 0.07117, 0.04159];
cov_lit.PD = [0.13670, 0.03881, 0.12127, 0.04635];

% initialise cell variables
vol_info = cell(mc,1);
%vol_data = cell(mc,1);
mnames = cell(mc,1);

%% loop over given path
for pn = 1:npth

    % change working directory to current subject folder
    cd(char(pnames(pn,:)));
    if contains(cd,'UNICORT')
        mpm_quality_w_def_atlas_ROI_UNICORT_newTB_ext_v4b(cd);
        continue
    end
    curdir = cd;
    outdir = fullfile(curdir,'MPM_Quality4');
    if ~exist(outdir,'dir')
        mkdir(outdir);
    elseif exist(fullfile(outdir,'measures_ext_v4.mat'),'file') == 2
        continue % skip further calculation, since measurment was completed
    end

    % extract basename of files
    bn_key = '_PD.nii';
    bnf = dir(strcat('*',bn_key));
    bn = strrep(bnf.name,bn_key,'');

    % loop over maps to be read in
    for mn = 1:mc
        if ~isempty(map.name_prefix{mn})
            if ~isempty(dir('../MPMCalc'))
                % create file names of segmentation files
                mnames{mn} = fullfile(strrep(curdir,'Results','MPMCalc'), ...
                    strcat(map.name_prefix{mn},bn,'_MTsat_outer_suppressed.nii'));
                idefcheck; % call function to check/create inverse deformation field
            elseif isempty(dir('c1*')) % create segmentations when missing
                create_seg_from_map(fullfile(curdir,strcat(bn,'_MTsat.nii')));
            elseif ~isempty(dir('c1*')) 
                % basename of files
                bnf = dir('c1*');
                bnc = strrep(bnf.name,'c1','');
                mnames{mn} = fullfile(curdir,strcat(map.name_prefix{mn},bnc)); 
                idefcheck;
            end
        else
            % create file names of maps
            mnames{mn} = fullfile(curdir,map.subfolder{mn}, ...
                strcat(map.name_prefix{mn},bn,'_',map.name_suffix{mn},'.nii'));
        end
        % create SPM handles
        vol_info{mn,1} = spm_vol(mnames{mn});
    end
    
    cd(outdir);
    
    % catch incomplete run
    if ~exist('temp.nii','file')
        
        %% binarise segmented WM and GM to get masks
        ROI_info.GM = vol_info{pos.c1,1};
        ROI_info.GM.fname = fullfile(outdir,'GM.nii');
        ROI_info.GM.descrip = 'GM mask';
        ROI_info.GM = spm_imcalc(vol_info{pos.c1,1},ROI_info.GM,'i1>0.99');
        
        ROI_info.WM = vol_info{pos.c2,1};
        ROI_info.WM.fname = fullfile(outdir,'WM.nii');
        ROI_info.WM.descrip = 'WM mask';
        ROI_info.WM = spm_imcalc(vol_info{pos.c2,1},ROI_info.WM,'i1>0.99');

        ROI_info.CSF = vol_info{pos.c3,1};
        ROI_info.CSF.fname = fullfile(outdir,'CSF.nii');
        ROI_info.CSF.descrip = 'CSF mask';
        ROI_info.CSF = spm_imcalc(vol_info{pos.c3,1},ROI_info.CSF,'i1>0.99');

        ROI_info.wholebrain = vol_info{pos.c3,1};
        ROI_info.wholebrain.fname = fullfile(outdir,'wholebrain.nii');
        ROI_info.wholebrain.descrip = 'whole brain mask';
        ROI_info.wholebrain = spm_imcalc([ROI_info.GM ROI_info.WM ROI_info.CSF], ...
            ROI_info.wholebrain,'(i1>.9)+(i2>.9)+(i3>.9)');
        
%         %% create inverse deformation field
%         % prepare options for spm_deformations
%         clear job;
%         
%         savedef.savedef.savedir.savepwd = 1;
%         savedef.savedef.ofname = strcat('i',bn,'_MT.nii');
%         job.out = {savedef};
%         
%         def.def = {vol_info{pos.idef,1}.fname};
%         inv.inv.comp = {def};
%         inv.inv.space = {vol_info{5,1}.fname};
%         job.comp = {inv};
%         
%         % run spm_deformations and set reference to output
%         spm_deformations(job);
%         vol_info{pos.idef,1} = spm_vol(fullfile(outdir,strcat('y_i',bn,'_MT.nii')));
        
        % loop over ROIs, skip GM/WM/WB
        for nr = 3:numel(ROI.name)-1
            % get ROI name
            rname = ROI.name{nr};
            ROI_info.(rname).fname = fullfile(outdir,strcat(rname,'.nii'));
            if contains(rname,'WM')
                copyfile(fullfile(atlas_ROI_dir,strcat(strrep(rname,'WM_',''),'.nii')),ROI_info.(rname).fname);
            elseif contains(rname,'GM')
                copyfile(fullfile(atlas_ROI_dir,strcat(strrep(rname,'GM_',''),'.nii')),ROI_info.(rname).fname);
            else
                copyfile(fullfile(atlas_ROI_dir,strcat(rname,'.nii')),ROI_info.(rname).fname);
            end
            
            %% rewarp masks to bring them back into the space of the MPMs
            % prepare options for spm_run_norm
            clear job;
            P = spm_imatrix(vol_info{1,1}.mat);
            job.woptions.vox = abs(P(7:9)); % voxel resolution
            job.woptions.bb = spm_get_bbox(char(vol_info{1,1}.fname), 'fv');
            job.woptions.interp = 4;
            job.woptions.prefix = 'inv';
            job.subj.resample = {ROI_info.(rname).fname};
            job.subj.def = {vol_info{pos.idef,1}.fname};
            
            % run rewarping
            spm_run_norm(job);
            
            %% reorientate rewarped masks to fit MPMs
            % rewrite file info after rewarping
            ROI_info.(rname) = spm_vol(fullfile(outdir,strcat('inv',rname,'.nii')));
            ROI_info.(rname).descrip = char(strcat(rname,' mask'));
            % prepare and run reslice
            flags.which  = 1;
            flags.mean = false;
            spm_reslice(char(vol_info{1,1}.fname,ROI_info.(rname).fname), flags);
            % rewrite file info after resclicing
            ROI_info.(rname) = spm_vol(fullfile(outdir,strcat('rinv',rname,'.nii')));
            
            %% rescale with low threshold to get rid of noise
            % and mask with GM/WM to avoid too big ROIs
            temp_info = ROI_info.(rname);
            %temp_info.fname = strcat(rname,'.nii');
            if strfind(rname,'WM') > 0
                ROI_info.(rname) = spm_imcalc([ROI_info.WM temp_info],...
                    ROI_info.(rname),'i1.*i2>0.1');
                cleanup(rname);
            end

            if strcmp(rname(1:2),'GM')
                ROI_info.(rname) = spm_imcalc([ROI_info.GM temp_info],...
                    ROI_info.(rname),'i1.*i2>0.1');
                cleanup(rname);
            end
        end % of loop over ROIs
        
    end % if no incomplete run
        
    %% populate masks after cleanup
    % loop over ROIs
    for nr = 1:numel(ROI.name)
        % get ROI name
        rname = ROI.name{nr};

        % create SPM handles
        ROI_info.(rname) = spm_vol(fullfile(outdir,strcat(rname,'.nii')));
        
        % load data and count voxels
        ROI_data.(rname) = spm_read_vols(ROI_info.(rname));
        ROI_voxels.(rname) = sum(ROI_data.(rname)(ROI_data.(rname)>0));
    end
            
    
    %% calculate values
    %loop over maps and ROIs
    for cm = 1:numel(MPM.name)
        mpname = MPM.name{cm};
        for cn = 1:numel(ROI.name)
            
            rname =  ROI.name{cn};
            
            fprintf('%s %s\n',mpname,rname);
            
            % initialise temporary data
            temp_info = vol_info{1,1};
            temp_info.fname = fullfile(outdir,'temp.nii');
            
            % masking current map with current ROI
            temp_info = spm_imcalc([vol_info{pos.(mpname),1} ROI_info.(rname)], ...
                temp_info,'i1.*i2'); %#ok<NASGU>
            temp_vol = spm_read_vols(spm_vol('temp.nii')); % spm... instead of temp_info
            
%             % write current mask for checking
%             temp_info.fname = char(strcat(pnames(pn,:),'/',mpname,'_',rname,'.nii'));
%             spm_write_vol(temp_info,temp_vol);
            
            % calculate values (mean, standard deviation, Coefficient of
            % Variance and afterwards [outer loop] Contrast to Noise Ratio)
            temp_vol = nonzeros(temp_vol); % get rid of vovels = 0
            temp_vol(isnan(temp_vol)) = []; % get rid of voxels = NaN

            if cn == pos.wb % whole brain mask -> only draw histogram, skip rest
                h = figure('Position',[50 50 600 400]);
                set(gcf,'Visible', 'off');
                hold on;
                ax1 = subplot(1,1,1);
                histogram(ax1,temp_vol);
                saveas(h,char(strcat(rname,'_',strcat(MPM.name{cm}),'_hist.png')));
                close(h);
              
               continue
            end
            
            mean_value.(mpname).(rname) = mean(temp_vol(:));
            med_value.(mpname).(rname) =  median(temp_vol(:));
            std_value.(mpname).(rname) = std(temp_vol(:));
            per05.(mpname).(rname) = prctile(temp_vol(:),.05);
            per95.(mpname).(rname) = prctile(temp_vol(:),.95);
            kurt.(mpname).(rname) =  kurtosis(temp_vol(:));
            cov_value.(mpname).(rname) = std_value.(mpname).(rname)/...
                mean_value.(mpname).(rname);
            
            %% calculate difference to literature
            switch cn
                case {1,6,8,9,10,12,14}
                    cp = 1;
                case {2,5,7,11,13,15,16}
                    cp = 3;
                case 3
                    cp = 2;
                case 4
                    cp = 4;
            end

            if cn < 5
            diff_mean.(mpname).(rname) = (mean_value.(mpname).(rname) - ...
                mean_lit.(mpname(1:2))(cp))/mean_lit.(mpname(1:2))(cp);
            diff_cov.(mpname).(rname) = (cov_value.(mpname).(rname) - ...
                cov_lit.(mpname(1:2))(cp))/cov_lit.(mpname(1:2))(cp);
            end

            %% create histogram of data and write to file
            try
                h=figure; set(gcf,'Visible', 'off'); hold on;
                ax1 = subplot(1,1,1);
                histfit(temp_vol(:)); 
                l1 = line([mean_value.(mpname).(rname) mean_value.(mpname).(rname)],ax1.YAxis.Limits,'Color','r');
                if cn < 5
                l2 = line([mean_lit.(mpname(1:2))(cp) mean_lit.(mpname(1:2))(cp)],ax1.YAxis.Limits,'Color','k');
                legend(ax1,[l1,l2],{char(strcat(rname,' mean')),...
                'THS'},'Location','southoutside',...
                'Orientation','horizontal','Interpreter','none');
                end
                saveas(h,char(strcat(rname,'_',strcat(MPM.name{cm}),'_histfit.png')));
                h1 = histogram(ax1,temp_vol(:));
                hist.Values = h1.Values; 
                hist.BinEdges = h1.BinEdges;
                hist.BinWidth = h1.BinWidth;
                hist.BinLimits = h1.BinLimits;
                close(h);
            catch error  %#ok<NASGU>
                %fprintf('no histfit for %s\n',char(strcat(ROI.name(cn),'_',strcat(MPM.name{cm}))));
                h = figure('Position',[50 50 600 400]);
                set(gcf,'Visible', 'off');
                hold on;

                ax1 = subplot(1,1,1);
                h1 = histogram(ax1,temp_vol(:));
                hist.Values = h1.Values; 
                hist.BinEdges = h1.BinEdges;
                hist.BinWidth = h1.BinWidth;
                hist.BinLimits = h1.BinLimits;
                l1 = line([mean_value.(mpname).(rname) mean_value.(mpname).(rname)],ax1.YAxis.Limits,'Color','r');
                if cn < 5
                l2 = line([mean_lit.(mpname(1:2))(cp) mean_lit.(mpname(1:2))(cp)],ax1.YAxis.Limits,'Color','k');
                legend(ax1,[l1,l2],{char(strcat(rname,' mean')),...
                    'THS'},'Location','southoutside',...
                    'Orientation','horizontal','Interpreter','none');
                end
                saveas(h,char(strcat(rname,'_',strcat(MPM.name{cm}),'_hist.png')));
                close(h);
            end
            %% calculate FWHM
            % code adapted from https://www.mathworks.com/matlabcentral/answers/77032-how-to-calculate-histogram-width-at-the-half-height
            [fullHeight, indexOfMax] = max(hist.Values);
            halfHeight = fullHeight / 2;
            % Initialize
            index1 = indexOfMax;
            index2 = indexOfMax;
            % Search dark side until values fall below half height.
            for k = indexOfMax-1 : -1 : 1
                if hist.Values(k) < halfHeight
                    break;
                end
                index1 = k;
            end
            % Search bright side until values fall below half height.
            for k = indexOfMax+1 : numel(hist.Values)
                if hist.Values(k) < halfHeight
                    break;
                end
                index2 = k;
            end
            FWHM.(mpname).(rname) = hist.BinEdges(index2)-hist.BinEdges(index1);
        end % of loop over ROIs
        
        % full GM/WM CNR
        cnr_value.(mpname).GMWM = abs(mean_value.(mpname).GM - ...
            mean_value.(mpname).WM)/sqrt(std_value.(mpname).GM^2 + ...
            std_value.(mpname).WM^2);
    end    
    
    %% create figure of first ROI values
    h = figure('Position',[50 50 1600 800]);
    set(gcf,'Visible', 'off');
    hold on;
 
    ax1 = subplot(2,5,1);
    ax2 = subplot(2,5,2);
    ax3 = subplot(2,5,3);
    ax4 = subplot(2,5,4);
    ax5 = subplot(2,5,5);
    ax6 = subplot(2,5,6);
    ax7 = subplot(2,5,7);
    ax8 = subplot(2,5,8);
    ax9 = subplot(2,5,9);
    ax10 = subplot(2,5,10);
        
    mv1 = mean_values('R1');
    bar(ax1,[mv1,mean_lit.R1']);
    ax1.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    legend(ax1,{'R1', 'THS'},'Location','southoutside','Orientation','horizontal');
    dm1 = diff_means('R1'); 
    for i=1:4
        if abs(dm1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax1,i-.2,mv1(i)*1.02,num2str(dm1(i)*100,ft),'rotation', 90);
    end
         
    mv1 = mean_values('R2s_OLS');
    bar(ax2,[mv1,mean_lit.R2']);
    ax2.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    legend(ax2,{'R2*', 'THS'},'Location','southoutside','Orientation','horizontal');
	dm1 = diff_means('R2s_OLS'); 
    for i=1:4
        if abs(dm1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax2,i-.2,mv1(i)*1.02,num2str(dm1(i)*100,ft),'rotation', 90);
    end
    
    mv = mean_values('MT');
    label = 'MTsat';
    bar(ax3,[mv,mean_lit.MT']);
    ax3.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dm = diff_means('MT');
    for i=1:4
        if abs(dm(i)*100)>=99, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax3,i-.2,mv(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
    end
    legend(ax3,{label,'THS'},'Location','southoutside','Orientation','horizontal');

     mv = mean_values('PD');
     bar(ax4,[mv,mean_lit.PD']);
     ax4.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
     legend(ax4,{'PD','THS'},'Location','southoutside','Orientation','horizontal');
     dm = diff_means('PD');
     for i=1:4
         if abs(dm(i)*100)>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
         text(ax4,i-.2,mv(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
     end
     
     mv = mean_values('T1w');
     bar(ax5,[mv,mean_lit.T1']);
     ax5.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
     legend(ax5,{'T1','THS'},'Location','southoutside','Orientation','horizontal');
     dm = diff_means('T1w');
     for i=1:4
         if abs(dm(i)*100)>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
         text(ax5,i-.2,mv(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
     end
     
    cv1 = cov_values('R1');
    bar(ax6,[cv1,cov_lit.R1'].*100);
    ax6.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    ax6.YAxis.TickLabelFormat = '%g%%';
    dc1 = diff_covs('R1');
    for i=1:4
        if abs(dc1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax6,i-.2,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
    end
    legend(ax6,{'CoV R1','THS'},'Location','southoutside','Orientation','horizontal');
        
    cv1 = cov_values('R2s_OLS');
    bar(ax7,[cv1,cov_lit.R2'].*100);
    ax7.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    ax7.YAxis.TickLabelFormat = '%g%%';
    dc1 = diff_covs('R2s_OLS');
    for i=1:4
        if abs(dc1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax7,i-.2,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
    end
    legend(ax7,{'CoV R2*','THS'},'Location','southoutside','Orientation','horizontal');
    
    cv1 = cov_values('MT');
    b=bar(ax8,[cv1,cov_lit.MT'].*100);
    ax8.YAxis.TickLabelFormat = '%g%%';
    ax8.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dc1 = diff_covs('MT');
    for i=1:4
        if abs(dc1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax8,i-.2,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
    end
    legend(ax8,{'CoV MT','THS'},'Location','southoutside','Orientation','horizontal');
    
    cv1 = cov_values('PD');
    b=bar(ax9,[cv1,cov_lit.PD'].*100);
    ax9.YAxis.TickLabelFormat = '%g%%';
    ax9.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dc1 = diff_covs('PD');
    for i=1:4
        if abs(dc1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax9,i-.2,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
    end
    legend(ax9,{'CoV PD','THS'},'Location','southoutside','Orientation','horizontal');
    
    cv1 = cov_values('T1w');
    b=bar(ax10,[cv1,cov_lit.T1'].*100);
    ax10.YAxis.TickLabelFormat = '%g%%';
    ax10.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dc1 = diff_covs('T1w');
    for i=1:4
        if abs(dc1(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax10,i-.2,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
    end
    legend(ax10,{'CoV T1','THS'},'Location','southoutside','Orientation','horizontal');
 
%     cnr1 = cnr_values('GMWM');
%     bar(ax8,cnr1);
%     ax8.XTickLabel = {'R1','PD','R2*','R2* E','MT' 'T1w'};
%     legend(ax8,'wholebrain GM/WM CNR','Location','southoutside','Orientation','horizontal');
    
    saveas(h,'MPMqlt_full.png');
    close(h);

    %% write data to file
    measures = struct('ROI_names',struct('ROI_names',ROI.name), ...
        'ROI_voxels',ROI_voxels,'mean',mean_value, ...
        'median',med_value,'std',std_value,'cov',cov_value,'CNR',cnr_value, ...
       'FWHM',FWHM,'p05',per05,'p95',per95,'Kurtosis',kurt); %#ok<NASGU>
    save('measures_ext_v4.mat','measures','-v7.3');
    
  
    % cleanup
    delete(fullfile(outdir,'temp.nii'));
    %delete(fullfile(outdir,strcat('y_i',bn,'_MT.nii')));

end    

function cleanup(rname)
    delete(strcat(rname,'.nii'),strcat('inv',rname,'.nii'));
    movefile(strcat('rinv',rname,'.nii'),strcat(rname,'.nii'));
end

function mv = mean_values(mpname)
    mv = [mean_value.(mpname).GM;mean_value.(mpname).GM_CN; ...
        mean_value.(mpname).WM;mean_value.(mpname).WM_CC];
end

function cv = cov_values(mpname)
    cv = [cov_value.(mpname).GM;cov_value.(mpname).GM_CN; ...
        cov_value.(mpname).WM;cov_value.(mpname).WM_CC];
end

function dm = diff_means(mpname)
    dm = [diff_mean.(mpname).GM;diff_mean.(mpname).GM_CN; ...
        diff_mean.(mpname).WM;diff_mean.(mpname).WM_CC];
end

function dc = diff_covs(mpname)
    dc = [diff_cov.(mpname).GM;diff_cov.(mpname).GM_CN; ...
        diff_cov.(mpname).WM;diff_cov.(mpname).WM_CC];
end

function cv = cnr_values(ROI_name)
    cv = [cnr_value.('R1').(ROI_name);cnr_value.('PD').(ROI_name); ...
        cnr_value.('R2s').(ROI_name);cnr_value.('R2s_OLS').(ROI_name); ...
        cnr_value.('MT').(ROI_name);cnr_value.('T1w').(ROI_name)];
end

function idefcheck % to check/create inverse deformation field
% check, if inverse deformation field exists, otherwise
% create one...
if (mn == pos.idef) && ~exist(char(mnames{mn}),'file')
    seg8file = strrep(char(mnames{mn}),'.nii','_seg8.mat');
    seg8file = strrep(seg8file,'iy_','');
    create_seg_from_mat(seg8file)
end
end

end
