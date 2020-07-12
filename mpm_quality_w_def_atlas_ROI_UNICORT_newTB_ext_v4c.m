function mpm_quality_w_def_atlas_ROI_UNICORT_newTB_ext_v4b(pnames)

% This function calculates some important measures and compares them to 
% values from the publication "Quantitative multi-parameter mapping of R1, 
% PD*, MT, and R2* at 3T: a multi-center validation" by N. Weiskopf et
% al. (Front. Neurosci., 10 June 2013, doi: 10.3389/fnins.2013.00095).
% Tobias Leutritz, 18.08.2016
% 2018-08-23: add switch from 2018-06-25 for OLS fit of averages + use of segmentations from MPMCalc
% 2018-11-16: adapt to equivalent B1maps version

% get files to be analyzed (segmented data from subjects, each in a folder)
%pnames = cfg_getfile(Inf,'dir','Select folders with segmented MPM images to analyse');
% get number of selected path 
[npth,~] = size(pnames);

% Define maps to be read (including two {normalised} segmentations [c1, c2,
% and c3] and the inverse deformation field [iy_*] to unwarp the ROI masks 
% from normalised space back to the maps space in order to avoid manual ROI 
% definition for each subject.
map.name_suffix = {'PD' 'T1w' 'R1' 'R1_UNICORT' 'R2s' 'R2s_OLS' 'MTsat' 'MTsat' 'MTsat' 'MTsat' 'MTsat'}; % _OLSfit_TEzero might be added to T1w
map.name_prefix = {'' '' '' '' '' '' '' 'c1' 'c2' 'c3' 'iy_'};
%                                        8     9   10    11
map.subfolder = {'' 'Supplementary' 'Supplementary' '' 'Supplementary' '' '' '' '' '' ''};
MPM.name = {'PD' 'T1w' 'R1' 'R1_UNICORT' 'R2s' 'R2s_ESTATICS' 'MT'};

[cp,~,~] = fileparts(mfilename('fullpath'));
atlas_ROI_dir = fullfile(cp,'ROIs');
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
mc = numel(map.name_suffix);
parf = 'hMRI_map_creation_mpm_params.json';

% define ROIs and MPM names for evaluation
ROI.name = {'GM' 'WM' 'GM_CN'  'WM_CC' ...
            'WM_CST' 'GM_S1' 'GM_M1' ...
            'GM_cerebellum' 'WM_cerebellum' ...
            'GM_Hippocampi' 'wholebrain'};
pos.wb = numel(ROI.name);


% literature values (for reference see preamble)
% order of values: GM, CN, WM, CC
mean_lit.R1 = [609 683 1036 1158]./1000;
mean_lit.R2 = [0.0152 0.0182 0.021 0.025].*1000;
mean_lit.MT = [794 836 1764 1978]./1000;
mean_lit.T1 = [319 321 398 414];
mean_lit.PD = [84.44 82.67 68.35 64.65];
mean_lit.A_ = [84.44 82.67 68.35 64.65];
cov_lit.R1 = [6.0 4.7 4.6 4.6]./100;
cov_lit.R2 = [20.3 12.0 11.4 12.1]./100;
cov_lit.MT = [7.9 7.6 6.1 7.4]./100;
cov_lit.T1 = [15.2 13.0 15.1 14.7]./100;
cov_lit.PD = [3.6 2.7 2.7 2.7]./100;
cov_lit.A_ = [3.6 2.7 2.7 2.7]./100;
sd_lit.R1 = [8 22 36 50];
sd_lit.R2 = [.4 1.2 .8 .5];
sd_lit.MT = [14 27 66 85]./1000;
sd_lit.T1 = [24 15 15 35];
sd_lit.PD = [1.87 1.64 0.06 0.86];
sd_lit.A_ = [1.87 1.64 0.06 0.86];

% initialise cell variables
vol_info = cell(mc,1);
%vol_data = cell(mc,1);
mnames = cell(mc,1);

%% loop over given path
for pn = 1:npth
    % change working directory to current subject folder
    cd(char(pnames(pn,:)));
    curdir = cd;
    outdir = fullfile(curdir,'MPM_Quality3');
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
    for mn = 1:mc %(mc-1) reuse def from former run
        if mn == pos.T1w % T1w
            mpm_param = get_metadata(fullfile(curdir,'Supplementary',parf));
            if mpm_param{1,1}.fullOLS
                map.name_suffix{mn} = 'T1w_OLSfit_TEzero';
            else
                map.name_suffix{mn} = 'T1w';
            end
        end
        if mn > 7
            if ~isempty(dir('../MPMCalc')) % reuse files from processing
                % create file names of maps
                mnames{mn} = fullfile(strrep(curdir,'Results','MPMCalc'), ...
                    strcat(map.name_prefix{mn},bn,'_MTsat_outer_suppressed.nii'));
                % check, if inverse defronamtion field exists, otherwise
                % create one...
                if (mn == pos.idef) && ~exist(char(mnames{mn}),'file')
                    seg8file = strrep(char(mnames{mn}),'.nii','_seg8.mat');
                    seg8file = strrep(seg8file,'iy_','');
                    create_seg_from_mat(seg8file)
                end
            end
        else
            % create file names of maps
            mnames{mn} = fullfile(curdir,map.subfolder{mn}, ...
                strcat(map.name_prefix{mn},bn,'_',map.name_suffix{mn},'.nii'));
        end
        % create SPM handles
        vol_info{mn,1} = spm_vol(mnames{mn});
    end
    
    % switch to output directory
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
%         inv.inv.space = {vol_info{7,1}.fname};
%         job.comp = {inv};
%         
%         % run spm_deformations and set reference to output
%         spm_deformations(job);
%         vol_info{pos.idef,1} = spm_vol(fullfile(outdir,strcat('y_i',bn,'_MT.nii')));
        
        % loop over ROIs in order to calculate or load predefiened atlas ROIs
        for nr = 3:numel(ROI.name)-1 % skip GM/WM/WB
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
            job.woptions.vox = [abs(P(7:9))];
            job.woptions.bb = spm_get_bbox(char(vol_info{1,1}.fname), 'fv');
            job.woptions.interp = 4;
            job.woptions.prefix = 'inv';
            job.subj.resample = {ROI_info.(rname).fname};
            job.subj.def = {vol_info{11,1}.fname};
            
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
            
            % initialise temporary data
            temp_info = vol_info{1,1};
            temp_info.fname = fullfile(outdir,'temp.nii');
            
            % masking current map with current ROI
            temp_info = spm_imcalc([vol_info{cm,1} ROI_info.(rname)], ...
                temp_info,'i1.*i2');
            temp_vol = spm_read_vols(spm_vol('temp.nii')); % spm... instead of temp_info
            
%             % write current mask for checking (for debugging)
%             temp_info.fname = char(strcat(pnames(pn,:),'/',mpname,'_',rname,'.nii'));
%             spm_write_vol(temp_info,temp_vol);
            
            % calculate values (mean, standard deviation, Coefficient of
            % Variance and afterwards [outer loop] Contrast to Noise Ratio)
            temp_vol = nonzeros(temp_vol); % get rid of vovels = 0
            temp_vol(isnan(temp_vol)) = []; % get rid of voxels = NaN
            mean_value.(mpname).(rname) = mean2(temp_vol);
            med_value.(mpname).(rname) =  median(temp_vol(:));
            std_value.(mpname).(rname) = std2(temp_vol);
            per05.(mpname).(rname) = prctile(temp_vol(:),.05);
            per95.(mpname).(rname) = prctile(temp_vol(:),.95);
            kurt.(mpname).(rname) =  kurtosis(temp_vol(:));
            cov_value.(mpname).(rname) = std_value.(mpname).(rname)/...
                mean_value.(mpname).(rname);
            
            % calculate difference to literature
            switch cn
                case {1,6,7,8,10}
                    cp = 1;
                case {2,5,9}
                    cp = 3;
                case 3
                    cp = 2;
                case 4
                    cp = 4;
            end
                
            diff_mean.(mpname).(rname) = (mean_value.(mpname).(rname) - ...
                mean_lit.(mpname(1:2))(cp))/mean_lit.(mpname(1:2))(cp);
            diff_cov.(mpname).(rname) = (cov_value.(mpname).(rname) - ...
                cov_lit.(mpname(1:2))(cp))/cov_lit.(mpname(1:2))(cp);
            
            %% create histogram of data and write to file
            try
                h=figure; set(gcf,'Visible', 'off'); hold on;
                ax1 = subplot(1,1,1);
                histfit(temp_vol(:)); 
                l1 = line([mean_value.(mpname).(rname) mean_value.(mpname).(rname)],ax1.YAxis.Limits,'Color','r');
                l2 = line([mean_lit.(mpname(1:2))(cp) mean_lit.(mpname(1:2))(cp)],ax1.YAxis.Limits,'Color','k');
                legend(ax1,[l1,l2],{char(strcat(rname,' mean')),...
                'literature'},'Location','southoutside',...
                'Orientation','horizontal','Interpreter','none');
                saveas(h,char(strcat(rname,'_',strcat(MPM.name{cm}),'_histfit.png')));
                h1 = histogram(ax1,temp_vol(:));
                hist.Values = h1.Values; 
                hist.BinEdges = h1.BinEdges;
                hist.BinWidth = h1.BinWidth;
                hist.BinLimits = h1.BinLimits;
                close(h);
            catch error 
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
                l2 = line([mean_lit.(mpname(1:2))(cp) mean_lit.(mpname(1:2))(cp)],ax1.YAxis.Limits,'Color','k');
                legend(ax1,[l1,l2],{char(strcat(rname,' mean')),...
                    'literature'},'Location','southoutside',...
                    'Orientation','horizontal','Interpreter','none');
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
    
%     %% create figure of regional values
%     
%     h = figure('Position',[50 50 1600 800]);
%     set(gcf,'Visible', 'off');
%     hold on;
%  
%     ax1 = subplot(1,2,1);
%     ax2 = subplot(1,2,2);
% 
%     bar(ax1,[mean_value('R1')',mean_value('R1_UNICORT')',mean_lit.R1']);
%     ax1.XTickLabel = ROI.name;
%     legend(ax1,'R1','R1 UNICORT', 'literature','Location','southoutside','Orientation','horizontal');
%     dm1 = diff_means('R1'); 
%     dm2 = diff_means('R1_UNICORT'); 
%     for i=1:4
%         if dm1(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax1,i-.3,mv1(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
%         if dm2(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax1,i,mv2(i)*1.02,num2str(dm2(i)*100,ft),'rotation', 90);
%     end
%          
%     mv1 = mean_values('R2s');
%     mv2 = mean_values('R2s_ESTATICS');
%     bar(ax2,[mv1,mv2,mean_lit.R2']);
%     ax2.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     legend(ax2,'R2*','R2* ESTATICS', 'literature','Location','southoutside','Orientation','horizontal');
% 	dm1 = diff_means('R2s'); 
%     dm2 = diff_means('R2s_ESTATICS'); 
%     for i=1:4
%         if dm1(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax2,i-.3,mv1(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
%         if dm2(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax2,i,mv2(i)*1.02,num2str(dm2(i)*100,ft),'rotation', 90);
%     end
%     
%     mv = mean_values('MT');
%     scale = floor(max(mean_lit.MT) / max(mv));
%     label = strcat('MT (x',int2str(scale),')');
%     bar(ax3,[mv.*scale,mean_lit.MT']);
%     ax3.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     dm = diff_means('MT');
%     for i=1:4
%         if abs(dm(i)*100)>=99 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax3,i-.2,mv(i)*scale*1.02,num2str(dm(i)*100,ft),'rotation', 90);
%     end
%     legend(ax3,label,'literature','Location','southoutside','Orientation','horizontal');
% 
%     mv = mean_values('T1w');
%     bar(ax4,[mv,mean_lit.T1]);
%     ax4.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     legend(ax4,'T1w','literature','Location','southoutside','Orientation','horizontal');
%     dm = diff_means('T1w');
%     for i=1:4
%         if abs(dm(i)*100)>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax4,i-.2,mv(i)*1.02,num2str(dm(i)*100,ft),'rotation', 90);
%     end
%     
%     cv1 = cov_values('R1');
%     cv2 = cov_values('R1_UNICORT');
%     bar(ax5,[cv1,cv2,cov_lit.R1'].*100);
%     ax5.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     ax5.YAxis.TickLabelFormat = '%g%%';
%     dc1 = diff_covs('R1');
%     dc2 = diff_covs('R1_UNICORT');
%     for i=1:4
%         if dc1(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax5,i-.3,cv1*102,num2str(dc1(i)*100,ft),'rotation', 90);
%         if dc2(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax5,i,cv(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
%     end
%     legend(ax5,'CoV R1','CoV R1 U','literature','Location','southoutside','Orientation','horizontal');
%         
%     cv1 = cov_values('R2s');
%     cv2 = cov_values('R2s_ESTATICS');
%     bar(ax6,[cv1,cv2,cov_lit.R2'].*100);
%     ax6.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     ax6.YAxis.TickLabelFormat = '%g%%';
%     dc1 = diff_covs('R2s');
%     dc2 = diff_covs('R2s_ESTATICS');
%     for i=1:4
%         if dc1(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax6,i-.3,cv1(i)*102,num2str(cv1(i)*100,ft),'rotation', 90);
%         if dc2(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax6,i,cv2(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
%     end
%     legend(ax6,'CoV R2*','CoV R2* E','literature','Location','southoutside','Orientation','horizontal');
%     
%     cv1 = cov_vlaues('MT');
%     cv2 = cov_values('T1w');
%     b=bar(ax7,[cv1,cov_lit.MT',cv2',cov_lit.T1'].*100);
%     b(2).FaceColor = 'y';
%     ax7.YAxis.TickLabelFormat = '%g%%';
%     ax7.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
%     dc1 = diff_covs('MT');
%     dc2 = diff_covs('T1w');
%     for i=1:4
%         if dc1(i)*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax7,i-.3,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
%         if diff_cov.T1w(i)>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
%         text(ax7,i+.1,cv2(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
%     end
%     legend(ax7,[b(1),b(3),b(4)],'CoV MT','CoV T1w','literature','Location','southoutside','Orientation','horizontal');
%     
%     cv = cnr_values(1);
%     bar(ax8,cv);
%     ax8.XTickLabel = {'R1','R1 U','R2*','R2* E','MT','T1w'};
%     legend(ax8,'GM/WM CNR','Location','southoutside','Orientation','horizontal');
%     
%     saveas(h,'MPMqlt_regions.png');
%     close(h);
    
    %% create figure of other values
    h = figure('Position',[50 50 1600 800]);
    set(gcf,'Visible', 'off');
    hold on;
 
    ax1 = subplot(2,4,1);
    ax2 = subplot(2,4,2);
    ax3 = subplot(2,4,3);
    ax4 = subplot(2,4,4);
    ax5 = subplot(2,4,5);
    ax6 = subplot(2,4,6);
    ax7 = subplot(2,4,7);
    ax8 = subplot(2,4,8);
        
    mv1 = mean_values('R1');
    mv2 = mean_values('R1_UNICORT');
    bar(ax1,[mv1,mv2,mean_lit.R1']);
    ax1.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    legend(ax1,'R1','R1 UNICORT', 'literature','Location','southoutside','Orientation','horizontal');
    dm1 = diff_means('R1'); 
    dm2 = diff_means('R1_UNICORT'); 
    for i=1:4
        if abs(dm1(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax1,i-.3,mv1(i)*1.02,num2str(dm1(i)*100,ft),'rotation', 90);
        if abs(dm2(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax1,i,mv2(i)*1.02,num2str(dm2(i)*100,ft),'rotation', 90);
    end
         
    mv1 = mean_values('R2s');
    mv2 = mean_values('R2s_ESTATICS');
    bar(ax2,[mv1,mv2,mean_lit.R2']);
    ax2.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    legend(ax2,'R2*','R2* ESTATICS', 'literature','Location','southoutside','Orientation','horizontal');
	dm1 = diff_means('R2s'); 
    dm2 = diff_means('R2s_ESTATICS'); 
    for i=1:4
        if abs(dm1(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax2,i-.3,mv1(i)*1.02,num2str(dm1(i)*100,ft),'rotation', 90);
        if abs(dm2(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax2,i,mv2(i)*1.02,num2str(dm2(i)*100,ft),'rotation', 90);
    end
    
    mv = mean_values('MT');
    scale = floor(max(mean_lit.MT) / max(mv));
    if scale == 0, scale = 1; end
    label = strcat('MT (x',int2str(scale),')');
    bar(ax3,[mv.*scale,mean_lit.MT']);
    ax3.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dm = diff_means('MT');
    for i=1:4
        if abs(dm(i)*100)>=99 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax3,i-.2,mv(i)*scale*1.02,num2str(dm(i)*100,ft),'rotation', 90);
    end
    legend(ax3,label,'literature','Location','southoutside','Orientation','horizontal');

    mv1 = mean_values('T1w');
    mv2 = mean_values('PD');
    bar(ax4,[mv1,mean_lit.T1',mv2,mean_lit.PD']);
    ax4.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    legend(ax4,{'T1w','lit','PD','lit'},'Location','southoutside','Orientation','horizontal');
    dm1 = diff_means('T1w');
    dm2 = diff_means('PD');
    for i=1:4
        if abs(dm1(i)*100)>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax4,i-.2,mv1(i)*1.02,num2str(dm1(i)*100,ft),'rotation', 90);
        if abs(dm2(i))*100>=100, ft = '%+.0f%%'; else, ft = '%+.2g%%'; end
        text(ax4,i,mv2(i)*102,num2str(dm2(i)*100,ft),'rotation', 90);
    end
    
    cv1 = cov_values('R1');
    cv2 = cov_values('R1_UNICORT');
    bar(ax5,[cv1,cv2,cov_lit.R1'].*100);
    ax5.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    ax5.YAxis.TickLabelFormat = '%g%%';
    dc1 = diff_covs('R1');
    dc2 = diff_covs('R1_UNICORT');
    for i=1:4
        if abs(dc1(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax5,i-.3,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
        if abs(dc2(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax5,i,cv2(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
    end
    legend(ax5,'CoV R1','CoV R1 U','literature','Location','southoutside','Orientation','horizontal');
        
    cv1 = cov_values('R2s');
    cv2 = cov_values('R2s_ESTATICS');
    bar(ax6,[cv1,cv2,cov_lit.R2'].*100);
    ax6.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    ax6.YAxis.TickLabelFormat = '%g%%';
    dc1 = diff_covs('R2s');
    dc2 = diff_covs('R2s_ESTATICS');
    for i=1:4
        if abs(dc1(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax6,i-.3,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
        if abs(dc2(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax6,i,cv2(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
    end
    legend(ax6,'CoV R2*','CoV R2* E','literature','Location','southoutside','Orientation','horizontal');
    
    cv1 = cov_values('MT');
    cv2 = cov_values('T1w');
    b=bar(ax7,[cv1,cov_lit.MT',cv2,cov_lit.T1'].*100);
    b(2).FaceColor = 'y';
    ax7.YAxis.TickLabelFormat = '%g%%';
    ax7.XTickLabel = {'GM' 'CN' 'WM' 'CC'};
    dc1 = diff_covs('MT');
    dc2 = diff_covs('T1w');
    for i=1:4
        if abs(dc1(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax7,i-.3,cv1(i)*102,num2str(dc1(i)*100,ft),'rotation', 90);
        if abs(dc2(i))*100>=100 ft = '%+.0f%%'; else ft = '%+.2g%%'; end
        text(ax7,i+.1,cv2(i)*102,num2str(dc2(i)*100,ft),'rotation', 90);
    end
    legend(ax7,[b(1),b(3),b(4)],'CoV MT','CoV T1w','literature','Location','southoutside','Orientation','horizontal');
    
    cnr1 = cnr_values('GMWM');    
    bar(ax8,cnr1);
    ax8.XTickLabel = {'R1','R1 U','R2*','R2* E','MT','T1w','PD'};
    legend(ax8,'GM/WM CNR','Location','southoutside','Orientation','horizontal');
    
    saveas(h,'MPMqlt_full.png');
    close(h);

    %% write data to file
    measures = struct('ROI_names',struct('ROI_names',ROI.name),'ROI_voxels',ROI_voxels,...
        'mean',mean_value,'median',med_value,'std',std_value,'cov',cov_value,...
        'CNR',cnr_value,'FWHM',FWHM,'p05',per05,'p95',per95,'Kurtosis',kurt);
    save('measures_ext_v4.mat','measures','-v7.3');
    
    % cleanup
    delete(fullfile(outdir,'temp.nii'));
%     compr = gzip('*.nii');
%     compr = strrep(compr,'.gz','');
%     for ii = 1:numel(compr)
%         delete(fullfile(outdir,compr{ii}));
%     end
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

function cv = cnr_values(nr)
    cv = [cnr_value.('R1').(nr);cnr_value.('R1_UNICORT').(nr); ...
        cnr_value.('R2s').(nr);cnr_value.('R2s_ESTATICS').(nr); ...
        cnr_value.('MT').(nr);cnr_value.('T1w').(nr);cnr_value.('PD').(nr)];
end


end

