function mpm_quality_multiecho_data_atlas_ROI_v4(pnames)

% This function calculates some important measures from the input files for
% MPM maps, i.e. the multi-echo data of T1-, PD- and MT-weighted contrasts.
% MPM_quality script had to be run first in order to create ROIs in
% individual space. Therefore the variable mpm_qlt_name has to be adapted,
% according to the folder containing these MPM ROIs (and results).
% Tobias Leutritz, 16.11.2018
% extended by measures summarized in Esteban et al: MRIQC — The MRI quality 
% control tool. PLOS one 2017, doi:10.1371/journal.pone.0184661
% 2019-08-01: fix CNR calc., which was only done for last echo

mpm_qlt_name = 'MPM_Quality4';
med_qlt_name = 'MultiEchoData_Qlt';
% flags.dtype=16;%'float32';

% get files to be analyzed (segmented data from subjects, each in a folder)
 if nargin < 1
     pnames = cfg_getfile(Inf,'dir','Select Results folders with MPM Quality already run at maps');
 end

% get number of selected path 
[npth,~] = size(pnames);

% define ROIs for evaluation
ROI.name = {'GM' 'WM' 'GM_CN'  'WM_CC' ...
            'WM_CST' 'GM_S1' 'GM_M1' ...
            'GM_cerebellum' 'WM_cerebellum' ...
            'GM_Hippocampi' 'wholebrain'};
mask_name = {'headmask' 'airmask' 'bgghostfree' 'fgghostfree' 'bgghost' 'fgghost'};
% whereas bg = background, fg = foreground

%% loop over given path
for pn = 1:npth
    % change working directory to current subject folder
    cd(char(pnames(pn,:)));
    curdir = cd;
    qltdir = fullfile(curdir,mpm_qlt_name);
    outdir = fullfile(curdir,med_qlt_name);
	if ~exist(outdir,'dir') 
       mkdir(outdir);  
    elseif exist(fullfile(outdir,'measures.mat'),'file') == 2
        continue % skip further calculation, since measurment was completed
    end
    
    fprintf(1,cd)
    
    % retrieve meta info from json-file
    ref_file = dir(strcat('*R2s*.nii'));
    metadata = get_metadata(ref_file.name);
    inputs = metadata{1,1}.history.procstep.procpar.input;
    if ~isfield(metadata{1,1}.history.procstep.procpar.proc.RFsenscorr,'RF_us')
        inputs = metadata{1,1}.history.input(1).history.procstep.procpar.input;
    end
    
    for nc = 1:numel(inputs) % loop over input contrasts
        
        c_tag = inputs(nc).tag; % contrast tag
        
        files.(c_tag) = inputs(nc).fnam;
        
        nr_e = numel(files.(c_tag)); % number of multi-echo data for contrast

        % initialise cell variables
        vol_info.(c_tag) = cell(nr_e,1);
        
        % loop over multi-echo data to be read in
        for nme = 1:nr_e
            % create new file name
            med_fn = fullfile(outdir,strcat(c_tag,num2str(nme),'.nii'));
            
            % copy data to output folder for reorientation etc.
            % Must strip the ',1' (at the end of the file extension '.nii,1') 
            % if the name has been collected using e.g. spm_select: [from
            % EBs get_metadata script]
            cfname = spm_file(files.(c_tag){nme},'number','');
            copyfile(cfname,med_fn);
            
            % check for headmask and create one, if needed:
            if contains(c_tag,'T1') && (nme == 1)
                headmaskfile = strrep(cfname,'.nii','_headmask.nii');
                ghostmaskfile = strrep(cfname,'.nii','_ghostmask.nii');
                % create mask if needed
                if ~exist(headmaskfile,'file')
                    create_head_mask2(cfname,0.4,10,false,true)
                end
                % assign masks
                ROI_info.headmask = spm_vol(headmaskfile);
                ROI_data.headmask = spm_read_vols(ROI_info.headmask);
                gm_max = max(ROI_data.headmask(:));
                ROI_data.headmask(ROI_data.headmask == gm_max) = 1;
                ROI_data.airmask = double(~ROI_data.headmask);
                gm_max = max(ROI_data.airmask(:));
                ROI_data.airmask(ROI_data.airmask == gm_max) = 1;
                ROI_info.ghostmask = spm_vol(ghostmaskfile);
                ROI_data.ghostmask = spm_read_vols(ROI_info.ghostmask);
                gm_max = max(ROI_data.ghostmask(:));
                gm_sc = gm_max/3; % to fetch spm rescaling on mask
                ROI_data.bgghostfree = zeros(size(ROI_data.ghostmask));
                ROI_data.bgghostfree(ROI_data.ghostmask == 1*gm_sc) = 1;
                ROI_data.fgghostfree = zeros(size(ROI_data.ghostmask));
                ROI_data.fgghostfree(ROI_data.ghostmask == 3*gm_sc) = 1;
                ROI_data.bgghost = zeros(size(ROI_data.ghostmask));
                ROI_data.bgghost(ROI_data.ghostmask == 0) = 1;
                ROI_data.fgghost = zeros(size(ROI_data.ghostmask));
                ROI_data.fgghost(ROI_data.ghostmask == 2*gm_sc) = 1;
            end
            
            % create SPM handles
            vol_info.(c_tag){nme,1} = spm_vol(med_fn);
        end
    end
    
    cd(outdir);
        
    %% populate masks
    % loop over ROIs
    for nr = 1:numel(ROI.name)
        % get ROI name
        rname = ROI.name{nr};

        % create SPM handles
        ROI_info.(rname) = spm_vol(fullfile(qltdir,strcat(rname,'.nii')));
        
        % load data and count voxels
        ROI_data.(rname) = spm_read_vols(ROI_info.(rname));
        ROI_voxels(nr) = sum(ROI_data.(rname)(ROI_data.(rname)>0)); %#ok<AGROW>
    end

    % initialise temporary data & masks
    temp_info = vol_info.MT{1,1};
    temp_info.fname = fullfile(outdir,'temp.nii');
    for mn = mask_name
        mname = char(mn);
        % initialize spm_info
        mask.(mname) = vol_info.T1{1,1};
        mask.(mname).pinfo(1) = 1; % reset scale
        mask.(mname).fname = fullfile(outdir,strcat(mname,'.nii'));
        spm_write_vol(mask.(mname),ROI_data.(mname));
    end

    %% calculate values
    %loop over contrasts, echoes and ROIs
    for nc = 1:numel(inputs) % loop over contrasts
        c_tag = inputs(nc).tag;
        for ne = 1:numel(files.(c_tag)) % loop over echoes
            for nr = 1:numel(ROI.name) % loop over ROIs
                rname =  ROI.name{nr};

                % masking current echo data with current ROI
                temp_info = spm_imcalc([vol_info.(c_tag){ne,1} ROI_info.(rname)], ...
                    temp_info,'i1.*i2');
                temp_vol = spm_read_vols(spm_vol('temp.nii')); % spm... instead of temp_info

    %             % write current mask for checking
    %             temp_info.fname = char(strcat(pnames(pn,:),'/',mpname,'_',rname,'.nii'));
    %             spm_write_vol(temp_info,temp_vol);

                % calculate values (mean, standard deviation, Coefficient of
                % Variance and afterwards [outer loop] Contrast to Noise Ratio)
                temp_vol = nonzeros(temp_vol); % get rid of vovels = 0
                temp_vol(isnan(temp_vol)) = []; % get rid of voxels = NaN
                mean_value.(rname).(c_tag)(ne) = mean2(temp_vol);
                med_value.(rname).(c_tag)(ne) = median(temp_vol(:));
                std_value.(rname).(c_tag)(ne) = std2(temp_vol);
                per05.(rname).(c_tag)(ne) = prctile(temp_vol(:),.05);
                per95.(rname).(c_tag)(ne) = prctile(temp_vol(:),.95);
                cov_value.(rname).(c_tag)(ne) = std_value.(rname).(c_tag)(ne)/...
                    mean_value.(rname).(c_tag)(ne);
                kurt.(rname).(c_tag)(ne) =  kurtosis(temp_vol(:));
                data{ne} = temp_vol(:); %#ok<AGROW>

                %% create histogram of data and write to file
                try
                    h=figure; set(gcf,'Visible', 'off'); hold on;
                    ax1 = subplot(1,1,1);
                    histfit(temp_vol(:));
                    l1 = line([mean_value.(rname).(c_tag)(ne) mean_value.(rname).(c_tag)(ne)],ax1.YAxis.Limits,'Color','b');
                    l2 = line([med_value.(rname).(c_tag)(ne) med_value.(rname).(c_tag)(ne)],ax1.YAxis.Limits,'Color','r');
                    l3 = legend(ax1,[l1,l2],' mean','median', ...
                        'Location','southoutside','Orientation','horizontal');
                    set(l3, 'Interpreter', 'none')
                    saveas(h,char(strcat(rname,'_',c_tag,'_',num2str(ne),'_histfit.png')));
                    % retrieve normal distribution fit values
                    tmp = fitdist(temp_vol(:),'Normal');
                    fit.(rname).(c_tag)(ne) = struct('mu',tmp.mu, ...
                        'sigma',tmp.sigma,'ParameterCovariance', tmp.ParameterCovariance);
                    % create histogram in order to calculate FWHM
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
                    l1 = line([mean_value.(rname).(c_tag)(ne) mean_value.(rname).(c_tag)(ne)],ax1.YAxis.Limits,'Color','b');
                    l2 = line([med_value.(rname).(c_tag)(ne) mean_value.(rname).(c_tag)(ne)],ax1.YAxis.Limits,'Color','r');
                    l3 = legend(ax1,[l1,l2],' mean','median', ...
                        'Location','southoutside','Orientation','horizontal');
                    set(l3, 'Interpreter', 'none')
                    saveas(h,char(strcat(rname,'_',c_tag,'_',num2str(ne),'_histplot.png')));
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
                FWHM.(rname).(c_tag)(ne) = hist.BinEdges(index2)-hist.BinEdges(index1);
                
                %% calculate values from head/ghostmask
                if nr > 1 
                    continue % only need to be run once per input
                else
                    for mn = mask_name
                        mname = char(mn);
                        % initialize masked data
                        mask_file = fullfile(outdir,strcat(c_tag,'_',mname,'.nii'));
%                         copyfile(mask.(mname).fname,mask_file);
%                         mask_info.(mname).fname = spm_vol(mask_file);
%                         mask_info.(mname).pinfo(1) = 0; % reset scale
                        flags.dtype = vol_info.(c_tag){ne,1}.dt(1);
                        % masking current echo data with current mask
                        spm_imcalc([vol_info.(c_tag){ne,1} mask.(mname)], ...
                            mask_file,'i1.*i2',flags);
                         mask_data.(mname) = spm_read_vols(spm_vol(mask_file));
                    end
                    % Foreground to Background Energy Ratio
                    %mask_name = {'headmask' 'airmask' 'bgghostfree' 'fgghostfree' 'bgghost' 'fgghost'};
                    FBER.(c_tag).whole(ne) = var(mask_data.headmask(:))/var(mask_data.airmask(:));
                    FBER.(c_tag).ghostfree(ne) = var(mask_data.fgghostfree(:))/var(mask_data.bgghostfree(:));
                    % Entropy Focus Criterion
                    temp_vol = spm_read_vols(vol_info.(c_tag){ne,1});
                    Bmax = sum(temp_vol(:));
                    EFC.(c_tag)(ne) = 0; % figure; hold on;
                    for n = 1:numel(temp_vol)
                        EFC.(c_tag)(ne) = EFC.(c_tag)(ne) - (abs(temp_vol(n)+1i)/Bmax)*log(abs((temp_vol(n)+1i)/Bmax));
                        %efs(n) = EFC.(c_tag)(ne);
                    end
                    %plot(1:numel(temp_vol),efs)
                    SNR.(c_tag)(ne) = mean(mask_data.fgghostfree(:))/std(mask_data.bgghostfree(:));
                    QI1.(c_tag)(ne) = numel(nonzeros(mask_data.bgghost(:)))/numel(nonzeros(mask_data.airmask(:)));
                    [h,p,stats] = chi2gof(nonzeros(mask_data.bgghost(:)));
                    QI2.(c_tag)(ne).h = h;
                    QI2.(c_tag)(ne).p = p;
                    QI2.(c_tag)(ne).stats = stats;
                end 
            end % of loop over ROIs
             
            % full GM/WM CNR
            CNR.(rname).(c_tag)(ne) = abs(mean_value.GM.(c_tag)(ne) - ...
                mean_value.WM.(c_tag)(ne))/sqrt(std_value.GM.(c_tag)(ne)^2 + ...
                std_value.WM.(c_tag)(ne)^2);
            CJV.(rname).(c_tag)(ne) =  (std_value.GM.(c_tag)(ne) + std_value.WM.(c_tag)(ne))./...
                (mean_value.WM.(c_tag)(ne) - mean_value.GM.(c_tag)(ne));
        end % of loop over echoes
        % boxplots over echoes
            h=figure; set(gcf,'Visible', 'off'); hold on;
            ax1 = subplot(1,1,1);
            boxplot(ax1,padcat(data{:}),'Notch','on','Whisker',1.5); % 'OutlierSize',0
            %xticklabels(s1,{'1','2','3','4'});
            saveas(h,char(strcat(rname,'_',c_tag,'_boxplot.png')));
            close(h);
    end % of loop over contrasts

    %% write data to file
    measures = struct('ROI_names',struct('ROI_names',ROI.name), ...
        'ROI_voxels',ROI_voxels,'inputs',struct('MT',inputs(1),'PD',inputs(2),'T1',inputs(3)), ...
        'mean',mean_value,'median',med_value,'std',std_value, ...
        'cov',cov_value,'CNR',CNR,'CJV',CJV,'SNR',SNR,'EFC',EFC,'FBER',FBER, ...
        'QI1',QI1,'QI2',QI2,'FWHM',FWHM,'p05',per05,'p95',per95,...
        'Kurtosis',kurt,'histfit',fit); %#ok<NASGU>
    save('measures.mat','measures','-v7.3');

delete(fullfile(outdir,'temp.nii'));
 for nc = 1:numel(inputs)
    c_tag = inputs(nc).tag;
    for ne = 1:numel(files.(c_tag))
        med_fn = fullfile(outdir,strcat(c_tag,num2str(ne),'.nii'));
        delete(med_fn);
    end
    for mn = mask_name
        mname = char(mn);
        temp_fn = fullfile(outdir,strcat(c_tag,'_',mname,'.nii'));
        delete(temp_fn);
        if nc > 1
            continue
        else
            mask_fn = fullfile(outdir,strcat(mname,'.nii'));
            delete(mask_fn);
        end
    end
 end

end % of loop over path
end % of function
