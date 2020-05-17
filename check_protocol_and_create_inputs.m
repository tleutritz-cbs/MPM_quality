% PURPOSE
% Script to check MRI protocol parameters and completeness of measurement
% procedure according to a set of sequences and settings within a sepearate
% file, which has to be selected with the variable 'protocol'.
% REQUIREMENTS
% Requires SPM and the hMRI-toolbox, unanonymized NIFTI files converted by
% SPM and sorted into folders with ProtocolName, alongside with JSON metadata.
% =========================================================================
% Written by Tobias Leutritz (2017-05-10)
% =========================================================================
% 2018-02-02: time optimization after checking with profsave:
%             * fixed search paths for sequence parameters (>50% of time)
%             * read in jsons and store in struct permamently during run (30%)
%             * only TE checked for different echoes - rest skipped
% 2018-04-10: made compatible with Philips data
% 2019-02-19: extend output to input structure for further processing
% to do: express tolerance at least partially in percentage?
% 2019-08-06: rework protocol to HTML output
function [inputs,vcd] = check_protocol_and_create_inputs(PathName)

%% make several variables globally available
global protocol_settings
global protocol_standard

%% spm dialog to gather the folder(s) including subfolders with converted NIFTI files
if nargin == 0
    PathName = cfg_getfile(Inf,'dir','Select folders containing NiFTI data to be checked');
end

%% protocol and other presets
protocol_list.s = {'rf_map' 'head' 'body' 'MT' 'PD' 'T1' ...
    't1_tse_sag' 't2_tse_sag' 't2_tse_tra' 'MEDIC'};
protocol_list.p = {'B1Map' 'head' 'body' 'MTw' 'PDw' 'T1w' ...
    'T1w_TSE' 'T2W_mDixon' 'T2w_TSE' 'mFFE' 'B1map120' 'B1map60'};
tolerance = 0.0001; % tolerance for checking numbers
map_folder = 'maps'; % folder for map creation output
[cp,~,~] = fileparts(mfilename('fullpath'));
protocol_standard_file = fullfile(cp,'protocol_standard.mat');

if ~exist(protocol_standard_file,'file')
    % read-in json files with standard protocol settings:
    protpath = fullfile(cp,'protocol');
    for vk = {'s' 'p'}
        switch cell2mat(vk)
            case 's'
                for seq_name = protocol_list.s
                    protocol_standard.(cell2mat(vk)).(cell2mat(seq_name)) = get_metadata(strcat(protpath,filesep,cell2mat(vk),'_',cell2mat(seq_name),'.json'));
                end
            case 'p'
                for curr_seq = 1:numel(protocol_list.p)
                    seq_name = protocol_list.p{curr_seq};
                    switch curr_seq
                        case {1 2 3 7 8 9 10 11 12}
                            protocol_standard.(cell2mat(vk)).(seq_name) = ...
                                get_metadata(strcat(protpath,filesep,cell2mat(vk),'_',seq_name,'.json'));
                        case {4 5 6}
                            for echo = 1:6
                                protocol_standard.(cell2mat(vk)).(strcat(seq_name,num2str(echo))) = ...
                                    get_metadata(strcat(protpath,filesep,cell2mat(vk),'_',seq_name,int2str(echo),'.json'));
                            end
                    end
                end
        end
    end
    save(protocol_standard_file,'protocol_standard');
else
    % read-in condensed json file with standard protocol settings:
    load(protocol_standard_file,'protocol_standard');
end

npth = size(PathName,1);
inputs = cell(8,npth);

%% loop over given folders
for pn = 1:npth
    cpn = char(PathName(pn,:));
    
    % prepare report file
    fid = fopen(fullfile(cpn,'protocol_check.htm'),'w');
    str = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"' ...
                '"http://www.w3.org/TR/html4/strict.dtd">\n' ...
           '<html>\n' ...
           '<head>\n' ...
           ' <title>NISCI Protocol check results for folder %s</title>\n' ...
           '</head>\n' ...
           '<h1>Protocol check results for folder %s</h1>\n'
          ];
    fprintf(fid,str,cpn,cpn); % cpn should be replaced with unified patient ID!
    fprintf(fid,['<h2>Table of contents</h2>\n' ...
        '<a href="#protcheck">1. Checks per sequence and parameter</a><br>\n' ...
        '<a href="#summary">2. Summary of protocol check and output</a>\n' ...
        '<h2><a id="protcheck">1. Checks per sequence and parameter</a></h2>\n']);
    
    ADir = dir(cpn); % get folder directories
    ADir = ADir([ADir.isdir]); % make sure to only have folders in the list
    AName = {ADir.name}; % convert to string array
    AName(strncmp(AName, '.', 1)) = []; % get rid of folder substitutions '.' and '..'
    
    %% preset of protocol settings, which have to be checked
    firstfile = dir(fullfile(cpn,AName{1},'*.nii'));
    firstfile = fullfile(firstfile(1).folder,firstfile(1).name);
    metadata = get_metadata(firstfile);
    vendor =  metadata{1,1}.acqpar.Manufacturer;
    vcd = lower(vendor(1));
    cvl = protocol_list.(vcd);
    nop = numel(cvl);
    
    %% check which protocol is present
    PName = string(AName); % convert to string array
    nos = numel(AName); % number of apparent sequence folders
    nms = 0; % counter for missing sequences
    all_ok = true; % initialize overall check variable
    slc = [1:6,11:12]; % counter for input creation
    % find locations of the protocol names
    for prn = 1:nop
        curr_seq_name = cvl{prn}; % short name of protocol to check for
        POccur.(curr_seq_name) = contains(lower(PName),lower(curr_seq_name));
        if contains(lower(curr_seq_name),'b1map') && (vcd == 'p')
            % go through different folders
            B1maps = contains(lower(PName),'b1map');
            for nb1 = 1:numel(PName)
                if ~B1maps(nb1) % skip non B1maps
                    continue
                end
                firstfile = dir(fullfile(cpn,AName{nb1},'*.nii'));
                firstfile = fullfile(firstfile(1).folder,firstfile(1).name);
                prep_typ = get_metadata_val(firstfile,'PrepulseType');
                if (prn == 1) && (contains(prep_typ,'NO'))
                    POccur.(curr_seq_name)(nb1) = 1;
                    break
                elseif (prn >= 11) && contains(prep_typ,'SAT') && ...
                        contains(lower(PName(nb1)),lower(curr_seq_name))
                    POccur.(curr_seq_name)(nb1) = 1;
                else
                    POccur.(curr_seq_name)(nb1) = 0;
                end
            end
        end
        if contains(lower(curr_seq_name),'t2w') && (vcd == 'p')
            % go through different folders
            T2ws = contains(lower(PName),'t2w');
            for nt2w = 1:numel(PName)
                if ~T2ws(nt2w) % skip non T2ws
                    continue
                end
                if contains(lower(PName),curr_seq_name)
                    POccur.(curr_seq_name)(nt2w) = 1;
                    break
                end
            end
        end
        if contains(lower(curr_seq_name),'t1') && (prn == 6) % make sure to have multi-echo T1
            T1ws = contains(lower(PName),'t1');
            for nt1w = 1:numel(PName)
                if ~T1ws(nt1w) % skip non T1ws
                    continue
                end
                if contains(lower(PName(nt1w)),'body') || ...
                        contains(lower(PName(nt1w)),'head') || ...
                        contains(lower(PName(nt1w)),'lowres') || ...
                        contains(lower(PName(nt1w)),'nosense')
                    POccur.(curr_seq_name)(nt1w) = 0;
                    continue
                end
                firstfile = dir(fullfile(cpn,AName{nt1w},'*.nii'));
                firstfile = fullfile(firstfile(1).folder,firstfile(1).name);
                if vcd == 'p'
                    sns = 'PulseSequenceName';
                else
                    sns = 'SequenceName';
                end
                sns = upper(get_metadata_val(firstfile,sns));
                if ~isempty(sns)
                    if ~contains(sns,'SE')
                        POccur.(curr_seq_name)(nt1w) = 1;
                    else
                        POccur.(curr_seq_name)(nt1w) = 0; % TSE/SE, i.e. clinical
                    end
                elseif get_metadata_val(firstfile,'RepetitionTime') == 18
                    POccur.(curr_seq_name)(nt1w) = 1;
                else
                    POccur.(curr_seq_name)(nt1w) = 0; 
                end                    
            end
        end
        if contains(lower(curr_seq_name),'pd') && (prn == 5) % make sure to have multi-echo PD
            PDws = contains(lower(PName),'pd');
            for pdw = 1:numel(PName)
                if ~PDws(pdw) % skip non PDws
                    continue
                end
                if contains(lower(PName(pdw)),'body') || ...
                        contains(lower(PName(pdw)),'head') || ...
                        contains(lower(PName(pdw)),'lowres') || ...
                        contains(lower(PName(pdw)),'nosense')
                    POccur.(curr_seq_name)(pdw) = 0;
                    continue
                end
            end
        end
        if contains(lower(curr_seq_name),'mt') && (prn == 4) % make sure to have multi-echo MT
            MTws = contains(lower(PName),'mt');
            for mtw = 1:numel(PName)
                if ~MTws(mtw) % skip non MTws
                    continue
                end
                if contains(lower(PName(mtw)),'body') || ...
                        contains(lower(PName(mtw)),'head') || ...
                        contains(lower(PName(mtw)),'lowres') || ...
                        contains(lower(PName(mtw)),'nosense') || ...
                        contains(lower(PName(mtw)),'medic')
                    POccur.(curr_seq_name)(mtw) = 0;
                    continue
                end
            end
        end
        if strcmp(curr_seq_name,'head')
            ign = contains(lower(PName),'localizer') + contains(lower(PName),'scout'); % ignore localizer
            for cc = 1:numel(ign)
                if ign(cc)
                    POccur.head(cc) = 0;
                end
            end
            if isempty(nonzeros(POccur.(curr_seq_name))) && vcd == 'p' % try NV16 instead for Philips data
                POccur.(curr_seq_name) = contains(lower(PName),'nv16');
            end
        end
        PNum.(curr_seq_name) = numel(nonzeros(POccur.(curr_seq_name))); % count occurences
        if ~PNum.(curr_seq_name) % count missing ones
            nms = nms + 1;
            missing{nms} = curr_seq_name; %#ok<AGROW>
            if (prn < 7) || (prn >= 11)
                slc(slc == prn) = []; % delete from sequence list counter
            end
        else
            %% check protocol parameters
            nps = 0; % counter for position of current sequence in list
            % finding all occurences
            for cix = 1:nos
                if ~isempty(POccur.(curr_seq_name)(cix)) && ~(POccur.(curr_seq_name)(cix) == 0)
                    nps = nps + 1;
                    PPos.(curr_seq_name){nps} = cix; % position of sequence in the list
                end
            end
            % going through all the appearances to check each and every p.
            if (prn > 3) && (prn < 7) % for multi-echo data:
                nps = 1; % restrict following loop to first occurence only (i.e. unfiltered data)
            end
            for cc = 1:nps
                cd(strcat(cpn,filesep,PName{PPos.(curr_seq_name){cc}}));
                % check if folder is empty
                if isempty(ls)
                    fprintf(fid,'Folder <a href="%s">%s</a> of protocol %s is empty.<br>\n', ...
                        cd,strrep(cd,[cpn filesep],''),curr_seq_name);
                    continue
                else
                    % read in NIFTI files
                    mf = dir('*.nii');
                end
                files = fullfile({mf.folder},{mf.name});
                
                fprintf(fid,'Checking sequence "%s" in folder <a href="%s">%s</a><br>\n',...
                    curr_seq_name,cd,PName{PPos.(curr_seq_name){cc}});
                for f_ind = 1:numel(files)
                    if ((prn == 1) || (prn > 6)) && (f_ind > 1) % skip multiple clinical & B1map files
                        break
                    end
                    cur_file = files{f_ind};
                    if (prn > 3) && (prn < 7) % for multi-echo data:
                        [mtch,mismtch] = pcheck(vcd,curr_seq_name,cur_file,f_ind,tolerance,fid);
                    else
                        [mtch,mismtch] = pcheck(vcd,curr_seq_name,cur_file,0,tolerance,fid);
                    end
                    if ~mismtch
                        fprintf(fid,'Checked %i parameters in total. All matched.<br><br>\n', mtch+mismtch);
                    else
                        fprintf(fid,'Checked %i parameters in total. Mismatch in %i case(s).<br><br>\n', mtch+mismtch,mismtch);
                        all_ok = false;
                    end
                end
            end
        end
    end
    %% give feedback about the presence of the protocols
    fprintf(fid,'<h2><a id="summary">2. Summary of protocol check and output</a></h2>\n');
    if PNum.(cvl{1}) && (PNum.(cvl{2}) >= 3) && (PNum.(cvl{3}) >= 3) && ...
            PNum.(cvl{4}) && PNum.(cvl{5}) && PNum.(cvl{6}) && ...
            PNum.(cvl{7}) && PNum.(cvl{8}) && PNum.(cvl{9}) && PNum.(cvl{10})
        fprintf(fid,'<font color="green">Data from all sequences present.</font><br>\n');
    elseif (PNum.(cvl{2}) < 3) && (PNum.(cvl{3}) < 3)
        fprintf(fid,'<font color="red">RF sensitivity maps not acquired for all 3 MPMs.<br>\n');
        fprintf(fid,'Switch to RF sensitivity correction with single measurement.</font><br>\n');
    end
    if nms
        for cix = 1:nms
            fprintf(fid,'<font color="red">The sequence "%s" is missing.</font><br>\n',missing{cix});
        end
    end
    if ~all_ok
        fprintf(fid,['<font color="red">There were protocol deviations, '...
            'which have to be checked - see red lines in section ' ...
            '<a href="#protcheck">Checks per sequence and parameter</a>.</font><br>\n']);
    end
    
    %% prepare inputs file
    inputs{1,pn} = {fullfile(fileparts(cpn),map_folder)};
    if (vcd == 's')
        slc(slc == 11) = [];
        slc(slc == 12) = [];
    end
    for prn = slc % loop over brain protocols to create input for map creation
        curr_seq_name = cvl{prn};
        if (prn < 4) || (prn >= 11) % loop only for B1 + RF sens
            nps = numel(PPos.(curr_seq_name));
        else
            nps = 1; % use only unfiltered Siemens inputs
        end
        for cc = 1:nps
            if isempty(ls) % check if folder is empty (rf_map!)
                continue
            end
            cd(strcat(cpn,filesep,PName{PPos.(curr_seq_name){cc}}));
            nf = dir('*.nii');
            nfi = fullfile({nf.folder},{nf.name});
            nfi = reshape(nfi,[numel(nfi),1]);
            if contains(lower(curr_seq_name),'map') && (vcd == 's')
                inputs{5,pn} = nfi;
                if numel(nfi) < 2 % files are spread about folders
                    nf = dir(fullfile(PName{PPos.(curr_seq_name){2}},'*.nii'));
                    nfi = fullfile({nf.folder},{nf.name});
                    inputs{5,pn}{2,1} = char(nfi);
                end
            elseif contains(lower(curr_seq_name),'map') && (vcd == 'p')
                for fcb1 = 1:numel(nfi)
                    if contains(get_metadata_val(nfi{fcb1},'PrepulseType'),'NO')
                        if contains(get_metadata_val(nfi{fcb1},'imtype'),'M_B1')
                            inputs{5,pn}{2,1} = char(nfi{fcb1});
                        elseif contains(get_metadata_val(nfi{fcb1},'imtype'),'M_FFE')
                            inputs{5,pn}{1,1} = char(nfi{fcb1});
                        end
                    elseif contains(get_metadata_val(nfi{fcb1},'PrepulseType'),'SAT')
                        if get_metadata_val(nfi{fcb1},'FlipAngle') == 60
                            inputs{5,pn}{2,1} = char(nfi{fcb1});
                        elseif get_metadata_val(nfi{fcb1},'FlipAngle') == 120
                            inputs{5,pn}{1,1} = char(nfi{fcb1});
                        end
                    end
                end
            elseif contains(curr_seq_name,'head') && (vcd == 'p')
                if PPos.(curr_seq_name){cc} == PPos.MTw{1} - 2
                    inputs{2,pn}{1,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.PDw{1} - 2
                    inputs{3,pn}{1,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.T1w{1} - 2
                    inputs{4,pn}{1,1} = char(nfi);
                end
            elseif contains(curr_seq_name,'head') && (vcd == 's')
                if PPos.(curr_seq_name){cc} == PPos.MT{1} - 2
                    inputs{2,pn}{1,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.PD{1} - 2
                    inputs{3,pn}{1,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.T1{1} - 2
                    inputs{4,pn}{1,1} = char(nfi);
                end
            elseif contains(curr_seq_name,'body') && (vcd == 'p')
                if PPos.(curr_seq_name){cc} == PPos.MTw{1} - 1
                    inputs{2,pn}{2,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.PDw{1} - 1
                    inputs{3,pn}{2,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.T1w{1} - 1
                    inputs{4,pn}{2,1} = char(nfi);
                end
            elseif contains(curr_seq_name,'body') && (vcd == 's')
                if PPos.(curr_seq_name){cc} == PPos.MT{1} - 1
                    inputs{2,pn}{2,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.PD{1} - 1
                    inputs{3,pn}{2,1} = char(nfi);
                end
                if PPos.(curr_seq_name){cc} == PPos.T1{1} -1
                    inputs{4,pn}{2,1} = char(nfi);
                end
            elseif contains(curr_seq_name,'MT')
                inputs{6,pn} = nfi;
            elseif contains(curr_seq_name,'PD')
                inputs{7,pn} = nfi;
            elseif contains(curr_seq_name,'T1') && ~contains(upper(curr_seq_name),'T1W_TSE')
                inputs{8,pn} = nfi;
            end
        end
    end
    
    %% write data to files
    prot_file_name = fullfile(cpn,'all_protocol_settings.mat');
    save(prot_file_name,'protocol_settings','-v7.3');
    inputs_file_name = fullfile(cpn,'inputs.mat');
    save(inputs_file_name,'inputs','-v7.3');
    fprintf(fid,['<p>Protocol settings and inputs for processing saved as HDF5 compatible Matlab files:<br>\n'...
        '<a href="%s">all_protocol_settings.mat</a><br>\n'...
        '<a href="%s">inputs.mat</a><br></p>\n'],prot_file_name,inputs_file_name);
    fprintf(fid,'</body></html>');
    fclose(fid); 
end

function tol = tolchk(actual,target,tolerance)
if abs(actual - target) <= tolerance
    tol = true;
else
    tol = false;
end
end


%% function for retrieving and checking the parameters for a distinct
% secquence, distinguishing the parameters to check and give feedback
% about mismatch of parameters; outputting the number of match/mismatch
function [check_okay,check_not_okay] = pcheck(vendor,sequence,mf,ind,tolerance,fid)
% initializing counters for match/mismatch
check_okay = 0;
check_not_okay = 0;

% set the json file with the current settings from the NIFTI file
cfname = spm_file(mf,'number','');
json = strrep(cfname,'.nii','.json');
echo_num = ind;% str2num(json(length(json)-5));

% set the json file with the standard settings and the current sequence
if contains(upper(sequence),'MT') || contains(upper(sequence),'PD') || ...
        contains(upper(sequence),'T1')  && ~contains(upper(sequence),'T1W_TSE')
    if vendor == 'p'
        cprn = sprintf('%s%i',sequence,echo_num);
        prot = protocol_standard.(vendor).(cprn);
    else
        prot = protocol_standard.(vendor).(sequence);
        cprn = sprintf('%s_%i',sequence,ind);
    end
else
    cprn = sequence;
    prot = protocol_standard.(vendor).(sequence);
end

prot = prot{1,1};

if ind % echo_num ~= 0, i.e. for multi-echo data
    fprintf(fid,'echo: %i, sequence: %s<br>\n',echo_num,sequence);
end

%% setting the search strings for get_metadata_val
% distinguish between numerical (num_val) and string values (str_val);
% switching the protocol between iPAT2 and iPAT3 only for MPMs
if echo_num > 1
    num_val = {'EchoTime'}; str_val = {};
else
    switch vendor
        case {'s'}
            num_val = {'EchoTime' 'RepetitionTime' 'FlipAngle' ...
                'Rows' 'Columns' 'dReadoutFOV' 'dPhaseFOV' 'lImagesPerSlab' ...
                'NumberOfAverages' 'dPhaseResolution' ...
                'NumberOfPhaseEncodingSteps' 'SliceThickness' ...
                'PercentSampling' 'PercentPhaseFieldOfView' ...
                'BandwidthPerPixelRO' 'PixelBandwidth' ...
                'AccelFactorPE' 'AccelFactor3D' 'SliceMeasurementDuration' ...
                'PhaseEncodingDirectionPositive' 'MT'};
            str_val = {'ImaCoilString' 'ScanningSequence' 'SequenceVariant' ...
                'SequenceName' 'PhaseEncodingDirection'};
            
        case {'p'}
            num_val = {'EchoTime' 'RepetitionTime' 'FlipAngle' ...
                'Rows' 'Columns' 'NumberofSlicesMR' 'PixelSpacing_x' 'PixelSpacing_y' ...
                'NumberOfAverages' 'NumberOfPhaseEncodingSteps' ...
                'MRAcquisitionPhaseEncodingStepsInPlane' ...
                'MRAcquisitionPhaseEncodingStepsOutOfPlane' 'SliceThickness' ...
                'PercentSampling' 'PercentPhaseFieldOfView' ...
                'BandwidthPerPixelRO' 'PixelBandwidth' 'AcquisitionDuration' ...
                'ParallelReductionFactorInPlane' 'ParallelReductionFactorOutOfPlane'};
            str_val = {'ReceiveCoilName' 'ScanningSequence' 'SequenceVariant' ...
                'PhaseEncodingDirection' 'InPlanePhaseEncodingDirection' 'PulseSequenceName' ...
                'MagnetizationTransfer'};
    end
    switch sequence
        case {'MEDIC' 't2_tse_sag' 't2_tse_tra' 't1_tse_sag'}
            num_val = {'EchoTime' 'RepetitionTime' 'FlipAngle' ...
                'Rows' 'Columns' 'NumberOfAverages' 'dPhaseResolution' ...
                'NumberOfPhaseEncodingSteps' 'SliceThickness' ...
                'PercentSampling' 'PercentPhaseFieldOfView' ...
                'BandwidthPerPixelRO' 'PixelBandwidth' ...
                'AccelFactorPE' 'AccelFactor3D' 'SliceMeasurementDuration' ...
                'PhaseEncodingDirectionPositive' 'MT'}; % 'PELines'='PELinesPAT'; 'PELinesPF'='NumberOfPhaseEncodingSteps'
            str_val = {'ImaCoilString' 'ScanningSequence' 'SequenceVariant' ...
                'SequenceName' 'PhaseEncodingDirection'};
        case 'rf_map'
            num_val = {'SliceThickness' 'RepetitionTime' 'NumberOfAverages' ...
                'NumberOfPhaseEncodingSteps' ...
                'PercentSampling' 'PercentPhaseFieldOfView' 'MT'};
            str_val = {'ImaCoilString' 'ScanningSequence' 'SequenceVariant' ...
                'SequenceName' 'PhaseEncodingDirection'};
        case {'head' 'body'}
            num_val = {'SliceThickness' 'RepetitionTime' 'NumberOfAverages' ...
                'NumberOfPhaseEncodingSteps' ...
                'PercentSampling' 'PercentPhaseFieldOfView' 'MT'};
            str_val = {'ImaCoilString' 'ScanningSequence' 'SequenceVariant' ...
                'SequenceName' 'PhaseEncodingDirection'};
            if vendor == 'p'
                num_val = {'SliceThickness' 'RepetitionTime' 'NumberOfAverages' ...
                    'NumberofSlicesMR' 'NumberOfPhaseEncodingSteps' ...
                    'PercentSampling' 'PercentPhaseFieldOfView' 'MT'};
                str_val = {'ReceiveCoilName' 'ScanningSequence' 'SequenceVariant' ...
                    'PhaseEncodingDirection'}; % 'PulseSequenceName'
            end
        case {'B1Map'}
            num_val = {'SliceThickness' 'RepetitionTime' 'NumberOfAverages' ...
                'NumberOfPhaseEncodingSteps' 'FlipAngle' ...
                'PercentSampling' 'PercentPhaseFieldOfView' 'MT'};
            str_val = {'ReceiveCoilName' 'ScanningSequence' 'SequenceVariant' ...
                'PhaseEncodingDirection'}; % 'PulseSequenceName'
        case {'T1w_TSE' 'T2w_sag' 'T2w_tra' 'mFFE'}
            num_val = {'SliceThickness' 'RepetitionTime' 'NumberOfAverages' ...
                'NumberofSlicesMR' 'NumberOfPhaseEncodingSteps' ...
                'PercentSampling' 'PercentPhaseFieldOfView' 'MT'};
            str_val = {'ReceiveCoilName' 'ScanningSequence' 'SequenceVariant' ...
                'PhaseEncodingDirection'}; % 'PulseSequenceName'
    end
end
%% checking the numerical values
for nvc = 1:numel(num_val)
    cval = num_val{nvc};
    if contains(cval,'dReadoutFOV') || contains(cval,'dPhaseFOV')
        [tv,tv_src] = get_metadata_val(prot,cval);
        if iscell(tv)
            tv = tv{contains(tv_src,'asSlice')};
        end
        [av, av_src] = get_metadata_val(json,cval);
        if iscell(av)
            av = av{contains(av_src,'asSlice')};
        end
    elseif (vendor == 'p')
        tv = get_metadata_val(prot,cval);
        if iscell(tv)
            tv = tv{1};
        end
        av = get_metadata_val(json,cval);
        if iscell(av)
            av = av{1};
        end
    elseif contains(cval,'PixelSpacing')
        tv = get_metadata_val(prot,'PixelSpacing');
        av = get_metadata_val(json,'PixelSpacing');
        if iscell(av)
            av = cell2mat(av);  % avoid problems with e.g. [1,1,]
        end
        if contains(cval,'x')
            tv = tv(1);
            av = av(1);
        else  % y
            tv = tv(2);
            av = av(2);
        end
    else
        tv = get_metadata_val(prot,cval);
        av = get_metadata_val(json,cval);
    end
    if contains(vendor,'s') && contains(cval,'EchoTime') && echo_num
        tv = get_metadata_val(prot,'EchoTimes');
        tv = tv(echo_num);
    end
    protocol_settings.(cprn).num.(cval).actual = av;
    protocol_settings.(cprn).num.(cval).reference = tv;
    if isempty(av)
        fprintf(fid,'<font color="red">Missing parameter ''%s'' in metadata - please check pseudonymization process!</font><br>\n',cval);
        continue
    end
    if ~tolchk(av,tv,tolerance)
        fprintf(fid,['<font color="red">Mismatch of %s:',...
            ' actual value is %5.2f instead of %5.2f (deviation of %3.2f%%).</font><br>\n'],...
            cval, av, tv,100*(av-tv)/tv);
        check_not_okay = check_not_okay + 1;
        protocol_settings.(cprn).num.(cval).match = false;
    else
        check_okay = check_okay + 1;
        fprintf(fid,'Actual value of "%s" (%5.2f): okay.<br>\n',cval,av);
        protocol_settings.(cprn).num.(cval).match = true;
    end
end

%% checking the string values
for svc = 1:numel(str_val)
    cval = str_val{svc};
    av = get_metadata_val(json,cval);
    if iscell(av)
        av = av{1};
    end
    if isempty(av)
        fprintf(fid,'<font color="red">Missing parameter ''%s'' in metadata - please check pseudonymization process!</font><br>\n',cval);
        continue
    end
    protocol_settings.(cprn).str.(cval).actual = av;
    tv = get_metadata_val(prot,cval);
    if contains(cval,'Coil') && (~contains(lower(tv),'b') && ~contains(lower(av),lower(tv)))
        tv = 'SENSE'; % alternative coil name
    end
    protocol_settings.(cprn).str.(cval).reference = tv;
    if ~contains(lower(av),lower(tv))
        fprintf(fid,'<font color="red">Mismatch of %s: actual value is "%s" instead of "%s".</font><br>\n',...
            cval, av, tv);
        check_not_okay = check_not_okay + 1;
        protocol_settings.(cprn).str.(cval).match = false;
    else
        check_okay = check_okay + 1;
        fprintf(fid,'Actual value of "%s" ("%s"): okay.<br>\n',cval,av);
        protocol_settings.(cprn).str.(cval).match = true;
    end
end

end
end
