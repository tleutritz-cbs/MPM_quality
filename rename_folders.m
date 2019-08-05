function rename_folders(varargin)
% script to rename folders from SPM DICOM2NIFTI standard structure 
% in form of 'ProtocolName_XXXX' to 'XXXX_ProtocolName'
% adapted from the code of Jan Simon from 2011-09-20 at
% https://de.mathworks.com/matlabcentral/answers/16283-renaming-a-lot-of-folders-automatically-by-matlab
% REQUIREMENT for a call without arguments: running SPM (for open dialog)
% by Tobias Leutritz (2016-08-02)

if nargin == 0
    % spm dialog to gather the folder including the folders to be renamed
    PathName = cfg_getfile(Inf,'dir',...
        'Select folders containing the folders to be renamed');
    [npth,~] = size(PathName);
else
    npth = numel(varargin);
    for n = 1:npth
        PathName(n,:) = varargin{1,n}; %#ok<AGROW>
    end
end


% loop over given directiories
for pn = 1:npth
    cpn = char(PathName(pn,:));

    % get folder directories
    ADir = dir(cpn);

    % make sure to only have folders in the list
    ADir = ADir([ADir.isdir]);

    % convert to string array
    AName = {ADir.name};

    % get rid of folder substitutions '.' and '..'
    AName(strncmp(AName, '.', 1)) = [];

    % loop over the list of folders
    for iFolder = 1:numel(AName)
      AFolder = AName{iFolder};
      APath   = fullfile(cpn, AFolder);
      len = length(AFolder);

      % check, if the folders are already in the correct format
      [x,sorted] = str2num(AFolder(1:4));
      if sorted 
          disp('Folders are already in the desired format.');
          break
      end

      % create new name, assuming the above mentioned format with 4 digits
      newName = sprintf('%s_%s', AFolder(len-3:len), AFolder(1:len-5));
      movefile(APath, fullfile(cpn, newName));
    end
end
end