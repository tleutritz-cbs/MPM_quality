% function to gzip or gunzip NIFTIs in order to save space
% required inputs:
% path - directory, where to gzip/gunzip NIFTIs
% modes - switch between gzip and gunzip: 'gz' or 'nii'
function gunzip_or_gzip(path,mode)
cd(path);

% look for nii and nii.gz files
nfts = dir('*.nii');
nigz = dir('*.nii.gz');

switch mode
    case 'gz'
        if isempty(nigz)
            for ng = 1:numel(nfts)
                gzip(fullfile(nfts(ng,1).folder,nfts(ng,1).name));
                delete(fullfile(nfts(ng,1).folder,nfts(ng,1).name));
            end
        else
            for ng = 1:numel(nfts)
                if exist(fullfile(nfts(ng,1).folder,strcat(nfts(ng,1).name,'.gz')),'file')
                    delete(fullfile(nfts(ng,1).folder,nfts(ng,1).name));
                end
            end
        end
    case 'nii'
        if isempty(nfts)
            for ng = 1:numel(nigz)
                gunzip(fullfile(nigz(ng,1).folder,nigz(ng,1).name));
                delete(fullfile(nigz(ng,1).folder,nigz(ng,1).name));
            end
        else
            for ng = 1:numel(nigz)
                if exist(fullfile(nigz(ng,1).folder,strrep(nigz(ng,1).name,'.gz','')),'file')
                    delete(fullfile(nigz(ng,1).folder,nigz(ng,1).name));
                end
            end
        end
end
end
