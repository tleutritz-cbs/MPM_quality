function hmri_local_defaults_NISCI_THS_PH_thrA10e8
% PURPOSE
% To set user-defined (site- or protocol-specific) defaults parameters
% which are used by the hMRI toolbox. Customized processing parameters can
% be defined, overwriting defaults from hmri_defaults. Acquisition
% protocols can be specified here as a fallback solution when no metadata
% are available. Note that the use of metadata is strongly recommended. 
%
% RECOMMENDATIONS
% Parameters defined in this file are identical, initially, to the ones
% defined in hmri_defaults.m. It is recommended, when modifying this file,
% to remove all unchanged entries and save the file with a meaningful name.
% This will help you identifying the appropriate defaults to be used for
% each protocol, and will improve the readability of the file by pointing
% to the modified parameters only.
%
% WARNING
% Modification of the defaults parameters may impair the integrity of the
% toolbox, leading to unexpected behaviour. Only recommended for expert
% users. 
%
% HOW DOES IT WORK?
% The modified defaults file can be selected using the "Configure toolbox"
% branch of the hMRI-Toolbox. For customization of B1 processing
% parameters, type "help hmri_b1_standard_defaults.m". 
%__________________________________________________________________________
% Written by E. Balteau, 2017.
% Cyclotron Research Centre, University of Liege, Belgium

% Global hmri_def variable used across the whole toolbox
global hmri_def

%==========================================================================
% Common processing parameters 
%==========================================================================

% Define the OLS fit as default. OLS fit at TE=0 is used instead of
% averaged contrast images for the map calculation
hmri_def.fullOLS = true;

% recommended TPM for segmentation and spatial processing
hmri_def.TPM = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'etpm','eTPM.nii');

%==========================================================================
% R1/PD/R2s/MT map creation parameters
%==========================================================================
%--------------------------------------------------------------------------
% Threshold values for qMRI maps
%--------------------------------------------------------------------------
hmri_def.qMRI_maps_thresh.MT       = 30;
hmri_def.qMRI_maps_thresh.A        = 10^8; % after rescaling
%hmri_def.qMRI_maps.QA          = 1; % 0 for testing only
hmri_def.cleanup = false;
%hmri_def.segment.tissue(1).warped = [0 1]; % output warped GM
%hmri_def.segment.tissue(2).warped = [0 1]; % output warped WM
hmri_def.segment.warp.write = [1 0]; % output inverse deformation field
end
