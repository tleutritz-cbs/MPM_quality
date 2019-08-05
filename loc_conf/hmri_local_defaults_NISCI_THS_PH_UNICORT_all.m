function hmri_local_defaults_NISCI_THS_PH_UNICORT_all
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

hmri_def.cleanup = false;

% recommended TPM for segmentation and spatial processing
hmri_def.TPM = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'etpm','eTPM.nii');

%==========================================================================
% R1/PD/R2s/MT map creation parameters
%==========================================================================
%--------------------------------------------------------------------------
% Threshold values for qMRI maps
%--------------------------------------------------------------------------
hmri_def.qMRI_maps_thresh.MT       = 15;
hmri_def.qMRI_maps_thresh.A        = 10^8; % after rescaling
hmri_def.UNICORT.PD = true;
hmri_def.UNICORT.MT = true;
end
