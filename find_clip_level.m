%% This script helps to find a clip level for getting rid of background
% noise in e.g. MRI data and is inspired by the AFNI tool 3dClipLevel:
% "Algorithm:
%   (a) Set some initial clip value using wizardry (AKA 'variance').
%   (b) Find the median of all positive values >= clip value.
%   (c) Set the clip value to 0.50 of this median.
%   (d) Loop back to (b) until the clip value doesn't change."
% [https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dClipLevel.html]
% Input: data as a matrix (n-dim), optional mfrac (as in the AFNI version)
% Tobias Leutritz, Neurophysik Leipzig, 04.04.2018

%% the function code starts here:
function clip_level = find_clip_level(data, mfrac)

%% initialization:
if nargin < 2
    mfrac = 0.5;
end
iter = 1;
clip_level = 0;
init(iter) = std(data(:)); % (a) changed from variance to standard deviation

%% the algorithm:
while clip_level ~= init(iter)
    iter = iter + 1;
    med_val = median(data(data>=init(iter-1))); % (b)
    init(iter) = mfrac * med_val; % (c)
    clip_level = init(iter-1);
end

%% the result:
clip_level = init(iter);

end