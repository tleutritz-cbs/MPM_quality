# Preface & important remarks
The QA.m script requires MATLAB and SPM12 with the [hMRI-toolbox](https://hmri.info/).

Please start SPM before running the script, which takes as input a single DICOM directory (i.e. call it by QA('full-path-to-DICOM-files')), which will then be converted to NiFTI files including acquisition parameter check, quantitative maps will be created and quality scripts on these maps and the original multi-echo data will be run as well.

Visual checks are still necessary for all input data, maps and segmentations!

# Adaptation possibilities & useful tools
Despite the fact that these scripts are dedicated for the [NISCI-trial](https://nisci-2020.eu/), it can be adopted to needs of other protocols run in the context of multi-parameter mapping (MPM) with the hMRI-toolbox.

Some hints will appear in the wiki or ask the author [T. Leutritz](mailto:tleutritz@cbs.mpg.de).

In the progress of conducting these scripts, some usefule scripts were created, which might be of interest also for other use-cases:
* create_head_mask2.m (creates a head from various input data as well as a ghost mask, currently just for saggital acquisition schemes)
* create_seg_from_mat.m (create segmentations and inverse deformation field just from the seg8.mat file from SPM's segmentation)
* gunzip_or_gzip.m (zip/gunzip multiple NiFTI files within a directory including deletion of files/archives afterwards to avoid multiple versions)
* rename_folders.m (renames NiFTI folders created by SPM's DICOM converter in order to have the numbering at the beginning of folder names)
