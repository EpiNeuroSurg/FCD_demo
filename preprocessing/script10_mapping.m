1. Setup your Freesurfer environment
bash
export FREESURFER_HOME=/usr/local/freesurfer_de
source $FREESURFER_HOME/SetUpFreeSurfer.sh
SUBJECTS_DIR=$PWD

2. Realign an SPM-map with Freesurfer’s average MNI scan
mri_convert /usr/local/freesurfer/subjects/fsaverage_sym/mri/orig.mgz /home/mo/Desktop/mniImage.nii

#####post.nii registered to pre_anat.nii
#####pre_anat.nii and rpost.nii normalized (east and 111)

pt='18xiaozhiguo'
side='lh'
# 3. Create a projection of your map onto a surface
mri_vol2surf --mov $SUBJECTS_DIR/$pt/spmT_0001.nii --o $SUBJECTS_DIR/$pt/mapping_TT.nii --hemi $side --regheader fsaverage_sym --projdist 1.5
# 4. convert nii2mgh
mri_convert $SUBJECTS_DIR/$pt/mapping_TT.nii $SUBJECTS_DIR/$pt/mapping_TT.mgh 

pt='18xiaozhiguo'
side='rh'
## surgical mapping
# 3. Create a projection of your map onto a surface
mri_vol2surf --mov $SUBJECTS_DIR/$pt/surg.nii --o $SUBJECTS_DIR/$pt/mapping_surg.nii --hemi $side --regheader fsaverage_sym --projdist 1.5
# 4. convert nii2mgh
mri_convert $SUBJECTS_DIR/$pt/mapping_surg.nii $SUBJECTS_DIR/$pt/mapping_surg.mgh 




